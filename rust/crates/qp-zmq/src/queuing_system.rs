//! The task queue of the task server (`ocaml/Queuing_system.ml`).
//!
//! The reference keeps two OCaml lists, `queued_front` (newest task first)
//! and `queued_back` (next task first), plus three maps/sets for the running
//! tasks, the task payloads and the connected clients. This port keeps the
//! same shape so that the observable order matches:
//!
//! * [`QueuingSystem::add_task`] pushes on the front list;
//! * [`QueuingSystem::pop_task`] moves the front list, reversed, into the
//!   back list when the back list is empty, then pops the back list, so the
//!   overall order is FIFO;
//! * `del_task` removes the payload from the task map **without** touching
//!   either queue list. A deleted task id can therefore still be popped; it
//!   then has no payload (the reference's `Not_found` branch), and a second
//!   `del_task` fails.
//!
//! # Counts
//!
//! The reference maintains four integer counters next to the collections.
//! They are updated exactly when the collections are, so this port derives
//! them from the collection lengths instead. In particular
//! `number_of_queued` is the total length of the two queue lists and is
//! *not* decremented by `del_task`; `number_of_tasks` is the map length and
//! is.
//!
//! # Failure semantics
//!
//! A failure carries one plain message: the reference's `Failure` payload or
//! its `file:line:column` assertion site. The task server puts that message
//! in its `error <message> : <raw request>` reply; the reply shape is the
//! same for every failure, so the queue does not classify them. The batch
//! methods ([`QueuingSystem::del_tasks`], [`QueuingSystem::end_tasks`])
//! pre-validate every identifier before mutating the queue, reproducing the
//! reference's functional fold: a failure leaves the whole queue untouched.

use std::collections::{HashMap, HashSet, VecDeque};
use std::fmt;

use crate::message::Id;

/// A failure of a queue operation.
///
/// The message is the text the task server prepends to the raw request in
/// its `error` reply (`<message> : <raw>`): the reference's `Failure`
/// payload or its `file:line:column` assertion site. The reply shape is the
/// same either way, so the queue carries no exception class.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct QueueError(String);

impl QueueError {
    /// Wrap a failure message.
    #[must_use]
    pub fn new(message: impl Into<String>) -> Self {
        Self(message.into())
    }
}

impl fmt::Display for QueueError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for QueueError {}

/// Result of one [`QueuingSystem::pop_task`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PoppedTask {
    /// Identifier of the popped task, or `None` when the queue was empty.
    pub task_id: Option<Id>,
    /// Payload of the popped task; `None` when the task was deleted while it
    /// was still queued, or when the queue was empty.
    pub task: Option<String>,
}

/// The task queue, with the reference's two-list / map / set layout.
#[derive(Debug)]
pub struct QueuingSystem {
    /// Oldest task first, so the reference's `queued_front` head (newest) is
    /// the last element here.
    queued_front: Vec<Id>,
    /// Next task to pop at the front.
    queued_back: VecDeque<Id>,
    /// Task id -> client id currently running it.
    running: HashMap<Id, Id>,
    /// Task id -> payload.
    tasks: HashMap<Id, String>,
    /// Connected clients.
    clients: HashSet<Id>,
    /// Identifier the next [`QueuingSystem::add_client`] will return.
    next_client_id: Id,
    /// Identifier the next [`QueuingSystem::add_task`] will return.
    next_task_id: Id,
}

impl Default for QueuingSystem {
    fn default() -> Self {
        Self::new()
    }
}

impl QueuingSystem {
    /// An empty queue with both identifier counters at 1
    /// (`Queuing_system.create`).
    #[must_use]
    pub fn new() -> Self {
        Self {
            queued_front: Vec::new(),
            queued_back: VecDeque::new(),
            running: HashMap::new(),
            tasks: HashMap::new(),
            clients: HashSet::new(),
            next_client_id: Id::new(1),
            next_task_id: Id::new(1),
        }
    }

    /// Append `task` and return its identifier (`Queuing_system.add_task`).
    ///
    /// The identifier is strictly increasing over the life of the queue and
    /// is never reused.
    pub fn add_task(&mut self, task: impl Into<String>) -> Id {
        let task_id = self.next_task_id;
        self.next_task_id = task_id.next();
        self.queued_front.push(task_id);
        self.tasks.insert(task_id, task.into());
        task_id
    }

    /// Register a new client and return its identifier
    /// (`Queuing_system.add_client`).
    pub fn add_client(&mut self) -> Id {
        let client_id = self.next_client_id;
        self.next_client_id = client_id.next();
        self.clients.insert(client_id);
        client_id
    }

    /// The identifier of the next task [`QueuingSystem::pop_task`] would
    /// return, without mutating the queue.
    #[must_use]
    pub fn next_queued_id(&self) -> Option<Id> {
        if let Some(task_id) = self.queued_back.front() {
            Some(*task_id)
        } else {
            // `queued_front` is oldest-first, so its first element is next.
            self.queued_front.first().copied()
        }
    }

    /// Whether `client_id` is currently connected
    /// (`ClientsSet.mem client_id q.clients`).
    ///
    /// The task server uses this to reproduce the up-front
    /// `assert (ClientsSet.mem client_id q.clients)` of
    /// `Queuing_system.pop_task` without consuming a task.
    #[must_use]
    pub fn has_client(&self, client_id: Id) -> bool {
        self.clients.contains(&client_id)
    }

    /// Payload of `task_id`, or `None` when it was deleted while queued or
    /// never existed (`TasksMap.find`).
    ///
    /// The task server's `get_tasks` needs the reference's "pop, see no
    /// payload, put the id back" branch; this accessor lets it detect that
    /// case before popping instead of cloning the whole queue.
    #[must_use]
    pub fn task_payload(&self, task_id: Id) -> Option<&str> {
        self.tasks.get(&task_id).map(String::as_str)
    }

    /// Pop one task for `client_id` (`Queuing_system.pop_task`).
    ///
    /// Asserts that `client_id` is a registered client
    /// (`Queuing_system.ml:68:2`).
    ///
    /// # Errors
    ///
    /// The reference's assert site when the client is unknown.
    pub fn pop_task(&mut self, client_id: Id) -> Result<PoppedTask, QueueError> {
        if !self.clients.contains(&client_id) {
            return Err(QueueError::new("Queuing_system.ml:68:2"));
        }
        if self.queued_back.is_empty() {
            // `List.rev queued_front` becomes the new back list. The
            // reference's `queued_front` is newest-first, so `rev` turns it
            // into oldest-first; this field already stores oldest-first
            // (`add_task` appends), so the drain goes in unreversed.
            self.queued_back.extend(self.queued_front.drain(..));
        }
        match self.queued_back.pop_front() {
            Some(task_id) => {
                self.running.insert(task_id, client_id);
                Ok(PoppedTask {
                    task_id: Some(task_id),
                    task: self.tasks.get(&task_id).cloned(),
                })
            }
            None => Ok(PoppedTask {
                task_id: None,
                task: None,
            }),
        }
    }

    /// Remove a client (`Queuing_system.del_client`).
    ///
    /// # Errors
    ///
    /// Asserts that the client is registered
    /// (`Queuing_system.ml:92:2`).
    pub fn del_client(&mut self, client_id: Id) -> Result<(), QueueError> {
        if !self.clients.remove(&client_id) {
            return Err(QueueError::new("Queuing_system.ml:92:2"));
        }
        Ok(())
    }

    /// Mark one running task as finished (`Queuing_system.end_task`).
    ///
    /// # Errors
    ///
    /// Asserts the client is registered (`Queuing_system.ml:103:2`),
    /// raises `Failure "Task already finished"` when the task is not running,
    /// and asserts that the task belongs to this client
    /// (`Queuing_system.ml:109:4`).
    pub fn end_task(&mut self, task_id: Id, client_id: Id) -> Result<(), QueueError> {
        if !self.clients.contains(&client_id) {
            return Err(QueueError::new("Queuing_system.ml:103:2"));
        }
        match self.running.get(&task_id) {
            None => Err(QueueError::new("Task already finished")),
            Some(&owner) if owner != client_id => Err(QueueError::new("Queuing_system.ml:109:4")),
            Some(_) => {
                self.running.remove(&task_id);
                Ok(())
            }
        }
    }

    /// Delete one task payload (`Queuing_system.del_task`).
    ///
    /// The task identifier stays in the queue lists (the reference does not
    /// remove it), and `number_of_queued` is unchanged.
    ///
    /// # Errors
    ///
    /// `Failure "Task <id> is already deleted"` when the payload is absent,
    /// including a repeated identifier.
    pub fn del_task(&mut self, task_id: Id) -> Result<(), QueueError> {
        self.del_tasks(std::slice::from_ref(&task_id))
    }

    /// Delete a batch of task payloads, as the reference's `List.fold_left`
    /// does.
    ///
    /// All-or-nothing, like the reference: `TaskServer.del_task` folds the
    /// *pure* `Queuing_system.del_task` over the identifier list
    /// (`ocaml/TaskServer.ml:376-378`), so a failure at the k-th identifier
    /// escapes before any part of the fold is committed. The mutable port
    /// reproduces that by validating every identifier before removing any.
    /// A repeated identifier fails on its second occurrence: the fold's
    /// first occurrence would have removed it already.
    ///
    /// # Errors
    ///
    /// `Task <id> is already deleted` for the first identifier that is
    /// absent, or that appears twice.
    pub fn del_tasks(&mut self, task_ids: &[Id]) -> Result<(), QueueError> {
        let mut seen = HashSet::with_capacity(task_ids.len());
        for &task_id in task_ids {
            if !self.tasks.contains_key(&task_id) || !seen.insert(task_id) {
                return Err(QueueError::new(format!(
                    "Task {} is already deleted",
                    task_id.get()
                )));
            }
        }
        for &task_id in task_ids {
            self.tasks.remove(&task_id);
        }
        Ok(())
    }

    /// Mark a batch of running tasks as finished, as the reference's
    /// `List.fold_left` does.
    ///
    /// All-or-nothing, for the same reason as [`QueuingSystem::del_tasks`]:
    /// the reference's `TaskServer.task_done` folds the pure
    /// `Queuing_system.end_task` (`ocaml/TaskServer.ml:587-590`), so the
    /// identifiers are validated here before any mutation.
    ///
    /// # Errors
    ///
    /// The same errors as [`QueuingSystem::end_task`], in the fold's order;
    /// a repeated identifier fails on its second occurrence with
    /// `Task already finished`.
    pub fn end_tasks(&mut self, task_ids: &[Id], client_id: Id) -> Result<(), QueueError> {
        if !self.clients.contains(&client_id) {
            return Err(QueueError::new("Queuing_system.ml:103:2"));
        }
        let mut seen = HashSet::with_capacity(task_ids.len());
        for &task_id in task_ids {
            match self.running.get(&task_id) {
                None => return Err(QueueError::new("Task already finished")),
                Some(&owner) if owner != client_id => {
                    return Err(QueueError::new("Queuing_system.ml:109:4"));
                }
                Some(_) if !seen.insert(task_id) => {
                    // The fold would have removed it on the first occurrence.
                    return Err(QueueError::new("Task already finished"));
                }
                Some(_) => {}
            }
        }
        for &task_id in task_ids {
            self.running.remove(&task_id);
        }
        Ok(())
    }

    /// Number of task payloads still stored (`Queuing_system.number_of_tasks`).
    #[must_use]
    pub fn number_of_tasks(&self) -> usize {
        self.tasks.len()
    }

    /// Number of ids in the two queue lists
    /// (`Queuing_system.number_of_queued`).
    ///
    /// Deleted payloads still count, as in the reference.
    #[must_use]
    pub fn number_of_queued(&self) -> usize {
        self.queued_front.len() + self.queued_back.len()
    }

    /// Number of tasks currently running (`Queuing_system.number_of_running`).
    #[must_use]
    pub fn number_of_running(&self) -> usize {
        self.running.len()
    }

    /// Number of connected clients (`Queuing_system.number_of_clients`).
    #[must_use]
    pub fn number_of_clients(&self) -> usize {
        self.clients.len()
    }

    /// The reference's debug dump (`Queuing_system.to_string`).
    ///
    /// Task identifiers are printed sorted, as `Map.bindings` does.
    #[must_use]
    pub fn dump(&self) -> String {
        let queued = self
            .queued_front
            .iter()
            .rev()
            .chain(self.queued_back.iter().rev())
            .map(Id::to_string)
            .collect::<Vec<_>>()
            .join(" ; ");

        let mut running = self.running.iter().collect::<Vec<_>>();
        running.sort_unstable_by_key(|(task_id, _)| **task_id);
        let running = running
            .into_iter()
            .map(|(task_id, client_id)| format!("({task_id}, {client_id})"))
            .collect::<Vec<_>>()
            .join(" ; ");

        let mut tasks = self.tasks.iter().collect::<Vec<_>>();
        tasks.sort_unstable_by_key(|(task_id, _)| **task_id);
        let tasks = tasks
            .into_iter()
            .map(|(task_id, task)| format!("({task_id}, \"{task}\")"))
            .collect::<Vec<_>>()
            .join(" ; ");

        format!(
            "{{\nTasks : {}   Queued : {}   Running : {}   Clients : {}\n\
             queued   : {{ {} }}\n\
             running  : {{ {} }}\n\
             tasks    : [ {}\n           ]\n}}",
            self.number_of_tasks(),
            self.number_of_queued(),
            self.number_of_running(),
            self.number_of_clients(),
            queued,
            running,
            tasks
        )
    }
}

impl fmt::Display for QueuingSystem {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.dump())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Golden output of `test_queuing_system.ml` (Core), produced by
    /// compiling the reference `Queuing_system.ml` + `Id.ml` with the same
    /// scenario and capturing `to_string`. See `implementation.md`.
    #[test]
    fn fifo_order_matches_the_reference_test() {
        let mut q = QueuingSystem::new();
        for i in 0..3 {
            q.add_task(format!("Task {i}"));
        }
        let client = q.add_client();
        let first = q.pop_task(client).unwrap();
        assert_eq!(first.task_id, Some(Id::new(1)));
        assert_eq!(first.task.as_deref(), Some("Task 0"));
        let second = q.pop_task(client).unwrap();
        assert_eq!(second.task_id, Some(Id::new(2)));
        assert_eq!(second.task.as_deref(), Some("Task 1"));
        let third = q.pop_task(client).unwrap();
        assert_eq!(third.task_id, Some(Id::new(3)));
        assert_eq!(third.task.as_deref(), Some("Task 2"));
        assert_eq!(
            q.pop_task(client).unwrap(),
            PoppedTask {
                task_id: None,
                task: None
            }
        );
    }

    #[test]
    fn dump_matches_the_reference_layout() {
        let mut q = QueuingSystem::new();
        q.add_task("First Task");
        q.add_task("Second Task");
        assert_eq!(
            q.dump(),
            "{\n\
             Tasks : 2   Queued : 2   Running : 0   Clients : 0\n\
             queued   : { 2 ; 1 }\n\
             running  : {  }\n\
             tasks    : [ (1, \"First Task\") ; (2, \"Second Task\")\n           ]\n\
             }"
        );
        let client = q.add_client();
        q.pop_task(client).unwrap();
        assert_eq!(
            q.dump(),
            "{\n\
             Tasks : 2   Queued : 1   Running : 1   Clients : 1\n\
             queued   : { 2 }\n\
             running  : { (1, 1) }\n\
             tasks    : [ (1, \"First Task\") ; (2, \"Second Task\")\n           ]\n\
             }"
        );
    }

    #[test]
    fn deleted_payload_stays_queued_and_pops_without_content() {
        let mut q = QueuingSystem::new();
        let task_id = q.add_task("gone");
        q.del_task(task_id).unwrap();
        assert_eq!(q.number_of_queued(), 1);
        assert_eq!(q.number_of_tasks(), 0);
        let client = q.add_client();
        let popped = q.pop_task(client).unwrap();
        // Reference: `Some task_id, None`; the server replies `get_task_reply 0`
        // but the id is moved to the running map.
        assert_eq!(popped.task_id, Some(task_id));
        assert_eq!(popped.task, None);
        assert_eq!(q.number_of_running(), 1);
        assert_eq!(
            q.del_task(task_id),
            Err(QueueError::new("Task 1 is already deleted"))
        );
    }

    #[test]
    fn repeated_delete_fails_on_the_second_occurrence() {
        let mut q = QueuingSystem::new();
        let task_id = q.add_task("x");
        assert_eq!(
            q.del_tasks(&[task_id, task_id]),
            Err(QueueError::new("Task 1 is already deleted"))
        );
        // Nothing was deleted: the fold's first step is not committed.
        assert_eq!(q.number_of_tasks(), 1);
    }

    #[test]
    fn end_task_checks_owner_and_membership() {
        let mut q = QueuingSystem::new();
        let task_id = q.add_task("x");
        let client = q.add_client();
        let other = q.add_client();
        q.pop_task(client).unwrap();
        assert_eq!(
            q.end_task(task_id, other),
            Err(QueueError::new("Queuing_system.ml:109:4"))
        );
        assert_eq!(
            q.end_tasks(&[task_id, task_id], client),
            Err(QueueError::new("Task already finished"))
        );
        assert_eq!(q.number_of_running(), 1);
        q.end_task(task_id, client).unwrap();
        assert_eq!(
            q.end_task(task_id, client),
            Err(QueueError::new("Task already finished"))
        );
        assert_eq!(
            q.end_task(task_id, Id::new(99)),
            Err(QueueError::new("Queuing_system.ml:103:2"))
        );
    }

    #[test]
    fn pop_requires_a_registered_client() {
        let mut q = QueuingSystem::new();
        assert_eq!(
            q.pop_task(Id::new(1)),
            Err(QueueError::new("Queuing_system.ml:68:2"))
        );
    }
}

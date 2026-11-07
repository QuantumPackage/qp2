//! `-g` in address mode, end to end.
//!
//! A fake remote REP on `P+9` answers the input protocol: `test` → `OK`,
//! `get_ezfio_filename` → `testdir`, `get_input` → the captured
//! `no_work.tar.gz` fixture. The child is a real process, so the whole
//! reference sequence is checked byte for byte:
//!
//! ```text
//! Connect to:
//! tcp://127.0.0.1:<lp>
//! Communication [ OK ]
//! Getting input... testdir ...done
//! <ssh block>
//! Ready
//! qp_tunnel exited properly.
//! ```
//!
//! The fixture has no `testdir/work/`: the reference fails there with
//! `Unix.EEXIST "mkdir"` because `Qpackage.ezfio_work` calls `Ezfio.set_file`
//! on a plain directory. This port creates `work/` and writes the address
//! (documented deviation), which is what the assertions below check.

mod support;

use std::fs;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::thread;
use std::time::Duration;

use qp_tunnel::ssh_instructions;
use support::{Running, free_remote_window, parse_local_port};
use zmq::{Context, SocketType};

/// The measured remote tarball, copied into `tests/data/`.
const TARBALL: &[u8] = include_bytes!("data/no_work.tar.gz");

/// Its one file, read out of the fixture with `tar -xzOf` when it was copied.
const INPUT_TXT: &str = "hello from remote\n";

/// A fake remote input peer on `port`.
struct FakeInput {
    stop: Arc<AtomicBool>,
    handle: Option<thread::JoinHandle<()>>,
}

impl FakeInput {
    fn start(port: u32) -> Self {
        let context = Context::new();
        let socket = context.socket(SocketType::REP).expect("REP");
        socket.set_linger(0).expect("linger");
        socket
            .bind(&format!("tcp://127.0.0.1:{port}"))
            .expect("bind the fake input peer");
        let stop = Arc::new(AtomicBool::new(false));
        let flag = Arc::clone(&stop);
        let handle = thread::spawn(move || {
            while !flag.load(Ordering::SeqCst) {
                match socket.poll(zmq::POLLIN, 100) {
                    Ok(0) => continue,
                    Ok(_) => {}
                    Err(zmq::Error::EINTR) => continue,
                    Err(error) => panic!("fake input poll failed: {error}"),
                }
                let request = socket.recv_bytes(0).expect("fake input recv");
                let reply: &[u8] = match request.as_slice() {
                    b"test" => b"OK",
                    b"get_ezfio_filename" => b"testdir",
                    b"get_input" => TARBALL,
                    other => panic!("unexpected request: {other:?}"),
                };
                socket.send(reply, 0).expect("fake input send");
            }
        });
        Self {
            stop,
            handle: Some(handle),
        }
    }

    fn stop(mut self) {
        self.stop.store(true, Ordering::SeqCst);
        if let Some(handle) = self.handle.take() {
            handle.join().expect("the fake input thread does not panic");
        }
    }
}

impl Drop for FakeInput {
    fn drop(&mut self) {
        self.stop.store(true, Ordering::SeqCst);
        if let Some(handle) = self.handle.take() {
            let _ = handle.join();
        }
    }
}

#[test]
fn get_input_downloads_unpacks_and_writes_the_address() {
    let directory = tempfile::tempdir().expect("tempdir");
    let remote = free_remote_window();
    let fake = FakeInput::start(remote + 9);

    let address = format!("tcp://127.0.0.1:{remote}");
    let mut child = Running::spawn(directory.path(), &["-g", &address], Some("/nonexistent"));
    child.wait_ready();

    let local = parse_local_port(&child.stdout());

    // The reference's exact -g block: no newline between the message and the
    // directory name, one space before `...done`.
    assert!(
        child.stdout().contains("Communication [ OK ]\n"),
        "stdout:\n{}",
        child.stdout()
    );
    assert!(
        child
            .stdout()
            .contains("Getting input... testdir ...done\n"),
        "stdout:\n{}",
        child.stdout()
    );

    // The tarball was unpacked in the child's working directory.
    let unpacked = directory.path().join("testdir");
    assert_eq!(
        fs::read_to_string(unpacked.join("input.txt")).expect("unpacked input.txt"),
        INPUT_TXT
    );
    // Plain directory: `work/` is created, but no EZFIO marker is invented.
    assert!(unpacked.join("work").is_dir());
    assert!(!unpacked.join("work").join("empty").exists());
    // The written address is the announced local window.
    assert_eq!(
        fs::read_to_string(unpacked.join("work").join("qp_run_address")).expect("address"),
        format!("tcp://127.0.0.1:{local}\n")
    );

    // The whole stdout before the signal, against the reference layout.
    let hostname = qp_cli::ffi::hostname().expect("hostname");
    let ssh = ssh_instructions(&hostname, "127.0.0.1", i64::from(remote), local);
    // `ssh_instructions` starts with a blank line and ends after the second
    // ` &`, so `Ready` follows it directly, as in the measured capture.
    let expected_running = format!(
        "Connect to:\ntcp://127.0.0.1:{local}\nCommunication [ OK ]\nGetting input... testdir ...done\n{ssh}Ready\n"
    );
    assert_eq!(child.stdout(), expected_running);

    stop_and_check(&mut child, expected_running);
    fake.stop();
}

/// Give `install_signal_handlers` time to run, send SIGUSR1, check the exit.
///
/// The 300 ms sleep closes the documented post-spawn/pre-install window: the
/// child prints `Ready` from its input thread, which may run before `main`
/// reaches `install_signal_handlers`.
fn stop_and_check(child: &mut Running, expected_running: String) {
    thread::sleep(Duration::from_millis(300));
    child.signal("USR1");
    let status = child.wait();
    assert_eq!(status.code(), Some(0), "stderr:\n{}", child.stderr());
    assert_eq!(
        child.stdout(),
        format!("{expected_running}qp_tunnel exited properly.\n")
    );
    assert_eq!(child.stderr(), "");
}

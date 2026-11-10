//! Fixtures shared by the `qp_run` integration tests.
//!
//! Every test builds a throw-away QP2 tree:
//!
//! ```text
//! <tmp>/root/data/executables   # generated table, $QP_ROOT already resolved
//! <tmp>/root/bin/<name>         # executable /bin/sh shims
//! <tmp>/ezfio/                  # an EZFIO tree (.version only)
//! <tmp>/logs/<n>.{out,err}      # captured output of one child
//! ```
//!
//! and runs the binary Cargo built (`CARGO_BIN_EXE_qp_run`) with `QP_ROOT`
//! pointing at it. The environment is built from scratch — every variable the
//! tool reads is set or removed explicitly — so the developer's shell cannot
//! influence a test, and `QP_NIC` is pinned to the loopback interface so the
//! advertised address is always `tcp://127.0.0.1:<port>`.
//!
//! # Why the output goes to files, not pipes
//!
//! `qp_run` runs its job through `/bin/sh -c` and the job inherits the
//! standard streams. A `qp_run` that hangs would leave the job holding a pipe
//! open, so a reader thread could wait forever even after `kill`. Files have
//! no such problem: on timeout the test kills the launcher, reaps it, and
//! reads whatever the two files contain.
#![allow(dead_code)] // each test binary uses a different subset of the helpers

use std::ffi::OsString;
use std::fs::{File, OpenOptions};
use std::io::Write;
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::process::{Child, Command, ExitStatus, Stdio};
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{Duration, Instant};

/// The binary under test, as built by Cargo for this test run.
pub const BIN: &str = env!("CARGO_BIN_EXE_qp_run");

/// Default deadline for a complete `qp_run`.
///
/// A successful run sleeps one second on purpose (`FINAL_SLEEP`); everything
/// else in a fixture run is milliseconds. The deadline exists so a regression
/// is a failed assertion after seconds, never a hung test.
pub const TIMEOUT: Duration = Duration::from_secs(30);

/// Name of the loopback interface for `QP_NIC`.
///
/// Pinning it makes `TaskServer.ip_address` return `127.0.0.1` without a
/// resolver call, so the advertised address is deterministic and the tests do
/// not depend on the machine's host name resolving.
pub fn loopback_interface() -> &'static str {
    if cfg!(target_os = "macos") {
        "lo0"
    } else {
        "lo"
    }
}

/// Counter that keeps the captured-output file names unique in one process.
static LOG_NUMBER: AtomicU64 = AtomicU64::new(0);

/// A temporary QP2 root, an EZFIO tree and a place for captured output.
pub struct Fixture {
    /// Owns every path below; removed when the test ends.
    dir: tempfile::TempDir,
    /// `$QP_ROOT`: contains `data/executables` and `bin/`.
    root: PathBuf,
    /// EZFIO directory given to `qp_run` as its second argument.
    ezfio: PathBuf,
}

impl Fixture {
    /// Create the tree, the two `qp_edit` shims (`exit 0` / `exit 1`) and an
    /// EZFIO directory.
    ///
    /// Panics when the temporary path contains a space: the
    /// `data/executables` format and the reference's job command both split on
    /// spaces, so such a path would test the parser, not `qp_run`.
    pub fn new() -> Self {
        let dir = tempfile::tempdir().expect("tempdir");
        assert!(
            !dir.path().to_string_lossy().contains(' '),
            "the fixture path must not contain a space: {:?}",
            dir.path()
        );
        let root = dir.path().join("root");
        let ezfio = dir.path().join("ezfio");
        std::fs::create_dir_all(root.join("data")).expect("create data/");
        std::fs::create_dir_all(root.join("bin")).expect("create bin/");
        std::fs::write(root.join("data").join("executables"), "").expect("create data/executables");
        let fixture = Self { dir, root, ezfio };
        fixture.shim("qp_edit_ok", "exit 0");
        fixture.shim("qp_edit_fail", "exit 1");
        drop(qp_ezfio::EzFio::create(&fixture.ezfio).expect("create the EZFIO tree"));
        fixture
    }

    /// The temporary directory that owns everything.
    pub fn dir(&self) -> &Path {
        self.dir.path()
    }

    /// `$QP_ROOT` of this fixture.
    pub fn root(&self) -> &Path {
        &self.root
    }

    /// The EZFIO directory passed to `qp_run`.
    pub fn ezfio(&self) -> &Path {
        &self.ezfio
    }

    /// [`Self::ezfio`] as the `&str` the command line takes.
    pub fn ezfio_str(&self) -> &str {
        self.ezfio.to_str().expect("UTF-8 temporary path")
    }

    /// `<ezfio>/work`, where `qp_run` and the fixtures' shims drop their
    /// markers.
    pub fn work(&self) -> PathBuf {
        self.ezfio.join("work")
    }

    /// Directory that receives the captured output files.
    pub fn logs(&self) -> PathBuf {
        self.dir.path().join("logs")
    }

    /// Write an executable `#!/bin/sh` script `name` in `<root>/bin` and
    /// return its path.
    pub fn shim(&self, name: &str, body: &str) -> PathBuf {
        let path = self.root.join("bin").join(name);
        std::fs::write(&path, format!("#!/bin/sh\n{body}\n")).expect("write shim");
        let mut permissions = std::fs::metadata(&path)
            .expect("shim metadata")
            .permissions();
        permissions.set_mode(0o755);
        std::fs::set_permissions(&path, permissions).expect("chmod shim");
        path
    }

    /// Add `name` to `data/executables`, backed by a shim with `body`, and
    /// return the shim's path.
    ///
    /// This is `Qpackage.executables` + `Programs::find`: the table's second
    /// field is an absolute path with no space, exactly what the build
    /// generates.
    pub fn program(&self, name: &str, body: &str) -> PathBuf {
        let path = self.shim(name, body);
        let text = path.to_str().expect("UTF-8 temporary path");
        assert!(!text.contains(' '), "the table splits on spaces: {text}");
        let mut table = OpenOptions::new()
            .append(true)
            .open(self.root.join("data").join("executables"))
            .expect("open data/executables");
        writeln!(table, "{name} {text}").expect("append to data/executables");
        path
    }

    /// A `Command` for the binary with this fixture's environment and no
    /// standard input.
    pub fn command(&self, args: &[&str]) -> Command {
        let mut command = Command::new(BIN);
        command.args(args).stdin(Stdio::null());
        self.apply_env(&mut command);
        command
    }

    /// Like [`Self::command`], but the child reads one line of standard input
    /// and then `exec`s the binary.
    ///
    /// `exec` keeps the pid, so [`Running::id`] is the pid the port window is
    /// computed from; the test can therefore bind a candidate port *before*
    /// the binary computes its window (issue #393 tests). The child cannot
    /// reach the binary before [`Running::release_gate`] is called.
    pub fn gated_command(&self, args: &[&str]) -> Command {
        let mut command = Command::new("/bin/sh");
        command
            .arg("-c")
            .arg("read x || exit 1; exec \"$@\"")
            .arg("sh")
            .arg(BIN)
            .args(args)
            .stdin(Stdio::piped());
        self.apply_env(&mut command);
        command
    }

    /// Run the binary once and wait up to [`TIMEOUT`].
    pub fn run(&self, args: &[&str]) -> Outcome {
        Running::spawn(self.command(args), &self.logs()).wait(TIMEOUT)
    }

    /// Set or remove every environment variable `qp_run` reads, and prepend
    /// `<root>/bin` to `PATH` so a shim can stand in for `qp_stop`.
    fn apply_env(&self, command: &mut Command) {
        let path = match std::env::var_os("PATH") {
            Some(path) => {
                let mut entries = vec![self.root.join("bin")];
                entries.extend(std::env::split_paths(&path));
                std::env::join_paths(entries).expect("join PATH")
            }
            None => OsString::from(self.root.join("bin")),
        };
        command
            .env("QP_ROOT", &self.root)
            .env("QP_NIC", loopback_interface())
            .env("QP_GIT_SHA1", "0000000")
            .env("QP_GIT_DATE", "1970-01-01")
            .env("QP_GIT_MESSAGE", "integration test")
            .env("QP_RUN_QP_EDIT", self.root.join("bin").join("qp_edit_ok"))
            .env("PATH", path)
            .env_remove("QP_RUN_ADDRESS")
            .env_remove("QP_RUN_ADDRESS_MASTER")
            .env_remove("QP_RUN_MAX_SCANS")
            .env_remove("QP_TASK_DEBUG");
    }
}

/// One running `qp_run`, with its output redirected to files.
pub struct Running {
    /// The launcher (or the gating shell, which becomes the launcher).
    child: Child,
    /// File receiving standard output.
    stdout: PathBuf,
    /// File receiving standard error.
    stderr: PathBuf,
}

impl Running {
    /// Spawn `command`, redirecting its standard output and error to
    /// `logs/<n>.out` and `logs/<n>.err`.
    ///
    /// The caller keeps ownership of standard input: [`Self::release_gate`]
    /// uses it for [`Fixture::gated_command`].
    pub fn spawn(mut command: Command, logs: &Path) -> Self {
        std::fs::create_dir_all(logs).expect("create logs/");
        let number = LOG_NUMBER.fetch_add(1, Ordering::Relaxed);
        let stdout = logs.join(format!("{number}.out"));
        let stderr = logs.join(format!("{number}.err"));
        let out = File::create(&stdout).expect("create the stdout file");
        let err = File::create(&stderr).expect("create the stderr file");
        command.stdout(Stdio::from(out)).stderr(Stdio::from(err));
        let child = command.spawn().expect("spawn qp_run");
        Self {
            child,
            stdout,
            stderr,
        }
    }

    /// Pid of the launcher, also after the gating shell `exec`s it.
    pub fn id(&self) -> u32 {
        self.child.id()
    }

    /// Let a [`Fixture::gated_command`] proceed: send one line on its standard
    /// input and close the pipe.
    ///
    /// # Panics
    ///
    /// When the command was not built by [`Fixture::gated_command`].
    pub fn release_gate(&mut self) {
        let mut stdin = self
            .child
            .stdin
            .take()
            .expect("release_gate needs a gated_command");
        stdin.write_all(b"\n").expect("write the gate line");
        // Dropping `stdin` closes the pipe, so the exec'd binary sees EOF.
    }

    /// Wait for the process to exit, at most `timeout`.
    ///
    /// # Panics
    ///
    /// When the deadline expires: the process is killed and reaped first, so
    /// the test fails with the captured output instead of hanging.
    pub fn wait(mut self, timeout: Duration) -> Outcome {
        let deadline = Instant::now() + timeout;
        loop {
            match self.child.try_wait() {
                Ok(Some(status)) => return self.collect(status),
                Ok(None) if Instant::now() < deadline => {
                    std::thread::sleep(Duration::from_millis(10));
                }
                Ok(None) => {
                    let _ = self.child.kill();
                    let status = self.child.wait().expect("reap the killed qp_run");
                    let outcome = self.collect(status);
                    panic!("qp_run did not exit within {timeout:?}\n{outcome}");
                }
                Err(error) => {
                    let _ = self.child.kill();
                    panic!("try_wait on qp_run failed: {error}");
                }
            }
        }
    }

    /// Kill the process, reap it, and return what it wrote so far.
    ///
    /// Used by tests that must unblock their fixture before failing.
    pub fn kill(mut self) -> Outcome {
        let _ = self.child.kill();
        let status = self.child.wait().expect("reap the killed qp_run");
        self.collect(status)
    }

    /// Read both output files and pair them with `status`.
    fn collect(self, status: ExitStatus) -> Outcome {
        Outcome {
            status,
            stdout: read_or_empty(&self.stdout),
            stderr: read_or_empty(&self.stderr),
        }
    }
}

/// Captured result of one `qp_run` process.
pub struct Outcome {
    /// How the process ended.
    pub status: ExitStatus,
    /// Standard output.
    pub stdout: String,
    /// Standard error.
    pub stderr: String,
}

impl Outcome {
    /// The exit code, or `-1` when the process was killed by a signal.
    pub fn code(&self) -> i32 {
        self.status.code().unwrap_or(-1)
    }

    /// Whether the process exited with status 0.
    pub fn success(&self) -> bool {
        self.status.success()
    }
}

impl std::fmt::Display for Outcome {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "status: {}\n--- stdout ---\n{}--- stderr ---\n{}",
            self.status, self.stdout, self.stderr
        )
    }
}

/// Poll for `path` to exist, up to `timeout`; returns whether it appeared.
pub fn wait_for_file(path: &Path, timeout: Duration) -> bool {
    let deadline = Instant::now() + timeout;
    while Instant::now() < deadline {
        if path.exists() {
            return true;
        }
        std::thread::sleep(Duration::from_millis(10));
    }
    path.exists()
}

/// The `key=value` lines a fixture shim wrote to `<work>/observed`.
///
/// Duplicate keys keep the last value, as a shell would.
pub fn observed(work: &Path) -> std::collections::BTreeMap<String, String> {
    let text = std::fs::read_to_string(work.join("observed")).expect("read work/observed");
    text.lines()
        .filter_map(|line| line.split_once('='))
        .map(|(key, value)| (key.to_string(), value.to_string()))
        .collect()
}

/// Shim body: record what the job inherits and sees, then run `extra`.
///
/// The keys are:
///
/// * `arg1` — the job's first argument, which must be the EZFIO directory;
/// * `address` — `QP_RUN_ADDRESS` as the job sees it in its own environment;
/// * `master` — `QP_RUN_ADDRESS_MASTER`, empty in master mode;
/// * `file` — the content of `<work>/qp_run_address` *while* the job runs, so
///   the test can compare it with the advertised address before `qp_run`
///   removes it.
pub fn observe_script(extra: &str) -> String {
    format!(
        "\
dir=\"$1\"
{{
printf 'arg1=%s\\n' \"$1\"
printf 'address=%s\\n' \"${{QP_RUN_ADDRESS-}}\"
printf 'master=%s\\n' \"${{QP_RUN_ADDRESS_MASTER-}}\"
if [ -f \"$dir/work/qp_run_address\" ]; then
  printf 'file=%s\\n' \"$(cat \"$dir/work/qp_run_address\")\"
fi
}} > \"$dir/work/observed\"
{extra}"
    )
}

/// Read a captured file, substituting a description when it is missing.
fn read_or_empty(path: &Path) -> String {
    std::fs::read_to_string(path)
        .unwrap_or_else(|error| format!("<unreadable {}: {error}>\n", path.display()))
}

//! Fixtures shared by the `qp_tunnel` integration tests.
//!
//! Every test owns a throw-away directory (std [`tempfile::TempDir`]) that
//! holds the child's stdout/stderr, its `TMPDIR` (so the `qp_tunnel*.tar.gz`
//! scratch file never leaves the test) and, for the `-g` test, the unpacked
//! EZFIO directory.
//!
//! The child environment is built from scratch: no variable of the developer's
//! shell leaks in, `QP_NIC` pins the advertised address to `127.0.0.1` (loopback
//! interface name from [`loopback_interface`]), `argv[0]` is exactly
//! `qp_tunnel` (the help text contains it) and `QP_ROOT` is set per test.
//!
//! Output goes to files, not pipes: a hung tunnel cannot then leave a pipe
//! writer alive after the test kills it, and the partial output stays readable
//! for the failure message.
//!
//! # Ports
//!
//! [`free_remote_window`] probes the four ports a fake remote binds (`P`,
//! `P+1`, `P+2`, `P+9`) on `127.0.0.1` with std `TcpListener`s; the window is
//! released before it is returned, so a concurrent process can still steal it
//! (the same race the reference has). The local window is discovered by the
//! binary itself and read back from its tunnel-address announcement.
#![allow(dead_code)] // each test binary uses a different subset of the helpers

use std::fs;
use std::net::TcpListener;
use std::os::unix::process::CommandExt;
use std::path::{Path, PathBuf};
use std::process::{Child, Command, ExitStatus, Stdio};
use std::thread;
use std::time::{Duration, Instant};

/// The binary under test, as built by Cargo for this test run.
pub const BIN: &str = env!("CARGO_BIN_EXE_qp_tunnel");

/// Deadline for a child that is expected to run or fail quickly.
pub const TIMEOUT: Duration = Duration::from_secs(20);

/// Name of the loopback interface for `QP_NIC` (`lo0` on macOS, `lo` on
/// Linux), so the advertised address is `127.0.0.1` on every machine.
pub fn loopback_interface() -> &'static str {
    if cfg!(target_os = "macos") {
        "lo0"
    } else {
        "lo"
    }
}

/// A child `qp_tunnel` with its output files and its private scratch tree.
pub struct Running {
    /// `None` once the child has been reaped.
    child: Option<Child>,
    /// Test-owned directory: output files and `TMPDIR` live here.
    dir: PathBuf,
}

impl Running {
    /// Spawn the binary in a fresh output directory.
    ///
    /// `qp_root` is `None` to leave `QP_ROOT` unset (the guard test) and
    /// `Some(value)` otherwise, including `Some("")` for the empty value.
    pub fn spawn(dir: &Path, arguments: &[&str], qp_root: Option<&str>) -> Self {
        Self::spawn_with_arg0(dir, arguments, qp_root, "qp_tunnel")
    }

    /// [`Self::spawn`] with an explicit `argv[0]` (the help block re-wraps
    /// the usage line with it, so the re-wrap tests invoke the binary as a
    /// renamed reference copy was).
    pub fn spawn_with_arg0(
        dir: &Path,
        arguments: &[&str],
        qp_root: Option<&str>,
        arg0: &str,
    ) -> Self {
        let tmp = dir.join("tmp");
        fs::create_dir_all(&tmp).expect("create the child TMPDIR");
        let stdout = fs::File::create(dir.join("stdout")).expect("create stdout");
        let stderr = fs::File::create(dir.join("stderr")).expect("create stderr");

        let mut command = Command::new(BIN);
        command
            .arg0(arg0)
            .args(arguments)
            .current_dir(dir)
            .env_clear()
            // `-g` shells out to `tar`, which the child finds through PATH.
            .env("PATH", "/usr/bin:/bin:/usr/sbin:/sbin")
            .env("TMPDIR", &tmp)
            .env("QP_NIC", loopback_interface())
            .stdin(Stdio::null())
            .stdout(Stdio::from(stdout))
            .stderr(Stdio::from(stderr));
        if let Some(root) = qp_root {
            command.env("QP_ROOT", root);
        }
        Self {
            child: Some(command.spawn().expect("spawn qp_tunnel")),
            dir: dir.to_path_buf(),
        }
    }

    /// The child's process id (it is alive).
    pub fn pid(&self) -> u32 {
        self.child.as_ref().expect("child is still running").id()
    }

    /// Everything written to stdout so far.
    pub fn stdout(&self) -> String {
        self.read("stdout")
    }

    /// Everything written to stderr so far.
    pub fn stderr(&self) -> String {
        self.read("stderr")
    }

    /// Poll the output file until it contains `needle`.
    ///
    /// # Panics
    ///
    /// After [`TIMEOUT`], or earlier when the child exits without ever
    /// printing the needle (the message carries both streams).
    pub fn wait_for_output(&mut self, needle: &str) {
        let deadline = Instant::now() + TIMEOUT;
        loop {
            if self.stdout().contains(needle) {
                return;
            }
            if let Some(status) = self.try_wait() {
                panic!(
                    "qp_tunnel exited ({status:?}) before printing {needle:?}\nstdout:\n{}\nstderr:\n{}",
                    self.stdout(),
                    self.stderr()
                );
            }
            assert!(
                Instant::now() < deadline,
                "qp_tunnel did not print {needle:?} within {TIMEOUT:?}\nstdout:\n{}\nstderr:\n{}",
                self.stdout(),
                self.stderr()
            );
            thread::sleep(Duration::from_millis(20));
        }
    }

    /// Wait for `ready` (printed once the input thread is in its loop).
    pub fn wait_ready(&mut self) {
        self.wait_for_output("\nready\n");
    }

    /// Send one signal with `/bin/kill` (no `libc` in the tests).
    pub fn signal(&self, name: &str) {
        let status = Command::new("/bin/kill")
            .arg(format!("-{name}"))
            .arg(self.pid().to_string())
            .status()
            .expect("/bin/kill runs");
        assert!(status.success(), "/bin/kill -{name} failed: {status:?}");
    }

    /// Wait for the child to exit, killing it after [`TIMEOUT`].
    ///
    /// # Panics
    ///
    /// When the deadline passes; the child is killed and reaped first so a
    /// failing test never leaves a process behind.
    pub fn wait(&mut self) -> ExitStatus {
        let deadline = Instant::now() + TIMEOUT;
        loop {
            if let Some(status) = self.try_wait() {
                return status;
            }
            if Instant::now() >= deadline {
                let pid = self.pid();
                let _ = Command::new("/bin/kill")
                    .arg("-KILL")
                    .arg(pid.to_string())
                    .status();
                let status = self.child.as_mut().expect("child").wait().expect("reap");
                panic!(
                    "qp_tunnel did not exit within {TIMEOUT:?}; killed ({status:?})\nstdout:\n{}\nstderr:\n{}",
                    self.stdout(),
                    self.stderr()
                );
            }
            thread::sleep(Duration::from_millis(20));
        }
    }

    /// Reap the child if it has exited, returning its status.
    fn try_wait(&mut self) -> Option<ExitStatus> {
        let child = self.child.as_mut()?;
        match child.try_wait().expect("try_wait") {
            Some(status) => {
                self.child = None;
                Some(status)
            }
            None => None,
        }
    }

    /// Read one of the two output files lossily (all reference text is ASCII).
    fn read(&self, name: &str) -> String {
        let bytes = fs::read(self.dir.join(name)).unwrap_or_default();
        String::from_utf8_lossy(&bytes).into_owned()
    }
}

impl Drop for Running {
    fn drop(&mut self) {
        if let Some(child) = self.child.as_mut() {
            // A test that panicked before `wait` must not leak the child.
            let _ = Command::new("/bin/kill")
                .arg("-KILL")
                .arg(child.id().to_string())
                .status();
            let _ = child.wait();
        }
    }
}

/// Find four ports `P, P+1, P+2, P+9` that are free on `127.0.0.1`.
///
/// The search starts at a per-test offset and advances by 100, keeping every
/// candidate in `30000..60000`; the probe listeners are dropped before the
/// window is returned (the caller re-binds it in 0MQ).
pub fn free_remote_window() -> u32 {
    let start = 30_000 + (std::process::id() % 100) * 100;
    for base in (start..60_000).step_by(100) {
        if [0, 1, 2, 9]
            .iter()
            .all(|offset| TcpListener::bind(("127.0.0.1", (base + offset) as u16)).is_ok())
        {
            return base;
        }
    }
    panic!("no free remote window in 30000..60000");
}

/// Create an EZFIO directory whose `work/qp_run_address` is `address`.
///
/// The reference treats a directory as EZFIO without reading `.version`, so
/// the fixture has none; `work/` exists, which is what the reference's
/// `Qpackage.ezfio_work` needs.
pub fn make_ezfio(directory: &Path, address: &str) -> PathBuf {
    let work = directory.join("work");
    fs::create_dir_all(&work).expect("create work/");
    fs::write(work.join("qp_run_address"), format!("{address}\n")).expect("write address");
    directory.to_path_buf()
}

/// Extract the local port from the
/// `Connect to the tunnel with: tcp://<ip>:<port>\n` announcement.
///
/// # Panics
///
/// When the announcement is missing or its endpoint has no colon.
pub fn parse_local_port(announcements: &str) -> u32 {
    let endpoint = announcements
        .strip_prefix("Connect to the tunnel with: ")
        .expect("the first announcement is the tunnel address")
        .lines()
        .next()
        .expect("announcement line");
    let port = endpoint.rsplit(':').next().expect("endpoint has a colon");
    port.parse().expect("the local port is decimal")
}

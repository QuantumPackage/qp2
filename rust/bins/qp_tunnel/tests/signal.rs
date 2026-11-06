//! `SIGUSR1` and `SIGINT` end a running tunnel cleanly.
//!
//! The child is a real process in EZFIO mode whose remote address is
//! `tcp://127.0.0.1:1` (well formed, nothing listening). All three remote
//! `connect` calls are asynchronous and succeed, so every thread reaches its
//! poll loop and `Ready` is printed. Both signals are checked in one test so
//! the two children never race for the same local port window.

mod support;

use std::time::Duration;

use support::{Running, make_ezfio};

/// The signal names `/bin/kill` accepts, checked in one after the other.
const SIGNALS: [&str; 2] = ["USR1", "INT"];

#[test]
fn signals_end_the_tunnel_cleanly() {
    for signal in SIGNALS {
        let directory = tempfile::tempdir().expect("tempdir");
        let ezfio = make_ezfio(directory.path(), "tcp://127.0.0.1:1");
        let mut child = Running::spawn(
            directory.path(),
            &[ezfio.to_str().expect("UTF-8 path")],
            Some("/nonexistent"),
        );
        child.wait_ready();
        assert!(child.stdout().starts_with("Connect to:\ntcp://127.0.0.1:"));

        // `Ready` is printed by the input thread; `main` may not have reached
        // `install_signal_handlers` yet (documented window). Give it time.
        std::thread::sleep(Duration::from_millis(300));
        child.signal(signal);

        let status = child.wait();
        assert_eq!(
            status.code(),
            Some(0),
            "signal {signal}: stderr:\n{}",
            child.stderr()
        );
        assert!(
            child
                .stdout()
                .ends_with("Ready\nqp_tunnel exited properly.\n"),
            "signal {signal}: stdout:\n{}",
            child.stdout()
        );
        assert_eq!(child.stderr(), "", "signal {signal}");
    }
}

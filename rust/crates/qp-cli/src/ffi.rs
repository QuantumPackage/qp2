//! The four `unsafe` libc blocks `qp_run` needs and `std` cannot provide.
//!
//! PI decision (wave 3, 2026-11-05): these four calls are allowed **here and
//! nowhere else** in the workspace, each in its own `// SAFETY:`-documented
//! block:
//!
//! | function | reference site | why `std` cannot do it |
//! |---|---|---|
//! | [`local_time`] | `Unix.localtime` (`qp_run.ml:102`) | `std` has no local-time conversion (`chrono` is not a dependency) |
//! | [`hostname`] | `Unix.gethostname` (`TaskServer.ml:99`) | `std` exposes no host name |
//! | [`ipv4_for_interface`] | `get_ipv4_address_for_interface` C stub (`QP_NIC`) | `std` enumerates no interfaces |
//! | [`install_sigusr1`] | `Sys.set_signal Sys.sigusr1` (`qp_run.ml:88`) | `std` has no signal API |
//!
//! Everything else in this module is safe code: the IPv4 lookup of a *host
//! name* uses [`std::net::ToSocketAddrs`], so the resolver never needs
//! `unsafe`.
//!
//! # Signal safety
//!
//! The SIGUSR1 handler runs on an arbitrary thread, so it only stores `true`
//! in a [`std::sync::atomic::AtomicBool`] (async-signal-safe). The reference
//! runs `Sys.command ("qp_stop "^ezfio)` *inside* its handler; here the
//! caller polls [`take_sigusr1`] and runs the command from a normal
//! thread, because `system(3)` is not async-signal-safe and the task server
//! keeps other threads alive (see [`crate::signal`]).

use std::ffi::CStr;
use std::fmt;
use std::net::Ipv4Addr;
use std::sync::atomic::{AtomicBool, Ordering};

/// A failed OS call (`Unix.Unix_error` in the reference).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Error {
    /// Name of the C function that failed.
    call: &'static str,
    /// `errno` description, or a short explanation when `errno` is unset.
    detail: String,
}

impl Error {
    /// Build an error from a call name and a description.
    fn new(call: &'static str, detail: impl fmt::Display) -> Self {
        Self {
            call,
            detail: detail.to_string(),
        }
    }

    /// The C function that failed.
    #[must_use]
    pub fn call(&self) -> &'static str {
        self.call
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Unix_error({}, {:?})", self.call, self.detail)
    }
}

impl std::error::Error for Error {}

/// Broken-down local time (`Unix.tm`), with the fields `qp_run` prints.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Tm {
    /// Calendar year, e.g. `2026` (`Unix.tm_year + 1900`).
    pub year: i32,
    /// Month, `1..=12` (`Unix.tm_mon + 1`).
    pub mon: u32,
    /// Day of the month, `1..=31`.
    pub mday: u32,
    /// Hour, `0..=23`, **without** the daylight-saving adjustment.
    pub hour: u32,
    /// Minute, `0..=59`.
    pub min: u32,
    /// Second, `0..=59`.
    pub sec: u32,
    /// Whether daylight saving time is in effect.
    pub isdst: bool,
}

/// Local time for `epoch_secs` (`Unix.localtime`).
///
/// # Errors
///
/// [`Error`] when `localtime_r` returns a null pointer (an out-of-range
/// epoch).
pub fn local_time(epoch_secs: i64) -> Result<Tm, Error> {
    let time = epoch_secs as libc::time_t;
    // SAFETY: `time` and `broken` are owned locals. `localtime_r` is the
    // reentrant form: it only writes through `&mut broken` and reads `&time`,
    // both valid for the call. A null return means it did not write, and we
    // discard `broken` in that case.
    let broken = unsafe {
        let mut broken: libc::tm = std::mem::zeroed();
        if libc::localtime_r(&time, &mut broken).is_null() {
            return Err(Error::new("localtime_r", "out of range"));
        }
        broken
    };
    Ok(Tm {
        year: broken.tm_year + 1900,
        mon: (broken.tm_mon + 1).try_into().unwrap_or(0),
        mday: broken.tm_mday.try_into().unwrap_or(0),
        hour: broken.tm_hour.try_into().unwrap_or(0),
        min: broken.tm_min.try_into().unwrap_or(0),
        sec: broken.tm_sec.try_into().unwrap_or(0),
        isdst: broken.tm_isdst > 0,
    })
}

/// The machine's host name (`Unix.gethostname`).
///
/// # Errors
///
/// [`Error`] when `gethostname` fails, or when the name is not UTF-8.
pub fn hostname() -> Result<String, Error> {
    // `HOST_NAME_MAX` is 255 on Linux and 255 on macOS; one extra byte for the
    // terminator so a truncated name is still a C string.
    let mut buffer = [0_i8; 256];
    // SAFETY: `buffer` is a 256-byte writable array and the length passed is
    // its exact size, so `gethostname` writes at most that many bytes and only
    // through that pointer. POSIX does not say whether a truncated name is
    // terminated, so the terminator is forced below before any C-string read.
    let result = unsafe { libc::gethostname(buffer.as_mut_ptr(), buffer.len()) };
    if result != 0 {
        return Err(Error::new("gethostname", last_os_error()));
    }
    // A host name is at most `HOST_NAME_MAX` (255) bytes, so the last byte is
    // normally the terminator already; forcing it keeps the invariant true for
    // a pathologically long name too.
    buffer[buffer.len() - 1] = 0;
    // SAFETY: the array is NUL-terminated (see above) and lives until the end
    // of the statement; `CStr::from_ptr` only scans for the terminator.
    let name = unsafe { CStr::from_ptr(buffer.as_ptr()) };
    name.to_str()
        .map(str::to_owned)
        .map_err(|error| Error::new("gethostname", error))
}

/// First IPv4 address of the interface called `interface`
/// (`get_ipv4_address_for_interface`, the `QP_NIC` path of
/// `TaskServer.ml:106`).
///
/// The reference C stub returns a string starting with `"error"` on failure
/// and the dotted quad on success; [`crate::host::SystemHost`] keeps that
/// convention, so the error is reported as a string here and not as [`Error`].
#[must_use]
pub fn ipv4_for_interface(interface: &str) -> String {
    match find_ipv4(interface) {
        Ok(Some(address)) => address.to_string(),
        Ok(None) => format!("error: no IPv4 address on interface {interface}"),
        Err(error) => format!("error: {error}"),
    }
}

/// Walk the interface list for the first IPv4 address of `interface`.
fn find_ipv4(interface: &str) -> Result<Option<Ipv4Addr>, Error> {
    // SAFETY: `getifaddrs` returns a null-terminated linked list owned by the
    // C library, or a null pointer on failure (checked). Every node is only
    // read, `ifa_addr` is null-checked before it is dereferenced, and
    // `ifa_addr` is reinterpreted as `sockaddr_in` only after `sa_family` says
    // so, which is the layout contract of `getifaddrs(3)`. The list is
    // released with `freeifaddrs` on every path, exactly once.
    unsafe {
        let mut head: *mut libc::ifaddrs = std::ptr::null_mut();
        if libc::getifaddrs(&mut head) != 0 {
            return Err(Error::new("getifaddrs", last_os_error()));
        }
        let mut found = None;
        let mut node = head;
        while !node.is_null() {
            let entry = &*node;
            if !entry.ifa_addr.is_null()
                && (*entry.ifa_addr).sa_family as i32 == libc::AF_INET
                && !entry.ifa_name.is_null()
                && CStr::from_ptr(entry.ifa_name)
                    .to_str()
                    .is_ok_and(|name| name == interface)
            {
                let address = &*(entry.ifa_addr as *const libc::sockaddr_in);
                // `s_addr` is stored in network byte order; `to_ne_bytes`
                // returns the four bytes in memory order, i.e. a.b.c.d.
                found = Some(Ipv4Addr::from(address.sin_addr.s_addr.to_ne_bytes()));
                break;
            }
            node = entry.ifa_next;
        }
        libc::freeifaddrs(head);
        Ok(found)
    }
}

/// Set by the SIGUSR1 handler; see [`install_sigusr1`].
static SIGUSR1_RECEIVED: AtomicBool = AtomicBool::new(false);

/// The SIGUSR1 handler body: one relaxed store, nothing else.
extern "C" fn handle_sigusr1(_signal: libc::c_int) {
    SIGUSR1_RECEIVED.store(true, Ordering::Relaxed);
}

/// Install [`handle_sigusr1`] for SIGUSR1 (`Sys.set_signal Sys.sigusr1`).
///
/// Calling this twice replaces the handler; the flag is not cleared.
///
/// # Errors
///
/// [`Error`] when `sigaction` fails.
pub fn install_sigusr1() -> Result<(), Error> {
    // SAFETY: `action` is an owned local, zeroed first because `sigaction` is
    // a struct with padding and platform-specific trailing fields. The handler
    // is an `extern "C" fn(c_int)` matching `sighandler_t`, `sa_mask` is
    // initialised by `sigemptyset` on the same object, and `SA_RESTART` keeps
    // a pending blocking call (e.g. `waitpid`) from failing with EINTR. The old
    // action is not needed, hence the null third argument.
    unsafe {
        let mut action: libc::sigaction = std::mem::zeroed();
        action.sa_sigaction = handle_sigusr1 as libc::sighandler_t;
        libc::sigemptyset(&mut action.sa_mask);
        action.sa_flags = libc::SA_RESTART;
        if libc::sigaction(libc::SIGUSR1, &action, std::ptr::null_mut()) != 0 {
            return Err(Error::new("sigaction", last_os_error()));
        }
    }
    Ok(())
}

/// Consume the SIGUSR1 flag: `true` if a signal was received since the last
/// call, and the flag is reset in the same atomic step.
///
/// This is the only thing the handler touches, so it is safe to call from any
/// thread at any time. Consuming (rather than reading then clearing) means a
/// signal that arrives between the two calls is not lost.
#[must_use]
pub fn take_sigusr1() -> bool {
    SIGUSR1_RECEIVED.swap(false, Ordering::Relaxed)
}

/// The current `errno` as its description.
fn last_os_error() -> std::io::Error {
    std::io::Error::last_os_error()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Oracle: 2025-11-05T00:00:00Z is 1_762_300_800. The day/month/year are
    /// checked against the BSD/GNU `date` command, which reads the same
    /// time-zone database through a different code path.
    #[test]
    fn local_time_matches_the_system_clock() {
        let tm = local_time(1_762_300_800).expect("localtime_r");
        assert_eq!(tm.year, 2025);
        assert_eq!(tm.mon, 11);
        let formatted = format!("{:02}/{:02}/{:4}", tm.mday, tm.mon, tm.year);
        let mut date = std::process::Command::new("date");
        if cfg!(target_os = "macos") {
            date.args(["-r", "1762300800"]);
        } else {
            date.args(["-d", "@1762300800"]);
        }
        let expected = date.arg("+%d/%m/%Y").output().expect("date");
        assert_eq!(
            formatted,
            String::from_utf8_lossy(&expected.stdout).trim(),
            "date stderr = {}",
            String::from_utf8_lossy(&expected.stderr)
        );
    }

    /// The host name is non-empty and agrees with `hostname(1)`.
    #[test]
    fn hostname_matches_the_command() {
        let name = hostname().expect("gethostname");
        assert!(!name.is_empty());
        let expected = std::process::Command::new("hostname")
            .output()
            .expect("hostname");
        assert_eq!(name, String::from_utf8_lossy(&expected.stdout).trim());
    }

    /// A non-existent interface is an error string starting with `"error"`,
    /// as the reference C stub produces.
    #[test]
    fn unknown_interface_is_an_error_string() {
        let result = ipv4_for_interface("qpnosuchif0");
        assert!(result.starts_with("error"), "result = {result}");
    }

    /// `lo0`/`lo` carries 127.0.0.1 on every platform QP2 builds on.
    #[test]
    fn loopback_interface_has_the_loopback_address() {
        let name = if cfg!(target_os = "macos") {
            "lo0"
        } else {
            "lo"
        };
        assert_eq!(ipv4_for_interface(name), "127.0.0.1");
    }

    /// Installing the handler and receiving the signal sets the flag;
    /// consuming it resets the flag. The process is not killed, which is the
    /// point of installing the handler.
    ///
    /// The signal is sent by `kill(1)` rather than by an `unsafe` `raise`, so
    /// this module stays inside the four approved libc areas.
    #[test]
    fn sigusr1_sets_the_flag() {
        install_sigusr1().expect("sigaction");
        assert!(!take_sigusr1(), "the flag starts clear");
        let status = std::process::Command::new("kill")
            .args(["-USR1", &std::process::id().to_string()])
            .status()
            .expect("run kill(1)");
        assert!(status.success(), "kill -USR1 failed");
        // Delivery is asynchronous: the flag is set by whichever thread runs
        // the handler, so poll instead of reading once.
        let mut received = false;
        for _ in 0..200 {
            if take_sigusr1() {
                received = true;
                break;
            }
            std::thread::sleep(std::time::Duration::from_millis(10));
        }
        assert!(received, "the SIGUSR1 handler did not set the flag");
        assert!(!take_sigusr1(), "consuming clears the flag");
    }
}

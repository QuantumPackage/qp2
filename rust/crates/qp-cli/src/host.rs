//! [`qp_zmq::task_server::HostLookup`] on the real operating system.
//!
//! `TaskServer.ip_address` (`ocaml/TaskServer.ml:99-127`) needs three facts:
//! the machine's host name, the first IPv4 address of that name, and — when
//! `QP_NIC` is set — the first IPv4 address of a named interface. *What* to do
//! with them (the `"localhost"` fallback and the `"error"` prefix test) is
//! [`qp_zmq::task_server::resolve_address`](qp_zmq::task_server::resolve_address),
//! which is pure and already tested; [`SystemHost`] only supplies the answers
//! that touch the operating system.
//!
//! # Resolution order and cost
//!
//! [`SystemHost::ipv4_for_host`] asks the system resolver (NSS: `/etc/hosts`,
//! then DNS), so the first call can block for as long as the resolver takes.
//! `qp_run` resolves the address exactly once per process, from the main
//! thread, before the task server starts (issue #393), so the cost is paid
//! once and never concurrently.
//!
//! # Deviations
//!
//! * `Unix.gethostbyname` returns A records only. [`std::net::ToSocketAddrs`]
//!   returns both families, so the IPv6 entries are discarded and the first
//!   IPv4 entry is returned — the value `h_addr_list.(0)` holds.
//! * The reference's error text for a failed *host* lookup is never printed:
//!   `TaskServer.ip_address` replaces the host name by `"localhost"` on
//!   failure and turns a second failure into
//!   `Unable to find IP address from host name.` ([`SystemHost`] returns the
//!   resolver's message, which `resolve_address` discards).

use std::net::{IpAddr, ToSocketAddrs};

use qp_zmq::task_server::HostLookup;

use crate::ffi;

/// [`HostLookup`] backed by `gethostname`, the system resolver and
/// `getifaddrs`.
///
/// The value has no state: every call re-reads the operating system, exactly
/// like the reference's `Lazy.force hostname` / `Unix.gethostbyname` pair, but
/// without the lazy. A caller that must not resolve twice (issue #393) is
/// responsible for calling `resolve_address` once.
#[derive(Debug, Clone, Copy, Default)]
pub struct SystemHost;

impl HostLookup for SystemHost {
    /// The machine's host name (`Unix.gethostname`).
    ///
    /// # Errors
    ///
    /// The `gethostname` failure or non-UTF-8 name as a message; the caller
    /// substitutes `"localhost"`.
    fn hostname(&self) -> Result<String, String> {
        ffi::hostname().map_err(|error| error.to_string())
    }

    /// First IPv4 address of `host` (`Unix.gethostbyname host` then
    /// `h_addr_list.(0)`).
    ///
    /// A dotted quad is parsed without the resolver; anything else goes
    /// through NSS/DNS.
    ///
    /// # Errors
    ///
    /// The resolver's message when the name does not resolve, or a fixed
    /// message when it resolves to no IPv4 address at all.
    fn ipv4_for_host(&self, host: &str) -> Result<String, String> {
        let addresses = (host, 0_u16)
            .to_socket_addrs()
            .map_err(|error| format!("{host}: {error}"))?;
        addresses
            .filter_map(|address| match address.ip() {
                IpAddr::V4(address) => Some(address.to_string()),
                IpAddr::V6(_) => None,
            })
            .next()
            .ok_or_else(|| format!("{host} has no IPv4 address"))
    }

    /// First IPv4 address of the interface called `interface`
    /// (`get_ipv4_address_for_interface`).
    ///
    /// The C stub reports failure as a string starting with `"error"`, which
    /// is the convention `resolve_address` tests for; `ffi` keeps it.
    fn ipv4_for_interface(&self, interface: &str) -> String {
        ffi::ipv4_for_interface(interface)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::net::Ipv4Addr;
    use std::str::FromStr;

    use qp_zmq::task_server::{ResolveError, resolve_address};

    /// A dotted quad needs no resolver and must come back unchanged, which
    /// also exercises the IPv4 filter of the address list.
    #[test]
    fn numeric_host_needs_no_resolver() {
        assert_eq!(
            SystemHost.ipv4_for_host("127.0.0.1"),
            Ok("127.0.0.1".to_string())
        );
    }

    /// `localhost` is a loopback address on every host QP2 builds on; the
    /// exact quad is not asserted because `/etc/hosts` may list several.
    #[test]
    fn localhost_resolves_to_loopback() {
        let address = SystemHost
            .ipv4_for_host("localhost")
            .expect("localhost must resolve");
        let address = Ipv4Addr::from_str(&address).expect("dotted quad");
        assert!(address.is_loopback(), "address = {address}");
    }

    /// `resolve_address(Some(loopback interface))` is the `QP_NIC` path of
    /// `qp_run`, end to end through the trait.
    #[test]
    fn loopback_interface_is_the_loopback_address() {
        let name = if cfg!(target_os = "macos") {
            "lo0"
        } else {
            "lo"
        };
        assert_eq!(
            resolve_address(&SystemHost, Some(name)),
            Ok("127.0.0.1".to_string())
        );
    }

    /// A missing interface keeps the reference's failure message, which names
    /// the variable's value and not the C stub's text.
    #[test]
    fn unknown_interface_is_reported() {
        assert_eq!(
            resolve_address(&SystemHost, Some("qpnosuchif0")),
            Err(ResolveError::Interface("qpnosuchif0".to_string()))
        );
    }
}

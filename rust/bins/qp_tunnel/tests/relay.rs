//! End-to-end relay checks with a fake remote, in one process.
//!
//! The fake remote binds the three peer sockets the tunnel connects to in
//! EZFIO mode (`REP P`, `PUB P+1`, `REP P+2`) on `127.0.0.1`, echoes every
//! multipart request with an `echo:` prefix, and publishes `news` every
//! 50 ms. The REPs reproduce the real master: the local QP2 client is a REQ,
//! so its identity is carried in the REP envelope that the tunnel's stateless
//! ROUTER→DEALER hop preserves.
//!
//! # The envelope
//!
//! A local REQ sends `[empty, payload…]`; the tunnel's ROUTER prepends the
//! client identity (`[id, empty, payload…]`) and its DEALER forwards every
//! frame, so the fake REP sees the envelope `[id, empty]` and payload
//! `[payload…]`. The reply `[id, empty, echo:payload…]` goes back through the
//! DEALER and is routed by the local ROUTER using `id`; the REQ then strips
//! `[empty]`. This is why the fake is a REP and not a ROUTER. The tunnel's advertised address is passed explicitly
//! (`start_with_ip`), because `QP_NIC` cannot be set safely from a
//! multi-threaded test process.
//!
//! `STOP` is process-wide, so the one tunnel test holds [`SERIAL`]; no other
//! in-process tunnel may run at the same time.
//!
//! The input thread is exercised in EZFIO mode: `test`, `get_ezfio_filename`,
//! an unknown message and `get_input` (a real `tar` of the EZFIO directory).
//! `-g` is requested as well and must be a no-op in this mode.

mod support;

use std::fs;
use std::process::Command;
use std::sync::{Arc, Mutex, mpsc};
use std::thread;
use std::time::{Duration, Instant};

use qp_tunnel::{Options, Target, Tunnel, ssh_instructions};
use support::{free_remote_window, make_ezfio, parse_local_port};
use zmq::{Context, Socket, SocketType};

/// Serializes in-process tunnels: the stop flag is a process-wide static.
static SERIAL: Mutex<()> = Mutex::new(());

/// A fake remote peer: two echoing REPs and a chatty PUB.
struct FakeRemote {
    commands: mpsc::Sender<()>,
    handle: Option<thread::JoinHandle<()>>,
}

impl FakeRemote {
    /// Bind the three peer sockets on `base`, `base+1` and `base+2`.
    fn start(base: u32) -> Self {
        let context = Context::new();
        let ocaml = context.socket(SocketType::REP).expect("REP");
        ocaml
            .bind(&format!("tcp://127.0.0.1:{base}"))
            .expect("bind the ocaml peer");
        let publisher = context.socket(SocketType::PUB).expect("PUB");
        publisher
            .bind(&format!("tcp://127.0.0.1:{}", base + 1))
            .expect("bind the pub peer");
        let fortran = context.socket(SocketType::REP).expect("REP");
        fortran
            .bind(&format!("tcp://127.0.0.1:{}", base + 2))
            .expect("bind the fortran peer");

        let (commands, receiver) = mpsc::channel();
        let handle = thread::spawn(move || {
            let mut items = [
                ocaml.as_poll_item(zmq::POLLIN),
                fortran.as_poll_item(zmq::POLLIN),
            ];
            let mut last_publish: Option<Instant> = None;
            loop {
                match zmq::poll(&mut items, 50) {
                    Ok(_) => {}
                    Err(zmq::Error::EINTR) => continue,
                    Err(error) => panic!("fake remote poll failed: {error}"),
                }
                if items[0].is_readable() {
                    echo(&ocaml);
                }
                if items[1].is_readable() {
                    echo(&fortran);
                }
                // Re-publish until the tunnel's SUB has subscribed (PUB
                // silently drops messages without subscribers).
                if last_publish.is_none_or(|last| last.elapsed() >= Duration::from_millis(50)) {
                    publisher.send(b"news".as_slice(), 0).expect("publish");
                    last_publish = Some(Instant::now());
                }
                if receiver.try_recv().is_ok() {
                    return;
                }
            }
        });
        Self {
            commands,
            handle: Some(handle),
        }
    }

    /// Ask the thread to leave its loop and wait for it.
    fn stop(mut self) {
        let _ = self.commands.send(());
        if let Some(handle) = self.handle.take() {
            handle
                .join()
                .expect("the fake remote thread does not panic");
        }
    }
}

impl Drop for FakeRemote {
    fn drop(&mut self) {
        if let Some(handle) = self.handle.take() {
            let _ = self.commands.send(());
            let _ = handle.join();
        }
    }
}

/// REP echo: each payload frame gains an `echo:` prefix. The envelope is
/// stored and re-attached by the REP socket itself.
fn echo(socket: &Socket) {
    let frames = socket.recv_multipart(0).expect("fake remote recv");
    assert!(
        !frames.is_empty(),
        "a REP request has at least one payload frame"
    );
    let reply: Vec<Vec<u8>> = frames
        .iter()
        .map(|frame| {
            let mut echoed = b"echo:".to_vec();
            echoed.extend_from_slice(frame);
            echoed
        })
        .collect();
    socket.send_multipart(reply, 0).expect("fake remote send");
}

/// Wait until the announcement buffer contains `needle`, then return a copy.
fn wait_for_announcement(state: &Arc<Mutex<String>>, needle: &str) -> String {
    let deadline = Instant::now() + Duration::from_secs(10);
    loop {
        {
            let text = state.lock().expect("announcement lock");
            if text.contains(needle) {
                return text.clone();
            }
        }
        assert!(
            Instant::now() < deadline,
            "the tunnel never announced {needle:?}; got {:?}",
            state.lock().expect("announcement lock")
        );
        thread::sleep(Duration::from_millis(10));
    }
}

/// Connect a local client socket to the tunnel with a receive timeout.
fn connect(context: &Context, kind: SocketType, endpoint: &str) -> Socket {
    let socket = context.socket(kind).expect("client socket");
    socket.set_linger(0).expect("linger");
    socket.set_rcvtimeo(5000).expect("receive timeout");
    socket.connect(endpoint).expect("client connect");
    socket
}

#[test]
fn relays_every_channel_and_serves_the_input_protocol() {
    let _guard = SERIAL
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner());

    let directory = tempfile::tempdir().expect("tempdir");
    let remote = free_remote_window();
    let fake = FakeRemote::start(remote);
    let ezfio = make_ezfio(directory.path(), &format!("tcp://127.0.0.1:{remote}"));

    let announcements = Arc::new(Mutex::new(String::new()));
    let sink = Arc::clone(&announcements);
    let tunnel = Tunnel::start_with_ip(
        Options {
            // `-g` must do nothing in EZFIO mode, as in the reference.
            get_input: true,
            target: Target::Ezfio(ezfio.clone()),
        },
        "127.0.0.1",
        move |text| {
            sink.lock().expect("announcement lock").push_str(text);
            Ok(())
        },
    )
    .expect("the tunnel starts");

    let text = wait_for_announcement(&announcements, "\nready\n");
    let local = parse_local_port(&text);
    assert!(
        text.starts_with(&format!(
            "Connect to the tunnel with: tcp://127.0.0.1:{local}\n"
        )),
        "announcements:\n{text}"
    );
    assert!(
        !text.contains("downloading input"),
        "`-g` must be a no-op in EZFIO mode:\n{text}"
    );
    assert!(
        !text.contains("communication:"),
        "EZFIO mode has no remote input socket:\n{text}"
    );
    let hostname = qp_cli::ffi::hostname().expect("hostname");
    let ssh = ssh_instructions(&hostname, "127.0.0.1", i64::from(remote), local);
    assert!(text.contains(&ssh), "announcements:\n{text}\nssh:\n{ssh}");
    assert!(text.ends_with("ready\n"));

    let context = Context::new();
    let ocaml = connect(
        &context,
        SocketType::REQ,
        &format!("tcp://127.0.0.1:{local}"),
    );
    let fortran = connect(
        &context,
        SocketType::REQ,
        &format!("tcp://127.0.0.1:{}", local + 2),
    );
    let subscriber = connect(
        &context,
        SocketType::SUB,
        &format!("tcp://127.0.0.1:{}", local + 1),
    );
    subscriber.set_subscribe(b"").expect("subscribe");
    let input = connect(
        &context,
        SocketType::REQ,
        &format!("tcp://127.0.0.1:{}", local + 9),
    );

    // OCaml channel, both directions and multipart framing.
    ocaml
        .send_multipart([b"one".as_slice(), b"two".as_slice()], 0)
        .expect("send the ocaml multipart message");
    assert_eq!(
        ocaml.recv_multipart(0).expect("receive the ocaml reply"),
        vec![b"echo:one".to_vec(), b"echo:two".to_vec()]
    );

    // Fortran channel.
    fortran.send(b"F".as_slice(), 0).expect("send fortran");
    assert_eq!(
        fortran.recv_bytes(0).expect("receive the fortran reply"),
        b"echo:F"
    );

    // Pub channel: remote -> local only (a PUB is never readable).
    assert_eq!(
        subscriber
            .recv_multipart(0)
            .expect("receive the pub message"),
        vec![b"news".to_vec()]
    );

    // Input protocol: the three reference actions plus the fallback.
    for (request, expected) in [
        (b"test".as_slice(), b"OK".to_vec()),
        (
            b"get_ezfio_filename".as_slice(),
            ezfio.as_os_str().as_encoded_bytes().to_vec(),
        ),
        (
            b"whatever".as_slice(),
            b"Message 'whatever' not understood".to_vec(),
        ),
    ] {
        input.send(request, 0).expect("send an input request");
        assert_eq!(
            input.recv_bytes(0).expect("receive the input reply"),
            expected
        );
    }

    // `get_input`: one gzipped tar of the EZFIO directory.
    input
        .send(b"get_input".as_slice(), 0)
        .expect("ask for input");
    let tarball = input.recv_bytes(0).expect("receive the tarball");
    assert!(tarball.starts_with(&[0x1f, 0x8b]), "gzip magic");
    let archive = directory.path().join("input.tar.gz");
    fs::write(&archive, &tarball).expect("write the tarball");
    let listing = Command::new("tar")
        .arg("-tzf")
        .arg(&archive)
        .output()
        .expect("tar -tzf runs");
    assert!(listing.status.success());
    let listing = String::from_utf8_lossy(&listing.stdout);
    assert!(
        listing
            .lines()
            .any(|line| line.ends_with("/work/qp_run_address")),
        "listing:\n{listing}"
    );

    Tunnel::request_stop();
    tunnel.join().expect("no thread failed");
    fake.stop();
}

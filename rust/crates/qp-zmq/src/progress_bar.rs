//! The `qp_run` progress bar (`ocaml/Progress_bar.ml`).
//!
//! The task server owns one bar per job. `add_task` extends the end value by
//! one per queued task, `task_done` advances the current value by one per
//! finished task, and the main loop asks for a redraw each time a request
//! arrives. The bar is also the only user of the two output formats:
//!
//! * **tty** — `%s : [%s] %4.1f%% | %8.0f s\r`, one carriage-return line;
//! * **file** — `%5.2f %%  in  %20.0f seconds \n`, one line per update.
//!
//! Both are written to **stderr**; only the *choice* between them looks at
//! stdout (the reference tests `Unix.isatty Unix.stdout`, see
//! `ocaml/Progress_bar.ml:82`). Times are seconds since the epoch, kept as
//! `f64` exactly like OCaml's `Unix.time`.
//!
//! # Reproducibility
//!
//! The reference reads the clock more than once per redraw (`init` twice,
//! `display_file` twice). This port reads it once per entry point and reuses
//! that value, so a redraw is a pure function of `(bar, now, is_tty)`; the
//! delay before the next redraw is therefore exactly 0.1 s (tty) or 10 s
//! (file) after the instant used for the elapsed time. The difference with
//! the reference is bounded by one clock read (microseconds) and never
//! observable in the printed value.
//!
//! Non-finite percentages (only possible if a caller builds a bar whose end
//! value equals its start value) render as Rust's `NaN`/`inf`, whereas C
//! `printf` writes `nan`/`inf`. `qp_run` always uses `start=0`, `end>=1`, so
//! the case is unreachable there.
//!
//! # Counted lines
//!
//! `clear` emits the reference's 75 spaces followed by `\r`
//! (`ocaml/Progress_bar.ml:24`); [`ProgressBar::clear_to`] builds the spaces
//! from a constant so the count cannot drift.

use std::io::{self, IsTerminal, Write};
use std::time::{SystemTime, UNIX_EPOCH};

/// Number of spaces written by the reference's `clear`
/// (`ocaml/Progress_bar.ml:24`).
const CLEAR_WIDTH: usize = 75;

/// Seconds since the Unix epoch as an `f64`, the Rust equivalent of
/// `Unix.time ()`.
///
/// A clock set before 1970 yields `0.0`; the reference would return a
/// negative float. The task server never observes the difference.
fn now_seconds() -> f64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0.0, |elapsed| elapsed.as_secs_f64())
}

/// A progress bar with the reference's fields and redraw cadence.
///
/// The bar is `dirty` after [`ProgressBar::update`],
/// [`ProgressBar::increment_cur`] and [`ProgressBar::increment_end`] (the
/// last one clears it again, as the reference does); a redraw happens only
/// when it is dirty **and** the next-update instant has been reached.
#[derive(Clone, Debug)]
pub struct ProgressBar {
    /// Text printed before the opening bracket.
    title: String,
    /// Value mapped to 0 %.
    start_value: f64,
    /// Current value.
    cur_value: f64,
    /// Value mapped to 100 %.
    end_value: f64,
    /// Number of characters between the brackets.
    bar_length: usize,
    /// Epoch seconds at creation, subtracted from "now" for the elapsed time.
    init_time: f64,
    /// Whether a redraw is pending.
    dirty: bool,
    /// Epoch seconds before which no redraw may happen.
    next: f64,
}

impl ProgressBar {
    /// Create a bar with the reference defaults: `start=0`, `end=1`,
    /// `bar_length=20`, the current time as origin, not dirty.
    pub fn new(title: impl Into<String>) -> Self {
        Self::new_at(title, now_seconds())
    }

    /// [`ProgressBar::new`] with an explicit origin time, in epoch seconds.
    ///
    /// This is the seam used by tests to make the elapsed-time column
    /// deterministic; the task server calls [`ProgressBar::new`].
    pub fn new_at(title: impl Into<String>, now: f64) -> Self {
        Self::new_with(title, 0.0, 1.0, 20, now)
    }

    /// Full constructor: `start_value` maps to 0 %, `end_value` to 100 %,
    /// `bar_length` characters are printed, and `now` is the origin time in
    /// epoch seconds (also used as the next-update instant).
    ///
    /// The reference's optional arguments are `~start_value:0.`,
    /// `~end_value:1.`, `~bar_length:20`.
    pub fn new_with(
        title: impl Into<String>,
        start_value: f64,
        end_value: f64,
        bar_length: usize,
        now: f64,
    ) -> Self {
        Self {
            title: title.into(),
            start_value,
            cur_value: start_value,
            end_value,
            bar_length,
            init_time: now,
            dirty: false,
            next: now,
        }
    }

    /// Replace the current value and mark the bar dirty
    /// (`Progress_bar.update`, `ocaml/Progress_bar.ml:17`).
    ///
    /// The reference defines it but never calls it; the task server uses
    /// [`ProgressBar::increment_cur`]/[`ProgressBar::increment_end`] only.
    /// Kept for API parity and used by the unit tests to set arbitrary
    /// progress values that the incrementors cannot express.
    #[allow(dead_code)]
    pub fn update(&mut self, cur_value: f64) {
        self.cur_value = cur_value;
        self.dirty = true;
    }

    /// Advance the current value by one and mark the bar dirty
    /// (`Progress_bar.increment_cur`).
    pub fn increment_cur(&mut self) {
        self.cur_value += 1.0;
        self.dirty = true;
    }

    /// Extend the end value by one and **clear** the dirty flag
    /// (`Progress_bar.increment_end`).
    ///
    /// The reference clears the flag because a newly queued task must not
    /// trigger a redraw by itself; the next `task_done` sets it again.
    pub fn increment_end(&mut self) {
        self.end_value += 1.0;
        self.dirty = false;
    }

    /// Whether a redraw is pending (the reference's `bar.dirty`).
    #[must_use]
    pub fn is_dirty(&self) -> bool {
        self.dirty
    }

    /// Erase the bar from a tty: 75 spaces then `\r`, on stderr.
    ///
    /// Mirrors `Progress_bar.clear`, which returns `None` so its caller can
    /// store `progress_bar = None`; callers here drop the [`ProgressBar`]
    /// themselves.
    ///
    /// # Errors
    ///
    /// Propagates the `stderr` write error (the reference's `eprintf` would
    /// raise `Sys_error`).
    pub fn clear() -> io::Result<()> {
        let mut stderr = io::stderr();
        Self::clear_to(&mut stderr)
    }

    /// [`ProgressBar::clear`] into an arbitrary writer, for tests.
    ///
    /// # Errors
    ///
    /// Propagates the writer error.
    pub fn clear_to<W: Write + ?Sized>(out: &mut W) -> io::Result<()> {
        out.write_all(&[b' '; CLEAR_WIDTH])?;
        out.write_all(b"\r")?;
        out.flush()
    }

    /// Redraw on stderr if the bar is dirty and due.
    ///
    /// The tty/file decision is `stdout.is_terminal()`, as the reference
    /// tests `isatty stdout` even though it writes to stderr.
    ///
    /// # Errors
    ///
    /// Propagates the `stderr` write error.
    pub fn display(&mut self) -> io::Result<()> {
        let now = now_seconds();
        let is_tty = io::stdout().is_terminal();
        let mut stderr = io::stderr();
        self.display_at(&mut stderr, now, is_tty)
    }

    /// Redraw into `out` at time `now`, with the tty/file choice supplied by
    /// the caller.
    ///
    /// Does nothing when the bar is not dirty or when `now` is before the
    /// next-update instant. On a redraw, clears the dirty flag and schedules
    /// the next one at `now + 0.1` (tty) or `now + 10` (file), in epoch
    /// seconds.
    ///
    /// # Errors
    ///
    /// Propagates the writer error.
    pub fn display_at<W: Write + ?Sized>(
        &mut self,
        out: &mut W,
        now: f64,
        is_tty: bool,
    ) -> io::Result<()> {
        if !self.dirty || now < self.next {
            return Ok(());
        }
        if is_tty {
            let line = self.tty_line(now);
            out.write_all(line.as_bytes())?;
            self.next = now + 0.1;
        } else {
            let line = self.file_line(now);
            out.write_all(line.as_bytes())?;
            self.next = now + 10.0;
        }
        self.dirty = false;
        out.flush()
    }

    /// Completion ratio in percent, computed with the reference's operation
    /// order so the rounded digits match.
    fn percent(&self) -> f64 {
        100.0 * (self.cur_value - self.start_value) / (self.end_value - self.start_value)
    }

    /// The tty rendering of `Progress_bar.display_tty`, terminated by `\r`.
    fn tty_line(&self, now: f64) -> String {
        let percent = self.percent();
        let bar_length = self.bar_length as f64;
        // `Float.to_int` truncates toward zero; the `as` cast does too.
        let n_hashes = (bar_length * percent / 100.0) as i64;
        let mut hashes = String::with_capacity(self.bar_length);
        for i in 0..self.bar_length {
            hashes.push(if (i as i64) < n_hashes { '#' } else { ' ' });
        }
        let running_time = now - self.init_time;
        format!(
            "{} : [{}] {:4.1}% | {:8.0} s\r",
            self.title, hashes, percent, running_time
        )
    }

    /// The file rendering of `Progress_bar.display_file`, terminated by
    /// `\n`.
    fn file_line(&self, now: f64) -> String {
        let percent = self.percent();
        let running_time = now - self.init_time;
        format!("{:5.2} %  in  {:20.0} seconds \n", percent, running_time)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Oracle: the reference format strings, rendered by
    /// `Printf.sprintf` in OCaml 5 (see the work log); percent 25, elapsed 0.
    #[test]
    fn tty_line_matches_ocaml_golden() {
        let mut bar = ProgressBar::new_at("scf", 1000.0);
        bar.update(0.25);
        let mut out = Vec::new();
        bar.display_at(&mut out, 1000.0, true).unwrap();
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "scf : [#####               ] 25.0% |        0 s\r"
        );
    }

    /// Oracle: `Printf.printf "%5.2f %%  in  %20.0f seconds \n" 25. 0.`
    /// (the extra `%` in the OCaml source is the escaped literal).
    #[test]
    fn file_line_matches_ocaml_golden() {
        let mut bar = ProgressBar::new_at("scf", 1000.0);
        bar.update(0.25);
        let mut out = Vec::new();
        bar.display_at(&mut out, 1000.0, false).unwrap();
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "25.00 %  in                     0 seconds \n"
        );
    }

    /// The golden pair for 12.34 % at 5.6 s, where `%5.2f` pads and
    /// `%20.0f` rounds 5.6 to 6.
    #[test]
    fn rounding_matches_ocaml_golden() {
        let mut bar = ProgressBar::new_at("scf", 0.0);
        bar.update(0.1234);
        let mut tty = Vec::new();
        bar.display_at(&mut tty, 5.6, true).unwrap();
        assert_eq!(
            String::from_utf8(tty).unwrap(),
            "scf : [##                  ] 12.3% |        6 s\r"
        );
        let mut bar = ProgressBar::new_at("scf", 0.0);
        bar.update(0.1234);
        let mut file = Vec::new();
        bar.display_at(&mut file, 5.6, false).unwrap();
        assert_eq!(
            String::from_utf8(file).unwrap(),
            "12.34 %  in                     6 seconds \n"
        );
    }

    /// A fresh bar is not dirty; `increment_cur` dirties it and
    /// `increment_end` clears it again.
    #[test]
    fn dirty_flags() {
        let mut bar = ProgressBar::new_at("scf", 0.0);
        assert!(!bar.is_dirty());
        let mut out = Vec::new();
        bar.display_at(&mut out, 1.0, true).unwrap();
        assert!(out.is_empty(), "clean bar must not redraw");

        bar.increment_cur();
        assert!(bar.is_dirty());
        bar.increment_end();
        assert!(!bar.is_dirty());
        assert_eq!(bar.cur_value, 1.0);
        assert_eq!(bar.end_value, 2.0);
    }

    /// The next redraw is gated at +0.1 s on a tty and +10 s otherwise,
    /// measured from the `now` passed to the redraw.
    #[test]
    fn redraw_cadence() {
        let mut bar = ProgressBar::new_at("scf", 0.0);
        bar.update(0.5);
        let mut out = Vec::new();
        bar.display_at(&mut out, 100.0, true).unwrap();
        assert!(!out.is_empty());
        out.clear();
        bar.display_at(&mut out, 100.05, true).unwrap();
        assert!(out.is_empty());
        bar.update(0.6);
        bar.display_at(&mut out, 100.1, true).unwrap();
        assert!(!out.is_empty());

        let mut bar = ProgressBar::new_at("scf", 0.0);
        bar.update(0.5);
        let mut out = Vec::new();
        bar.display_at(&mut out, 100.0, false).unwrap();
        assert!(!out.is_empty());
        out.clear();
        bar.display_at(&mut out, 109.99, false).unwrap();
        assert!(out.is_empty());
        bar.update(0.6);
        bar.display_at(&mut out, 110.0, false).unwrap();
        assert!(!out.is_empty());
    }

    /// Hashes clamp naturally: over 100 % fills the bar, a negative ratio
    /// leaves it empty.
    #[test]
    fn hash_count_is_truncated_and_clamped() {
        let mut bar = ProgressBar::new_with("x", 0.0, 1.0, 4, 0.0);
        bar.update(1.5);
        assert_eq!(bar.tty_line(0.0), "x : [####] 150.0% |        0 s\r");
        bar.update(-0.5);
        assert_eq!(bar.tty_line(0.0), "x : [    ] -50.0% |        0 s\r");
    }

    /// `clear` is exactly 75 spaces and a carriage return, and flushes.
    #[test]
    fn clear_is_75_spaces() {
        let mut out = Vec::new();
        ProgressBar::clear_to(&mut out).unwrap();
        assert_eq!(out.len(), 76);
        assert_eq!(&out[..75], &[b' '; 75]);
        assert_eq!(out[75], b'\r');
    }
}

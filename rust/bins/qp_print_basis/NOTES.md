# P4 tools — dual-run gate recipe (`qp_print_basis` + `element_create_db`)

**Status (2026-11-05):** both ports are implemented and the gate below is
green — 14 compared cases byte-identical on `(exit code, stdout, stderr)`,
plus one informational deviation (broken pipe). This file is the recipe the
parent re-runs after merging; it is self-contained and uses only files that
are in the repository (no scratch captures required).

```
PASS noarg            exit 2  out      0 B  err   63 B
PASS nonexistent      exit 2  out      0 B  err   59 B
PASS emptystr         exit 2  out      0 B  err   59 B
PASS h2               exit 2  out   1103 B  err   59 B
PASS hcn              exit 0  out  62691 B  err    0 B
PASS extra            exit 0  out  62691 B  err    0 B
PASS emptydir         exit 2  out      0 B  err   73 B
PASS plainfile        exit 2  out      0 B  err   74 B
PASS noao             exit 2  out      0 B  err   58 B
PASS ecd-write        exit 0  out      0 B  err    0 B
PASS ecd-table             110 lines,     2987 B, md5 54a2ab1991fa42279a8ef75d0f8d91ed
PASS ecd-unset        exit 2  out      0 B  err  119 B
PASS ecd-nodata       exit 2  out      0 B  err   92 B
PASS ecd-emptyroot    exit 2  out      0 B  err   86 B
INFO broken-pipe      reference exit 141 (stderr 0 B), port exit 2
GATE: all compared cases identical
```

## Prerequisites

* Reference binaries live in `~/qp2/ocaml/`, **not** in `~/qp2/bin`:
  `~/qp2/ocaml/qp_print_basis` and `~/qp2/ocaml/element_create_db`.
* Build the ports from the worktree's `rust/` directory:

  ```sh
  cargo build -j4 -p qp_print_basis -p element_create_db
  ```

  (`cargo -j4 build` is rejected by cargo 1.92; the flag goes after the
  subcommand. A shared 8 GiB machine means `-j4`, always.)
* Both reference binaries instantiate `Qpackage` at module-init time, so they
  need `QP_ROOT` **even for `--help`-less, argument-less runs**. The
  `qp_print_basis` port never reads `QP_ROOT`; the gate sets it for both sides
  anyway so the environments match.
* `element_create_db` *writes* `$QP_ROOT/data/list_element.txt`. The gate
  therefore gives both sides a **fake, relative** root (`QP_ROOT=root`, one
  per side) and never points it at `~/qp2`. This is also why the recipe does
  not source `quantum_package.rc` in a way that exports the real root into the
  child: a stray `QP_ROOT` would overwrite `~/qp2/data/list_element.txt`.

## What is compared, and why

Each side runs in its own directory (`gate/ref`, `gate/port`) and receives
**relative** paths. The reference's `Unix.EEXIST` message quotes the argument
verbatim, so relative paths make the two stderr streams comparable without
normalisation; absolute paths would differ by the side's directory.

The EZFIO fixtures are copied **per side** before use. Opening an EZFIO tree
in the reference restamps `ao_basis/ao_md5` and materialises `ao_cartesian`,
so sharing one copy would let the first run mutate the second's input. The
copies under `gate/` are scratch; `rust/tests/golden/fixtures/` is never
written (verified with `git status` after the run).

| case | `argv[1]` | reference behaviour pinned |
| --- | --- | --- |
| `noarg` | — | `Sys.argv.(1)` → `Invalid_argument("index out of bounds")`, exit 2 |
| `nonexistent` | `does-not-exist` | `Failure("Error reading EZFIO file")`, exit 2 |
| `emptystr` | `""` | same as `nonexistent` (empty path does not exist) |
| `h2` | copy of fixture | basis block on stdout, then `Failure("Error reading the mo set")`, exit 2 |
| `hcn` | copy of fixture | 62 691 B (AO + MO blocks), exit 0 |
| `extra` | `hcn ignored1 ignored2` | later arguments are ignored |
| `emptydir` | existing empty dir | `Unix.Unix_error(Unix.EEXIST, "mkdir", "emptydir")`, exit 2 |
| `plainfile` | regular file | same `EEXIST` shape; `mkdir` is attempted on any non-EZFIO path |
| `noao` | fixture minus `ao_basis/` | `Failure("Error reading basis set")`, exit 2 |
| `ecd-write` | fake root with `data/` | silent, exit 0; table byte-diffed against both sides and the in-repo golden |
| `ecd-unset` | `QP_ROOT` unset | `Failure("QP_ROOT environment variable is not set.\nPlease source the quantum_package.rc file.")`, exit 2 |
| `ecd-nodata` | `QP_ROOT=nodata` (no `data/`) | `Sys_error("nodata/data/list_element.txt: No such file or directory")`, exit 2 |
| `ecd-emptyroot` | `QP_ROOT=""` | `Sys_error("/data/list_element.txt: No such file or directory")`, exit 2 |

Reference captures for the two full runs are pinned as repository goldens:
`bins/qp_print_basis/tests/data/h2_basis.rst` (1103 B, md5
`794c04208fee52bae8af21b9aeb9c74a`) and
`bins/element_create_db/tests/data/list_element.txt` (2987 B, 110 lines, md5
`54a2ab1991fa42279a8ef75d0f8d91ed`). `hcn` has no full golden in the repo; the
gate compares reference vs port and pins the length, while
`bins/qp_print_basis/tests/action.rs` pins its first 200 bytes
(`tests/data/hcn_basis_prefix.rst`).

## Known, deliberate deviations (reported, not failed)

* **Broken pipe.** `qp_print_basis hcn | head -c 1`: the reference dies of
  `SIGPIPE` (exit 141, empty stderr) because OCaml leaves the signal at its
  default disposition. Rust ignores `SIGPIPE`, so the port sees `EPIPE` and
  reports `Fatal error: exception Sys_error("Broken pipe (os error 32)")`,
  exit 2. Nothing in the reference contract distinguishes the two for a
  file/terminal stdout.
* **`.version` as a directory.** `qp_ezfio::EzFio::open` requires
  `<dir>/.version` to be a regular file and reports "not an EZFIO directory";
  the reference only calls `Sys.file_exists`. No EZFIO writer produces a
  `.version` directory, and the gate does not exercise it.
* **No `ao_md5` restamp.** The port reads but does not write `ao_md5` /
  `ao_cartesian`; only the opened tree differs, never stdout, stderr or the
  exit status. Covered by the port's action test
  `the_printed_tree_is_the_copy_and_stays_unchanged`.

## Gotchas that cost time

* `~/qp2/bin` has no copies of these two tools; the binaries are in
  `~/qp2/ocaml/`. Using `bin/` silently finds nothing.
* Run the script with **bash**. In zsh the pipeline status is
  `$pipestatus[1]`, not `${PIPESTATUS[0]}`; using the bash spelling expands to
  the empty string and hides the 141-vs-2 difference.
* `head -c 10` on a 62 691-byte write never closes the pipe: the whole block
  fits the kernel buffer and no `SIGPIPE` is delivered. Use a reader that
  stops after one byte.
* All scratch (this gate included) lives under
  `$HOME/OpenCode/tmp/w3-p4-tools/print-basis/`, never `/tmp`.

## The script (copy into scratch and run with `bash`)

```bash
#!/bin/bash
# Dual-run gate: OCaml reference vs Rust port, byte-for-byte.
# Scratch only; never writes into ~/qp2 (the fake roots are here).
set -u

G=$HOME/OpenCode/tmp/w3-p4-tools/print-basis/gate
WT=<port-checkout>/rust
BIN=$WT/target/debug
SRC=$WT/tests/golden/fixtures
H2GOLD=$WT/bins/qp_print_basis/tests/data/h2_basis.rst
ECOGOLD=$WT/bins/element_create_db/tests/data/list_element.txt

REF=$HOME/qp2/ocaml/qp_print_basis
PORT=$BIN/qp_print_basis
REFECD=$HOME/qp2/ocaml/element_create_db
PORTECD=$BIN/element_create_db

fail=0

report() { # name ref_exit port_exit
    local name=$1 r=$2 p=$3 verdict=OK
    cmp -s "$G/ref/$name.out" "$G/port/$name.out" || verdict=STDOUT
    cmp -s "$G/ref/$name.err" "$G/port/$name.err" || verdict="$verdict/STDERR"
    [ "$r" = "$p" ] || verdict="$verdict/EXIT($r!=$p)"
    if [ "$verdict" = OK ]; then
        printf 'PASS %-16s exit %d  out %6d B  err %4d B\n' "$name" "$r" \
            "$(wc -c <"$G/ref/$name.out")" "$(wc -c <"$G/ref/$name.err")"
    else
        printf 'FAIL %-16s -> %s\n' "$name" "$verdict"
        fail=1
    fi
}

rm -rf "$G/ref" "$G/port"
mkdir -p "$G/ref" "$G/port"

# ---------------------------------------------------------------- fixtures --
for side in ref port; do
    cp -R "$SRC/h2" "$G/$side/h2"
    cp -R "$SRC/hcn" "$G/$side/hcn"
    cp -R "$SRC/h2" "$G/$side/noao"
    rm -rf "$G/$side/noao/ao_basis"
    mkdir "$G/$side/emptydir"
    printf 'hello\n' >"$G/$side/plainfile"
done

# ------------------------------------------------------- qp_print_basis ----
# Same relative path on both sides, so paths inside messages are comparable.
run_case() { # name args...
    local name=$1
    shift
    local r p
    (cd "$G/ref" && QP_ROOT="$HOME/qp2" "$REF" "$@" \
        >"$G/ref/$name.out" 2>"$G/ref/$name.err")
    r=$?
    (cd "$G/port" && QP_ROOT="$HOME/qp2" "$PORT" "$@" \
        >"$G/port/$name.out" 2>"$G/port/$name.err")
    p=$?
    report "$name" "$r" "$p"
}

run_case noarg
run_case nonexistent does-not-exist
run_case emptystr ""
run_case h2 h2
run_case hcn hcn
run_case extra hcn ignored1 ignored2
run_case emptydir emptydir
run_case plainfile plainfile
run_case noao noao

# The byte-identical port output must also equal the in-repo goldens.
cmp -s "$G/ref/h2.out" "$H2GOLD" || {
    echo "FAIL h2 stdout != bins/qp_print_basis/tests/data/h2_basis.rst"
    fail=1
}
[ "$(wc -c <"$G/ref/h2.out")" -eq 1103 ] || {
    echo "FAIL h2 stdout is not 1103 B"
    fail=1
}
[ "$(wc -c <"$G/ref/hcn.out")" -eq 62691 ] || {
    echo "FAIL hcn stdout is not 62691 B"
    fail=1
}

# -------------------------------------------------- element_create_db ------
# Writing case: two fresh fake roots, same relative name, then byte-diff of
# the generated tables and of both against the reference capture.
for side in ref port; do
    rm -rf "$G/$side/root"
    mkdir -p "$G/$side/root/data"
done
(cd "$G/ref" && QP_ROOT=root "$REFECD" >"$G/ref/ecd-write.out" 2>"$G/ref/ecd-write.err")
r=$?
(cd "$G/port" && QP_ROOT=root "$PORTECD" >"$G/port/ecd-write.out" 2>"$G/port/ecd-write.err")
p=$?
report ecd-write "$r" "$p"
if cmp -s "$G/ref/root/data/list_element.txt" "$G/port/root/data/list_element.txt" &&
    cmp -s "$G/ref/root/data/list_element.txt" "$ECOGOLD"; then
    echo "PASS ecd-table        $(wc -l <"$G/ref/root/data/list_element.txt") lines, $(wc -c <"$G/ref/root/data/list_element.txt") B, md5 $(md5 -q "$G/ref/root/data/list_element.txt")"
else
    echo "FAIL ecd-table        generated tables differ"
    fail=1
fi
head -1 "$G/ref/root/data/list_element.txt"
tail -1 "$G/ref/root/data/list_element.txt"

# Failure cases: same relative QP_ROOT on both sides.
ecd_case() { # name mode value
    local name=$1 mode=$2 value=${3-} r p
    if [ "$mode" = unset ]; then
        (cd "$G" && env -u QP_ROOT "$REFECD" \
            >"$G/ref/$name.out" 2>"$G/ref/$name.err")
        r=$?
        (cd "$G" && env -u QP_ROOT "$PORTECD" \
            >"$G/port/$name.out" 2>"$G/port/$name.err")
        p=$?
    else
        (cd "$G" && QP_ROOT="$value" "$REFECD" \
            >"$G/ref/$name.out" 2>"$G/ref/$name.err")
        r=$?
        (cd "$G" && QP_ROOT="$value" "$PORTECD" \
            >"$G/port/$name.out" 2>"$G/port/$name.err")
        p=$?
    fi
    report "$name" "$r" "$p"
}

ecd_case ecd-unset unset
ecd_case ecd-nodata set nodata
ecd_case ecd-emptyroot set ""

# ------------------------------------------------ broken pipe (deviation) --
# `head -c 1` closes after one byte; the remaining 62 690 B cannot be written.
# The reference dies of SIGPIPE (exit 141, empty stderr); Rust ignores SIGPIPE
# and reports Sys_error. Documented in the module docs, so this check only
# prints. PIPESTATUS is read inside the command substitution that owns the
# pipeline.
r=$(cd "$G/ref" && QP_ROOT="$HOME/qp2" "$REF" hcn 2>"$G/ref/pipe.err" |
    head -c 1 >/dev/null
    echo "${PIPESTATUS[0]}")
p=$(cd "$G/port" && QP_ROOT="$HOME/qp2" "$PORT" hcn 2>"$G/port/pipe.err" |
    head -c 1 >/dev/null
    echo "${PIPESTATUS[0]}")
echo "INFO broken-pipe      reference exit $r (stderr $(wc -c <"$G/ref/pipe.err") B), port exit $p"
echo "     port stderr: $(cat "$G/port/pipe.err")"

echo
if [ "$fail" -eq 0 ]; then
    echo "GATE: all compared cases identical"
else
    echo "GATE: FAILURES"
fi
exit "$fail"
```

Note on `ecd-nodata` / `ecd-emptyroot`: `Sys_error` quotes the string the
reference concatenated (`root ^ "/data/list_element.txt"`), which is why the
empty root yields the doubled-slash-free `/data/list_element.txt` and why a
`Path::join`-style normalisation in the port would break both cases.

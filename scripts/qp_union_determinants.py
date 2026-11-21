#!/usr/bin/env python3
"""
qp_union_determinants.py -- replace the selected determinants of several
EZFIO files by their union.

Use case: potential energy curves with CIPSI. A small CIPSI calculation
is done at each geometry (each point in its own EZFIO directory), and the
selected determinant spaces differ from one geometry to the next, which
makes the energy curve noisy.  This script takes the union of the
determinant occupations found in all given EZFIO files and replaces, in
each one, the stored variational space (the `determinants` group:
n_det, psi_det, psi_coef) with that common set.  Determinants that were
not present in a given geometry keep their occupation, with a zero
coefficient as initial value.  Re-diagonalizing the Hamiltonian in each
geometry on the identical space (read_wf True, see below) then gives a
smooth curve.

Usage:
    python3 qp_union_determinants.py [--dry-run] EZFIO_DIR [EZFIO_DIR ...]

Each EZFIO_DIR is a Quantum Package EZFIO directory (e.g. mycalc.ezfio)
holding the result of a CIPSI calculation, i.e. it must contain
determinants/psi_det.gz and determinants/psi_coef.gz.

The files about to be modified are backed up first as <name>.orig
(existing .orig backups are left untouched, so the pre-union original is
preserved across repeated runs).

Next step (per directory):
    qp set_file mycalc.ezfio
    qp set determinants read_wf True
    qp run fci
    qp run pt2        # optional but useful: on the read_wf path fci does
                      # NOT write fci/energy (save_energy is only called
                      # inside the selection loop); pt2 re-diagonalizes
                      # the stored WF, computes the PT2 correction and
                      # saves both fci/energy and fci/energy_pt2.
The script raises determinants.n_det_max to the union size when it is
smaller, and since the CIPSI selection loop only runs while
N_det < n_det_max, the fci run then performs a pure re-diagonalization
of the common determinant space (no new determinants are selected).
Raising n_det_max again would resume the selection from the common
space.

The on-disk format used here is the native EZFIO text format:
  <group>/<name>.gz is a gzip-compressed text file holding
    line 1: the rank (number of dimensions)
    line 2: the dimensions (for psi_det: n_chunks, 2, n_det)
    then one value per line.  In psi_det the values are ordered by
    determinant (slowest), then spin, then chunk (fastest); in
    psi_coef by determinant, then state.  The chunks are integer(8)
    stored SIGNED: a determinant with an electron in orbital 63 of a
    chunk carries a NEGATIVE value (bit 63 set) - any bit counting
    must mask to the 64-bit unsigned pattern first.
Scalars (n_det, n_det_max, ...) are plain text integers.
"""

import argparse
import gzip
import os
import shutil
import sys
import tempfile


def die(msg):
    print(f"qp_union_determinants: error: {msg}", file=sys.stderr)
    sys.exit(1)


def read_int(path):
    try:
        with open(path) as f:
            return int(f.readline().strip())
    except (OSError, ValueError) as e:
        die(f"cannot read integer from {path}: {e}")


def write_int(path, value):
    tmp = tempfile.mkstemp(dir=os.path.dirname(path))[1]
    with open(tmp, "w") as f:
        f.write("%20d\n" % value)
    os.replace(tmp, path)


def read_array(path):
    """Read an EZFIO gzipped text array. Returns (dims, flat list of values)."""
    try:
        with gzip.open(path, "rt") as f:
            lines = f.read().splitlines()
    except OSError as e:
        die(f"cannot read {path}: {e}")
    if len(lines) < 2:
        die(f"{path}: not an EZFIO array file (no header)")
    try:
        rank = int(lines[0].strip())
        dims = [int(x) for x in lines[1].split()]
    except ValueError:
        die(f"{path}: malformed EZFIO array header")
    if len(dims) != rank:
        die(f"{path}: rank {rank} but {len(dims)} dimensions in header")
    expected = 1
    for d in dims:
        expected *= d
    values = lines[2:]
    if len(values) != expected:
        die(f"{path}: header says {expected} values, file has {len(values)}")
    return dims, values


def read_int_array(path):
    dims, values = read_array(path)
    try:
        return dims, [int(x) for x in values]
    except ValueError:
        die(f"{path}: non-integer value in array")


def read_double_array(path):
    dims, values = read_array(path)
    try:
        return dims, [float(x) for x in values]
    except ValueError:
        die(f"{path}: non-numeric value in array")


def write_array(path, dims, values, fmt):
    """Write an EZFIO gzipped text array (values in the same order the
    file stores them: for psi_det determinant, spin, chunk; for
    psi_coef determinant, state)."""
    tmp = tempfile.mkstemp(dir=os.path.dirname(path))[1]
    with gzip.open(tmp, "wt") as f:
        f.write("%3d\n" % len(dims))
        for d in dims:
            f.write("%20d " % d)
        f.write("\n")
        for v in values:
            f.write(fmt % v + "\n")
    os.replace(tmp, path)


def backup(path):
    """Copy path to path + '.orig', unless the backup already exists."""
    if not os.path.exists(path + ".orig"):
        shutil.copy2(path, path + ".orig")


def popcount(x):
    """Popcount of a 64-bit chunk value.  QP2 stores integer(8) chunks
    as SIGNED int64: a chunk with bit 63 set (an electron in orbital 63
    of that chunk) is negative in the file, so mask to the 64-bit
    unsigned pattern first (bit_kind_size = 64 is fixed in QP2)."""
    return bin(x & ((1 << 64) - 1)).count("1")


def read_ezfio_dir(name):
    """Read the determinant selection stored in one EZFIO directory."""
    det_dir = os.path.join(name, "determinants")
    if not os.path.isdir(det_dir):
        die(f"{name}: not an EZFIO directory (no 'determinants' group)")
    for f in ("n_det", "psi_det.gz", "psi_coef.gz"):
        if not os.path.exists(os.path.join(det_dir, f)):
            die(f"{name}: missing determinants/{f}; "
                f"did you run a CIPSI calculation in this directory?")

    n_det = read_int(os.path.join(det_dir, "n_det"))
    n_states = read_int(os.path.join(det_dir, "n_states"))
    n_int = read_int(os.path.join(det_dir, "n_int"))
    bit_kind = read_int(os.path.join(det_dir, "bit_kind"))
    n_det_max = read_int(os.path.join(det_dir, "n_det_max"))
    try:
        n_det_qp_edit = read_int(os.path.join(det_dir, "n_det_qp_edit"))
    except SystemExit:
        n_det_qp_edit = None

    mo_label = None
    mo_label_path = os.path.join(det_dir, "mo_label")
    if os.path.exists(mo_label_path):
        with open(mo_label_path) as f:
            mo_label = f.readline().strip()

    elec = {}
    for key in ("elec_num", "elec_alpha_num", "elec_beta_num"):
        path = os.path.join(name, "electrons", key)
        if os.path.exists(path):
            elec[key] = read_int(path)
        else:
            elec[key] = None

    if n_det <= 0:
        die(f"{name}: n_det = {n_det} (empty wave function)")

    # psi_det: the file header carries the logical dims (n_chunks, 2,
    # n_det), but the values are stored with the determinant index first
    # and the chunk index fastest: for determinant i, spin s (0 = up,
    # 1 = down), chunk c: flat[(i*2 + s)*n_chunks + c].  (Verified
    # against the OCaml `qp get determinants psi_det` output.)
    dims, flat = read_int_array(os.path.join(det_dir, "psi_det.gz"))
    n_chunks = n_int * bit_kind // 8
    if dims != [n_chunks, 2, n_det]:
        die(f"{name}: psi_det dims {dims} inconsistent with "
            f"n_int={n_int}, bit_kind={bit_kind}, n_det={n_det}")
    if n_chunks > 1:
        print(f"qp_union_determinants: note: {name}: n_int = {n_int} "
              f"(> 64 orbitals); verified with C2H6/aug-cc-pVDZ "
              f"(100 MOs, 2 chunks, union of 10 CIPSI selections)",
              file=sys.stderr)

    # psi_coef: dims (n_det, n_states), state index fastest
    cdims, coef_flat = read_double_array(os.path.join(det_dir, "psi_coef.gz"))
    if cdims != [n_det, n_states]:
        die(f"{name}: psi_coef dims {cdims} inconsistent with "
            f"n_det={n_det}, n_states={n_states}")

    # Build the occupation -> coefficients map.  The key of a determinant
    # is its occupation: spin-up chunks then spin-down chunks, in the
    # same order the file stores them.
    det_to_coef = {}
    nbad = 0
    for i in range(n_det):
        base = (2 * i) * n_chunks
        key = tuple(flat[base + j] for j in range(2 * n_chunks))
        if key in det_to_coef:
            nbad += 1
        else:
            det_to_coef[key] = [coef_flat[i * n_states + k]
                                for k in range(n_states)]
        if elec["elec_alpha_num"] is not None:
            if sum(popcount(flat[base + c]) for c in range(n_chunks)) \
                    != elec["elec_alpha_num"]:
                nbad += 1
            if sum(popcount(flat[base + n_chunks + c])
                   for c in range(n_chunks)) \
                    != elec["elec_beta_num"]:
                nbad += 1
    if nbad:
        die(f"{name}: {nbad} determinants with duplicates or electron "
            f"count mismatch (corrupt or inconsistent selection)")

    # Optional display copies used by qp_edit
    qe = None
    qe_path = os.path.join(det_dir, "psi_det_qp_edit.gz")
    if os.path.exists(qe_path) and os.path.exists(
            os.path.join(det_dir, "psi_coef_qp_edit.gz")):
        qdims, qflat = read_int_array(qe_path)
        if qdims == [n_chunks, 2, n_det_qp_edit]:
            _, qcoef = read_double_array(
                os.path.join(det_dir, "psi_coef_qp_edit.gz"))
            qe = (n_det_qp_edit, qflat, qcoef)

    return {
        "name": name,
        "det_dir": det_dir,
        "n_det": n_det,
        "n_states": n_states,
        "n_chunks": n_chunks,
        "n_int": n_int,
        "bit_kind": bit_kind,
        "n_det_max": n_det_max,
        "n_det_qp_edit": n_det_qp_edit,
        "mo_label": mo_label,
        "elec": elec,
        "dets": det_to_coef,
        "qe": qe,
    }


def write_union(d, union_keys):
    """Write the union wave function into one EZFIO directory."""
    det_dir = d["det_dir"]
    n_det = len(union_keys)
    n_chunks = d["n_chunks"]
    n_states = d["n_states"]
    dets = d["dets"]

    # psi_det: for each determinant, spin s (0 = up, 1 = down), chunk c
    # (chunk index fastest) -- the key tuple is already in that order.
    psi_det_flat = []
    for key in union_keys:
        psi_det_flat.extend(key)
    write_array(os.path.join(det_dir, "psi_det.gz"),
                [n_chunks, 2, n_det], psi_det_flat, "%20d")

    # psi_coef: (n_det, n_states), coefficient 0.0 where the determinant
    # was not part of this geometry's original selection
    psi_coef_flat = []
    for i in range(n_det):
        coef = dets.get(union_keys[i])
        if coef is None:
            psi_coef_flat.extend([0.0] * n_states)
        else:
            psi_coef_flat.extend(coef)
    write_array(os.path.join(det_dir, "psi_coef.gz"),
                [n_det, n_states], psi_coef_flat, "%24.15E")

    write_int(os.path.join(det_dir, "n_det"), n_det)

    if d["n_det_max"] < n_det:
        write_int(os.path.join(det_dir, "n_det_max"), n_det)
        d["n_det_max_raised_to"] = n_det

    # Refresh the qp_edit display copies so qp edit shows the new space
    if d["qe"] is not None:
        n_qe_old, _, _ = d["qe"]
        n_qe = min(n_qe_old, n_det)
        qflat = []
        for key in union_keys[:n_qe]:
            qflat.extend(key)
        write_array(os.path.join(det_dir, "psi_det_qp_edit.gz"),
                    [n_chunks, 2, n_qe], qflat, "%20d")
        qcoef = []
        for i in range(n_qe):
            coef = dets.get(union_keys[i])
            qcoef.extend(coef if coef is not None else [0.0] * n_states)
        write_array(os.path.join(det_dir, "psi_coef_qp_edit.gz"),
                    [n_qe, n_states], qcoef, "%24.15E")
        if d["n_det_qp_edit"] is not None and d["n_det_qp_edit"] != n_qe:
            write_int(os.path.join(det_dir, "n_det_qp_edit"), n_qe)
            d["n_det_qp_edit"] = n_qe


def verify(d, union_keys):
    """Re-read what was written and check it matches the union."""
    det_dir = d["det_dir"]
    n_det = read_int(os.path.join(det_dir, "n_det"))
    if n_det != len(union_keys):
        die(f"{d['name']}: verification failed: n_det = {n_det} "
            f"expected {len(union_keys)}")
    dims, flat = read_int_array(os.path.join(det_dir, "psi_det.gz"))
    if dims != [d["n_chunks"], 2, n_det]:
        die(f"{d['name']}: verification failed: psi_det dims {dims}")
    n_chunks = d["n_chunks"]
    keys = set()
    for i in range(n_det):
        base = (2 * i) * n_chunks
        keys.add(tuple(flat[base + j] for j in range(2 * n_chunks)))
    if keys != set(union_keys):
        die(f"{d['name']}: verification failed: determinant set changed "
            f"on disk (missing {len(set(union_keys) - keys)}, "
            f"unexpected {len(keys - set(union_keys))})")
    cdims, coef_flat = read_double_array(os.path.join(det_dir, "psi_coef.gz"))
    if cdims != [n_det, d["n_states"]]:
        die(f"{d['name']}: verification failed: psi_coef dims {cdims}")
    if len(coef_flat) != n_det * d["n_states"]:
        die(f"{d['name']}: verification failed: psi_coef has "
            f"{len(coef_flat)} values, expected {n_det * d['n_states']}")


def main():
    parser = argparse.ArgumentParser(
        description="Replace the selected determinants of several EZFIO "
                    "files by the union of all their determinants "
                    "(smooth CIPSI potential energy curves).")
    parser.add_argument("dirs", nargs="+",
                        help="EZFIO directories holding CIPSI results "
                             "(at least two)")
    parser.add_argument("--dry-run", action="store_true",
                        help="report the union without modifying anything")
    args = parser.parse_args()

    names = []
    for d in args.dirs:
        d = os.path.normpath(os.path.abspath(d))
        if d in names:
            print(f"qp_union_determinants: note: {d} given more than once; "
                  f"ignoring duplicates", file=sys.stderr)
            continue
        names.append(d)
    if len(names) < 2:
        die("at least two distinct EZFIO directories are required")

    data = [read_ezfio_dir(n) for n in names]

    # Consistency gates: the union only makes sense if all directories
    # describe the same electron system on the same MO space
    gates = {
        "n_states": (d["n_states"] for d in data),
        "n_int": (d["n_int"] for d in data),
        "bit_kind": (d["bit_kind"] for d in data),
        "electrons": (tuple(d["elec"].values()) for d in data),
        "mo_label": (d["mo_label"] for d in data),
    }
    for what, vals in gates.items():
        vals = list(vals)
        if any(v != vals[0] for v in vals):
            detail = "; ".join(
                f"{d['name']}: {v}" for d, v in zip(data, vals))
            die(f"directories are not compatible ({what} differs): {detail}")

    # Union of the occupations, in a canonical (sorted) order so that all
    # output files carry the determinants in the same order
    union_keys = sorted(set().union(*(d["dets"].keys() for d in data)))
    n_union = len(union_keys)

    print(f"qp_union_determinants: union of {len(data)} directories: "
          f"{n_union} determinants")
    for d in data:
        own = len(d["dets"])
        added = n_union - own
        note = ""
        if d["n_det_max"] < n_union:
            note = (f"  [n_det_max will be raised {d['n_det_max']} -> "
                    f"{n_union}]")
        print(f"  {d['name']}: {own} -> {n_union} determinants "
              f"(+{added} with zero coefficient){note}")

    if args.dry_run:
        print("qp_union_determinants: --dry-run, nothing modified")
        return

    for d in data:
        # Backup the files that are about to be modified
        for f in ("n_det", "psi_det.gz", "psi_coef.gz"):
            backup(os.path.join(d["det_dir"], f))
        if d["n_det_max"] < n_union:
            backup(os.path.join(d["det_dir"], "n_det_max"))
        if d["qe"] is not None:
            for f in ("psi_det_qp_edit.gz", "psi_coef_qp_edit.gz"):
                backup(os.path.join(d["det_dir"], f))
            if d["n_det_qp_edit"] is not None and \
                    min(d["n_det_qp_edit"], n_union) != d["n_det_qp_edit"]:
                backup(os.path.join(d["det_dir"], "n_det_qp_edit"))

        write_union(d, union_keys)
        verify(d, union_keys)
        print(f"  {d['name']}: updated and verified")

    print("\nNext, re-diagonalize each geometry on the common space, e.g.:")
    for d in data:
        print(f"  cd {os.path.dirname(d['name'])} && "
              f"qp set_file {os.path.basename(d['name'])} && "
              f"qp set determinants read_wf True && qp run fci && "
              f"qp run pt2")
    print("(the trailing 'qp run pt2' is optional but recommended: on the")
    print("read_wf path fci does not write fci/energy - save_energy is only")
    print("called inside the selection loop - while pt2 re-diagonalizes the")
    print("stored WF and saves fci/energy and fci/energy_pt2.)")
    print("n_det_max is set to the union size, so fci re-diagonalizes the")
    print("common space without selecting new determinants; raising")
    print("n_det_max again would resume the CIPSI selection from it.")


if __name__ == "__main__":
    main()

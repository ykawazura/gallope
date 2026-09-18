#!/usr/bin/env python3
"""Machine-check the TSUBAME4.0 scaling data before any of it is plotted.

Nothing here trusts a file name.  Every quantity that decides where a point
lands on a graph -- grid size, rank count, and above all the NODE PLACEMENT --
is read back out of the file the benchmark wrote, and compared against what the
path says it should be.  On this machine the placement IS the experiment (the
whole point is 4 GPUs in one node versus 4 GPUs spread over four), so a run that
landed differently than intended is not a small error: it is the wrong
measurement wearing the right file name.

That failure is not hypothetical.  In the Miyabi data,

    cuFFTMp-pencil-202607/2048^3/node1-2

contains "nx = 1024" -- cuFFTMp had quietly fallen back to a grid that fit, and
the directory kept saying 2048.  It had to be caught by hand at plot time.

Usage:
    verify_tsubame.py <dir> [<dir> ...]

where each <dir> is a retrieved "results-202607/tsubame" tree.  Exits non-zero
if any file fails, and prints exactly which check failed and what was found.
"""

import argparse
import pathlib
import re
import sys
from collections import defaultdict

# Ranks per node used for the CPU curves; both CPU codes were run on cpu_160.
CPU_PPN = 160

# Timed steps the campaign asks the two applications for (NSTEP in the job
# scripts).  A run that stopped short still writes a timer summary that parses
# cleanly, so the step count has to be checked rather than assumed.
NSTEP = 100

# The physics model both applications must be solving.  Not a formality: a
# gallope binary built from the KRMHD branch ran the whole 384^3 sweep to
# "Program completed", with a plausible timer block and the right grid in the
# header, and every check above passed -- the only tell was this line saying
# KRMHD.  Makefile.in said MHD_INCOMP; the build had failed and left the
# previous binary in place.  See lessons L28.
MODEL = "MHD_INCOMP"

# ---------------------------------------------------------------- parsing


def _num(pattern, text, cast=int):
    m = re.search(pattern, text)
    return cast(m.group(1)) if m else None


def _local_shape(text):
    """The real-space block one rank holds, as printed by the benchmarks.

    Slab prints "local_rshape", pencil "local_rshape_in", and both are

        (2*(nz/2+1),  ny_local,  nx_local)

    -- note the reversal: cufftmp_r2c.f90 declares "nx slowest", so nx is the
    distributed axis and nz is the contiguous one carrying the in-place r2c
    padding.  The slab splits nx only (nx_local = nx/nproc, ny_local = ny);
    the pencil splits both by nranks1d = sqrt(nproc).  Either way

        ny_local * nx_local * nproc == ny * nx

    which is the one printed quantity that proves the data really was spread
    over `nproc` ranks.  "GPUs per node" does not: see the note in check().
    """
    m = re.search(r"local_rshape(?:_in)?\s*:\s*(\d+)\s+(\d+)\s+(\d+)", text)
    return tuple(int(m.group(i)) for i in (1, 2, 3)) if m else None


def parse_cufftmp(text):
    """cuFFTMp slab and pencil (config/job.sh-tsubame).

    The two variants are told apart by their own banner rather than by which
    directory the file came from: they write the same relative paths under
    their own results tree, so a checker keyed on the path merges them and
    reports the slab/pencil difference as trial-to-trial scatter.  It did --
    768^3 on 4 GPUs came out as "6 trials, spread 90.8%".
    """
    return {
        "kind": "cufftmp-pencil" if "pencil test" in text else "cufftmp-slab",
        "nx": _num(r"nx\s*=\s*(\d+)", text),
        "ny": _num(r"ny\s*=\s*(\d+)", text),
        "nz": _num(r"nz\s*=\s*(\d+)", text),
        "nproc": _num(r"Total MPI processes:\s*(\d+)", text),
        "nodes": _num(r"Unique nodes:\s*(\d+)", text),
        "ndev": _num(r"GPUs per node:\s*(\d+)", text),
        "local": _local_shape(text),
        "correct": "Results are correct" in text,
        "time": _num(r"cpu time per loop\s*=\s*([0-9.EDed+-]+)", text, float),
    }


def parse_p3dfft(text):
    """P3DFFT (driver_rand.F90)."""
    return {
        "kind": "p3dfft",
        "nx": _num(r"nx\s*=\s*(\d+)", text),
        "ny": _num(r"ny\s*=\s*(\d+)", text),
        "nz": _num(r"nz\s*=\s*(\d+)", text),
        "nproc": _num(r"Total MPI processes:\s*(\d+)", text),
        "nodes": _num(r"Unique nodes:\s*(\d+)", text),
        "threads": _num(r"Running on\s*(\d+) threads", text),
        "correct": "Results are correct" in text,
        "time": _num(r"cpu time per loop\s*=\s*([0-9.EDed+-]+)", text, float),
    }


def _model(text):
    m = re.search(r"^Solving\s+(\S+)", text, re.M)
    return m.group(1) if m else None


def parse_gallope(text):
    """gallope (batch/job.sh-tsubame).  Prints its own placement.

    "time" is the Advance-steps timer, not "total from timer is:", because that
    is the quantity the figure plots (step_time = Advance steps / # of steps).
    The total also carries initialisation and the unconditional end-of-run
    save_restart, neither of which the figure shows; on the Miyabi reference
    runs save_restart cost 0.000 min only because restart/ did not exist and
    the unchecked MPI_FILE_OPEN failed silently.
    """
    return {
        "kind": "gallope",
        "nx": _num(r"nlx\s*=\s*(\d+)", text),
        "ny": _num(r"nly\s*=\s*(\d+)", text),
        "nz": _num(r"nlz\s*=\s*(\d+)", text),
        "nproc": _num(r"MPI processes:\s*(\d+)", text),
        "nodes": _num(r"Nodes:\s*(\d+)", text),
        "ndev": _num(r"GPUs per node:\s*(\d+)", text),
        "steps": _num(r"# of steps advanced\s*(\d+)", text),
        "model": _model(text),
        "namelist_failed": re.findall(r"Reading (\w+) failed", text),
        "correct": "Program completed" in text,
        "time": _num(r"Advance steps\s+([0-9.]+) min", text, float),
    }


def parse_calliope(text):
    """calliope (batch/job.sh-tsubame-scaling).

    calliope does not print a node count -- it only knows nproc, iproc, jproc.
    The node count is therefore checked indirectly, through nproc == 160*nodes,
    and 'nodes' is left None so the report says so rather than implying the
    code confirmed a placement it never saw.

    Nor does it print gallope's "Program completed" banner: it simply stops
    after the timer summary.  Reaching "Final dt" is the last thing a finished
    run writes, so that is what completion means here -- checked against the
    real Miyabi output rather than assumed to match gallope.
    """
    return {
        "kind": "calliope",
        "nx": _num(r"nlx\s*=\s*(\d+)", text),
        "ny": _num(r"nly\s*=\s*(\d+)", text),
        "nz": _num(r"nlz\s*=\s*(\d+)", text),
        "nproc": _num(r"nproc\s*=\s*(\d+)", text),
        "nodes": None,
        "threads": _num(r"Running on\s*(\d+) threads", text),
        "steps": _num(r"# of steps advanced\s*(\d+)", text),
        "model": _model(text),
        "namelist_failed": re.findall(r"Reading (\w+) failed", text),
        # nlxc < nlx means pruned = .true., i.e. only the dealiased band is
        # transformed -- about 0.3 of the work per step, and not what the
        # Miyabi reference measured.  See lessons L27.
        "nxc": _num(r"nlxc\s*=\s*(\d+)", text),
        "correct": "Final dt" in text,
        # Advance steps, for the reason given in parse_gallope.
        "time": _num(r"Advance steps\s+([0-9.]+) min", text, float),
    }


def sniff(text):
    if "cuFFTMp" in text:
        return parse_cufftmp(text)
    if "P3DFFT test" in text:
        return parse_p3dfft(text)
    if "Gallope started" in text:
        return parse_gallope(text)
    if "nproc" in text and "iproc" in text:
        return parse_calliope(text)
    return None


# ---------------------------------------------------------------- expectations

NAME_GPU = re.compile(r"^gpu(\d+)-node(\d+)-(\d+)$")
NAME_RANK = re.compile(r"^rank(\d+)-node(\d+)-(\d+)$")
NAME_NODE = re.compile(r"^node(\d+)-(\d+)$")


def expected_from_path(path, root):
    """What the file's location claims the run was.

    Returns (mode, base_grid, nproc, nodes, trial, label) or None if the path
    is not one this campaign writes.
    """
    rel = path.relative_to(root).parts
    if len(rel) < 3:
        return None
    mode = rel[0]                       # strong_scaling | weak_scaling
    if mode not in ("strong_scaling", "weak_scaling"):
        return None
    m = re.match(r"^(\d+)\^3$", rel[1])
    if not m:
        return None
    base = int(m.group(1))
    label = rel[2] if len(rel) > 3 else ""   # gallope's BUILD / "calliope"
    name = rel[-1]

    for regex, kind in ((NAME_GPU, "gpu"), (NAME_RANK, "rank")):
        mm = regex.match(name)
        if mm:
            return mode, base, int(mm.group(1)), int(mm.group(2)), int(mm.group(3)), label
    mm = NAME_NODE.match(name)
    if mm:                              # calliope: ranks implied by the node count
        nodes = int(mm.group(1))
        return mode, base, CPU_PPN * nodes, nodes, int(mm.group(2)), label
    return None


def expected_grid(mode, base, nproc):
    """The box the run should have used.

    Strong scaling holds the global box; weak scaling holds base^3 per rank and
    doubles the global box x -> y -> z, which is only exact for powers of two --
    the same rule the job scripts apply, restated here independently so a bug in
    one is not silently reproduced by the other.
    """
    if mode == "strong_scaling":
        return base, base, base
    k = nproc.bit_length() - 1
    if 1 << k != nproc:
        return None                     # not a power of two: no weak point exists
    return base << ((k + 2) // 3), base << ((k + 1) // 3), base << (k // 3)


# ---------------------------------------------------------------- checking


def check(path, root):
    """Return (record, [failures])."""
    exp = expected_from_path(path, root)
    if exp is None:
        return None, []
    mode, base, e_nproc, e_nodes, trial, label = exp

    text = path.read_text(errors="replace")
    got = sniff(text)
    if got is None:
        return None, [(path, "unrecognised output format (truncated or crashed?)")]

    fails = []

    def bad(msg):
        fails.append((path, msg))

    if not got["correct"]:
        bad("no correctness marker -- the run did not verify or did not finish")

    e_grid = expected_grid(mode, base, e_nproc)
    g_grid = (got["nx"], got["ny"], got["nz"])
    if e_grid is None:
        bad(f"weak point at nproc={e_nproc} is not a power of two; it should not exist")
    elif None in g_grid:
        bad("grid size missing from the header")
    elif tuple(e_grid) != g_grid:
        bad(f"grid is {g_grid[0]}x{g_grid[1]}x{g_grid[2]}, path says "
            f"{e_grid[0]}x{e_grid[1]}x{e_grid[2]}")

    if got["nproc"] != e_nproc:
        bad(f"nproc is {got['nproc']}, path says {e_nproc}")

    # The placement check.  This is the one the whole TSUBAME campaign exists for.
    if got["nodes"] is None:
        if got["kind"] != "calliope":
            bad("node count missing from the header")
    elif got["nodes"] != e_nodes:
        bad(f"ran on {got['nodes']} nodes, path says {e_nodes}")

    # "GPUs per node" is cudaGetDeviceCount in both codes (cufftmp_r2c.f90:71,
    # gallope src/cuFFTmp.F90:44): it is what the allocated node TYPE has, not
    # what this run used.  A 3-rank job on node_f prints 4.  So the only thing
    # it can be checked for is oversubscription -- more ranks on a node than it
    # has devices, which would put two ranks on one GPU and quietly halve the
    # measurement.  Both codes bind rank -> device as mod(rank, ndevices), so
    # ranks_per_node <= ndevices is exactly the condition for one rank per GPU.
    if got.get("ndev") is not None and e_nodes:
        ppn = -(-e_nproc // e_nodes)            # ceil, the packing the scripts use
        if ppn > got["ndev"]:
            bad(f"{ppn} ranks on a node with {got['ndev']} devices: "
                "mod(rank, ndevices) puts two ranks on one GPU")

    # What "GPUs per node" cannot tell us, the local block can: if the data
    # really was cut into e_nproc pieces, one rank's real-space block times the
    # rank count rebuilds the global ny*nz plane.  This is the check that would
    # have caught the Miyabi pencil run that silently fell back to nx = 1024.
    if got.get("local") is not None and e_grid is not None:
        pad, ny_l, nx_l = got["local"]
        if ny_l * nx_l * e_nproc != e_grid[0] * e_grid[1]:
            bad(f"local block {ny_l}x{nx_l} over {e_nproc} ranks covers "
                f"{ny_l * nx_l * e_nproc} of {e_grid[0] * e_grid[1]} points "
                "in the nx-ny plane: the decomposition is not the one asked for")
        if pad != 2 * (e_grid[2] // 2 + 1):
            bad(f"contiguous extent {pad}, expected {2 * (e_grid[2] // 2 + 1)} "
                f"for nz={e_grid[2]}: the transform was planned on another grid")

    if got.get("threads") is not None and got["threads"] != 1:
        bad(f"running on {got['threads']} OpenMP threads, expected 1 "
            "(the Miyabi reference data is single-threaded)")

    if got.get("steps") is not None and got["steps"] != NSTEP:
        bad(f"advanced {got['steps']} steps, expected {NSTEP}")

    # The application checks that "Program completed" cannot make.  A binary
    # solving the wrong model, or reading a namelist that is not its own,
    # finishes just as cleanly as a correct one.
    if got["kind"] in ("gallope", "calliope"):
        if got.get("model") != MODEL:
            bad(f"solving {got.get('model')}, expected {MODEL} "
                "-- stale binary or wrong branch (lessons L28)")
        # physical_parameters carries nu, eta, shear and q: if it did not read,
        # the dissipation and the shear are the code's defaults, not the run's.
        # force_parameters is expected to fail -- these runs are undriven.
        for nml in got.get("namelist_failed") or []:
            if nml != "force_parameters":
                bad(f"'Reading {nml} failed' -- the input does not match "
                    "this binary; its values are defaults, not the ones asked for")
        if got.get("nxc") is not None and got["nxc"] != got["nx"]:
            bad(f"nlxc={got['nxc']} < nlx={got['nx']}: pruned transform, "
                "which is not what the Miyabi reference measured (lessons L27)")

    if got["time"] is None:
        bad("no timing line -- nothing to plot from this file")

    rec = dict(got, mode=mode, base=base, trial=trial, label=label,
               nproc_exp=e_nproc, nodes_exp=e_nodes, path=path)
    return rec, fails


# ---------------------------------------------------------------- reporting


def report(records, fails):
    groups = defaultdict(list)
    for r in records:
        groups[(r["kind"], r["label"], r["mode"], r["base"], r["nproc_exp"],
                r["nodes_exp"])].append(r)

    print(f"{'code':<15} {'build':<9} {'mode':<7} {'base':>6} {'nproc':>6} "
          f"{'node':>5} {'n':>2}  {'median':>12}  {'min':>12}  {'max':>12}  spread")
    print("-" * 109)
    for key in sorted(groups):
        kind, label, mode, base, nproc, nodes = key
        ts = sorted(r["time"] for r in groups[key] if r["time"] is not None)
        if not ts:
            continue
        med = ts[len(ts) // 2] if len(ts) % 2 else 0.5 * (ts[len(ts) // 2 - 1] +
                                                          ts[len(ts) // 2])
        spread = (ts[-1] - ts[0]) / med * 100 if med else 0.0
        flag = "" if len(ts) == 3 else f"   <-- {len(ts)} trials, expected 3"
        print(f"{kind:<15} {label:<9} {mode[:6]:<7} {base:>6} {nproc:>6} "
              f"{nodes:>5} {len(ts):>2}  {med:12.6g}  {ts[0]:12.6g}  {ts[-1]:12.6g}  "
              f"{spread:5.1f}%{flag}")

    # Weak scaling only means anything if the work per rank really is constant.
    print()
    for (kind, label, mode, base, nproc, nodes), rs in sorted(groups.items()):
        if mode != "weak_scaling":
            continue
        r = rs[0]
        per = r["nx"] * r["ny"] * r["nz"] / nproc
        if abs(per - base ** 3) > 1:
            print(f"!! weak point {kind} {base}^3 nproc={nproc}: "
                  f"{per:.0f} points per rank, expected {base**3}")

    if fails:
        print(f"\n{len(fails)} FAILED CHECK(S):")
        for path, msg in fails:
            print(f"  {path}\n      {msg}")
    else:
        print(f"\nall {len(records)} files pass.")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("roots", nargs="+", type=pathlib.Path,
                    help='retrieved "results-202607/tsubame" trees')
    args = ap.parse_args()

    records, fails, skipped = [], [], []
    for root in args.roots:
        if not root.is_dir():
            sys.exit(f"not a directory: {root}")
        for path in sorted(root.rglob("*")):
            if not path.is_file() or path.name.startswith("."):
                continue
            rec, f = check(path, root)
            fails.extend(f)
            if rec is not None:
                records.append(rec)
            elif not f:
                skipped.append(path)

    if skipped:
        print(f"({len(skipped)} file(s) ignored: not campaign output)\n")
    if not records and not fails:
        sys.exit("no result files found -- is the tree the right one?")
    report(records, fails)
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()

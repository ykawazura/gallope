# runs/ — the Gallope half of the determinism comparison, as executed

Copied verbatim from Miyabi-G (`~/work/determinism/`) before the machine was
decommissioned on 2026-08-13.  The protocol, the reasoning and the result
tables are in `../README.md`; this directory is the raw material behind them.
The Calliope half is in `calliope_dev/tests/determinism/runs/`.

## Layout

| directory | what ran | count |
|---|---|---|
| `rg/` | MHD_INCOMP, `run_{shear,noshear}_n{1,10,100,1000}_g{1,4}` plus `run_shear_n2000_g{1,4}` | 18 |
| `probe/` | MHD_INCOMP convergence probes: `{probe,order,lin}_dt{1e3,5e4,25e5}_gallope` plus `probe_lin_gallope` | 10 |
| `rmhd/` | RMHD, `run_n{1,10,100,1000}_g{1,4}` plus `probe_dt{1e3,5e4,25e5}_gallope` | 11 |

Each run directory holds the input it was started from, its `out.std`, and its
`rms.dat` / `cfl.dat` time series.  The PBS job logs (`det-gallope.o*`,
`rmhd-gallope.o*`, `det-{shear,order,lin}probe.o*`) sit one level up, next to
the driver that submitted them.

## Results

| file | covers |
|---|---|
| `rg/compare_all_fixed.out` | the post-fix MHD_INCOMP comparison — the numbers quoted in `../README.md` |
| `rg/compare_all.out` | the same comparison run *before* the pressure fix, kept for contrast |
| `probe/probe_all_fixed.out` | self-convergence, shear and linear-order probes, post-fix |
| `rmhd/compare_rmhd.out` | the RMHD comparison (also committed at `../rmhd/compare_rmhd.out`) |

## Two versions of compare.py

`rg/compare.py` and `probe/compare.py` are byte-identical to each other but
**not** to `../compare.py`.  The copies here are the MHD_INCOMP-only version
that actually produced `compare_all_fixed.out` and `probe_all_fixed.out`; the
harness copy one level up is the later revision that added `--model RMHD` and
the `omg = -kprp^2*phi` identity check in place of the solenoidality check.
The MHD_INCOMP code path is unchanged between them, but the version of record
for those two output files is the one archived here.

`balance.py`, `divchk.py` and `modes.py` existed only on Miyabi; they have been
promoted to `../` so that they are version controlled, and are left here too so
this directory stays a faithful copy of what was on disk.

## What was not copied

`restart/*.dat` — the spectral field dumps the comparison actually reads — is
about 100 MB per run and was left on Miyabi.  Re-running `compare.py` against
this archive is therefore not possible; the verdicts it reached are the files
listed above.  `out2d/`, `out3d/`, `*.modes.out`, the NetCDF diagnostics and the
run binaries were dropped for the same reason.

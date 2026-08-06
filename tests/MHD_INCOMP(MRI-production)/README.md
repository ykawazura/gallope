# MHD_INCOMP(MRI-production) — nonlinear MRI production run, Gallope half

The shearing-box MRI run behind Fig. 4 and Section 5 of COMPHY-D-26-00110,
rerun after the shearing-box pressure fix (`49e52a2` in Gallope, `92f4047` in
Calliope).  The Calliope half of the same comparison, started from the same
physical and numerical parameters, is in
`calliope_dev/tests/MHD_INCOMP(MRI-production)/`.

This is the *production* run, not a test.  The small regression run that
exercises the same code path lives in `tests/MHD_INCOMP(MRI)/` (64x64x32).

## Why the rerun

The original Fig. 4 predates the pressure fix.  In shearing coordinates the
incompressibility constraint is imposed on the sheared wavenumber
`kxt = kx + q*tsc*ky`, whose time dependence contributes a `-q*ky*ux` term
that the old pressure solve dropped (`src/model/MHD_INCOMP/advance.F90`,
`get_nonlinear_terms`).  The sheared path was therefore effectively first
order in time, while the manuscript describes it as third-order SSP-IF-RK3.
Rerunning removes that inconsistency and puts every number in Section 5 on
post-fix data.

## Run definition

Everything is in `gallope.in`, which carries its own header; the essentials:

| | |
|---|---|
| grid | 512 x 512 x 256 |
| box | `(Lx, Ly, Lz) = 2*pi * (2, 4, 1)`, radial x azimuthal x vertical |
| shear | `shear = .true.`, `q = 1.5` (Keplerian) |
| net vertical flux | `b0 = (0, 0, 5e-2)` |
| seed | `init_type = 'random'`, `\|k_j\| <= 3`, component amplitude 5e-2 |
| scheme | `eSSPIFRK3`, `cfl = 2.0`, `reset_method = 'decrement'` |
| dissipation | `nu = eta = 0`; `nu_h = eta_h = 2.00e2`, `nu_h_exp = eta_h_exp = 2` |

`nu_h` is the damping **rate at the grid scale**, not a coefficient
multiplying `k^(2n)`: the codes apply `nu*(k^2/k2_max) + nu_h*(k^2/k2_max)^n`
(`advance.F90`, `get_imp_terms_tintg`).  The value 2.00e2 at `n = 2` for this
grid comes from the 2023 hyperdissipation scan
(`calliope_results/.../5lmd-2Lx4Lx1L_512x512x256_nexp2_2.00d+2`).  The
original production input had been lost, which is why Section 5 of the first
revision could only state the functional form; with this run the numeric
values can be quoted.

## Placement

`job.pbs-prod-gallope`: 16 nodes, one MPI rank and one GH200 per node, queue
`regular-g` (the scheduler routes it to `small-g` by node count).  The script
carries no `:ompthreads=1` — with it the 16-rank run hangs between
`Solving MHD_INCOMP` and `Initialization done`; without it the same binary,
grid and input reach the main loop in 10 s.  Measured cost at this grid and
node count is 0.087 s/step.

The run is wall-time bounded rather than step bounded (`max_wall_time = 47.0`),
and the script archives each finished leg into `legs/leg<N>/` before
resubmitting, because `src/model/MHD_INCOMP/io.F90` opens the NetCDF
diagnostics with `NF90_CLOBBER` and has no restart branch — a second leg would
otherwise destroy the first leg's time history.  There is no `legs/` here: the
run finished inside a single leg.

## Outcome

Reached `t = 920.0` in 842 889 steps, then **hung** and was killed by PBS
("walltime 172868 exceeded limit 172800").  Every output stream stops at the
same instant — 05:40:41 for `rms.dat`, `cfl.dat`, the NetCDF files and
`gallope.modes.out`, 05:40:42 for `out2d/`, 05:40:43 for `restart/` and
`out.std` — so the `t = 920` block completed and the hang is on the step after
it; 29 h of silence followed.  42 stale MPI-IO lock files were left in
Miyabi's `out2d/` (Calliope's had none) and the Gallope source contains no lock
logic of its own, so they come from the Lustre/MPI-IO layer.  No stack trace
exists, because the process was `SIGKILL`ed.

The data is unaffected: `gallope.out.nc` reads cleanly with 1829 records
covering `t` in `[0, 920]`, well past the `t = 900` end of the averaging
window, so the hang costs nothing the paper depends on.  It was not chased
further — reproducing it would cost 19+ h of GPU time on a machine that shuts
down on 2026-08-13, and neither the correctness nor the performance claim
rests on running past `t = 920`.

Section 5 quantities, time-averaged over `t` in `[200, 900]`
(`diagnostics/MHD_INCOMP/time_average.txt`, indices 393–1787):

| | Gallope | Calliope |
|---|---|---|
| `<W>` | 4.043 +- 1.150 | 3.792 +- 0.808 |
| `<W_dot>` | 5.055e-3 +- 1.724e-1 | 1.236e-3 +- 1.334e-1 |
| `<P>` | 1.987 +- 0.562 | 1.945 +- 0.404 |
| `<D>` | 1.980 +- 0.574 | 1.943 +- 0.409 |

`<W_dot>` is zero to within its scatter and `<P>` matches `<D>` to 0.4 %
(Calliope 0.1 %), so the window is genuinely in the statistically steady state
and the energy budget closes.  The two codes agree to 6.6 % in `<W>`, 2.2 % in
`<P>` and 1.9 % in `<D>`, all inside 1 sigma.  Averaging both runs over the
*same* window is a change from the first revision, which used `[200, 676]` for
Gallope and `[200, 573]` for Calliope.

## What is in this directory

| | |
|---|---|
| kept, tracked | `gallope.in`, `job.pbs-prod-gallope`, `out.std`, this file |
| kept, untracked | `rms.dat`, `cfl.dat`, `gallope.out*.nc`, `out2d/`, the figures under `diagnostics/` |
| not copied from Miyabi | `restart/` (3.1 GB), `gallope.modes.out` (2.1 GB), and all but the last snapshot of `out2d/` |

`out2d/` here is **only the final snapshot**, 48 MB out of the 10 GB on
Miyabi.  Fig. 4(b) shows the three orthogonal cuts at the last output time and
nothing else does, so each of the 36 `{field}_r_{x0,y0,z0}.dat` files was
truncated to its last record with `tail -c` (2 097 152 B for `z0`,
1 048 576 B for `x0` and `y0`) and `time.dat` to its last line.
`plot_fields.py` reads the result unchanged and reproduces the panel at
`t = 920`.  There is no `out3d/`: the input sets `write_intvl_3D = 0.0d0`.

`diagnostics/` is a copy of the analysis tree as it stood when the figures
were made, so the archive reproduces them without depending on the current
state of `diagnostics/` at the repository root.  Note that
`diagnostics/MHD_INCOMP/load.py` opens `<runname>.out.nc` unconditionally, so
the NetCDF file — untracked, but present here — is what makes this archive
self-contained.  `time_average.py` here has `avg_start_time`/`avg_end_time`
pinned to the window above rather than the repository defaults.

Verified self-contained: with Miyabi unreachable, `time_average.py` regenerates
`time_average.txt`, `fig_energy/` and `fig_kspectrum/` (Fig. 4(a) and (c)) and
`plot_fields.py` regenerates `fig_fields/b.pdf` (Fig. 4(b)) from what is in
this directory alone.  The `*_r`/`grid` MATLAB dumps `plot_fields.py` writes
alongside the PDFs are intermediates for external rendering and were deleted;
rerunning the script recreates them.

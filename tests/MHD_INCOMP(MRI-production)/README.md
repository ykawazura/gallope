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
otherwise destroy the first leg's time history.  If `legs/` is absent from this
archive, the run finished inside a single leg.

## What is in this directory

| | |
|---|---|
| kept, tracked | `gallope.in`, `job.pbs-prod-gallope`, `out.std`, this file |
| kept, untracked | `rms.dat`, `cfl.dat`, `gallope.out*.nc`, the figures under `diagnostics/` |
| not copied from Miyabi | `restart/` (3.1 GB), `out2d/` (8.0 GB), `gallope.modes.out` (2.1 GB) |

`diagnostics/` is a copy of the analysis tree as it stood when the figures
were made, so the archive reproduces them without depending on the current
state of `diagnostics/` at the repository root.  Note that
`diagnostics/MHD_INCOMP/load.py` opens `<runname>.out.nc` unconditionally, so
the NetCDF file — untracked, but present here — is what makes this archive
self-contained.

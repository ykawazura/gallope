# KRMHD-extension regression check — RMHD & MHD_INCOMP

**Verdict: PASS.** The Phase-1 KRMHD extension (which modified the shared
`mp` / `mpiio` / `grid` / `force_common` layers plus the RMHD model files) does
**not** break the existing `RMHD` and `MHD_INCOMP` models. One real compile-time
regression was found and fixed along the way (see below).

Date: 2026-07-19. Machine: Miyabi-G (JCAHPC), 4 GPUs (`select=4:mpiprocs=1`,
`debug-g`). Parallelization: default → `P_m = P_s = 1`, `P_fft = np = 4`
(exercises the distributed reduction / FFT / MPI-IO paths that the extension touched).

## Runs

| # | model | input | grid | nstep | dt₀ | wall | steps run | completion |
|---|-------|-------|------|-------|-----|------|-----------|------------|
| 1 | RMHD              | `input_example/RMHD.in`               | 64³      | 5000  | 5.0e-3 | 24 s | 5000/5000   | `Program completed! ^_^` |
| 2 | MHD_INCOMP forcing| `input_example/MHD_INCOMP(forcing).in`| 64³      | 10000 | 2.0e-3 | 69 s | 10000/10000 | `Program completed! ^_^` |
| 3 | MHD_INCOMP MRI    | `input_example/MHD_INCOMP(MRI).in`    | 64²×32   | 10000 | 1.0e-1 | 65 s | 10000/10000 | `Program completed! ^_^` |

All three jobs exited `rc=0` (`REGRESS_DONE rc=0`) and ran the full `nstep`.

## Judgment criteria & results

A pre-KRMHD baseline does not exist in this repo, and all three inputs are
driven with random/clock-seeded ICs (non-reproducible), so bit comparison is
impossible. The regression is therefore judged on three criteria:

1. **Clean compile** — ✅ after fixing the regression below.
2. **Run to completion, no NaN/Inf, no blow-up** — ✅ `np.isfinite` is `True`
   over *every* row of all three `rms.dat` files; total energy stays bounded.
3. **Physical validity** — ✅ (see `energy_timeseries.png`):
   - **RMHD (driven):** energy grows and builds toward a turbulent state,
     `E_tot ∈ [0.006, 14]`, bounded; the step transitions near `t≈13.5, 16`
     are `dt`-reset events, not divergence.
   - **MHD_INCOMP forcing:** starts from `u=0`, forcing grows turbulence which
     saturates (`E_b≈1.7`, `E_u≈1.1`), `E_tot ∈ [0.003, 2.9]`.
   - **MHD_INCOMP MRI:** `shear=.true., q=1.5` (Keplerian). Exponential growth
     followed by recurrent bursts — the MRI channel-mode ↔ turbulence cycle —
     `E_tot ∈ [0.5, 89]`, bounded (no runaway).

### Note on `rms.dat` sampling
`get_nonlinear_terms` (which writes the `rms` line) is called once per
**eSSPIFRK3 sub-stage**, so `rms.dat` holds **3 rows per timestep** sharing the
same `tt`. The plot collapses each step's 3 sub-stage rows to one point (mean
over stages) and overlays a moving-average trend; the raw per-step points are
drawn thin/faint. Regenerate with `python3 plot_regress.py`.

## Regression found & fixed (compile-time)

**Symptom:** `MODEL=MHD_INCOMP` failed to build —
`NVFORTRAN-S-0155 Could not resolve generic procedure mpiio_read_one`
(`fields.F90:857-862`, 6 severes).

**Root cause:** the KRMHD extension added a **required (non-optional) `comm`
argument** to the shared generic `mpiio_read_one` / `mpiio_write_one`. The RMHD
calls were updated to pass `comm_fft`, but the **12 MHD_INCOMP calls** (6 reads
in `fields.F90`, 6 writes in `io.F90`) were not — i.e. the refactor genuinely
broke MHD_INCOMP. This is exactly the class of regression this check exists to catch.

**Fix (working tree, uncommitted):** `src/model/MHD_INCOMP/{fields.F90, io.F90}`
— pass `comm_fft` to all 12 calls and import it via `use mp`, matching the RMHD
pattern. MHD_INCOMP always runs `P_m = P_s = 1`, so `comm_fft == world`,
consistent with the `proc_id` offset — the physically correct communicator.
`mpiio_write_var` / `_var_2d` (fh/disp based) take no `comm` and were untouched.

## Files

Tracked (per `.gitignore`): each run's `gallope.in`, `job.pbs-regress`,
`out.std`, and this `README.md`.
Untracked/local (reproducible): `rms.dat`, `cfl.dat`, `gallope.modes.out`, the
PBS `regress.o<jobid>` logs, `energy_timeseries.png`, `plot_regress.py`.

See `tasks/todo_regression.md` for the full worklog and `tasks/lessons.md` (L21)
for the lesson recorded from this finding.

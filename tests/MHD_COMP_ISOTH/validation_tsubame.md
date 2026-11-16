# MHD_COMP_ISOTH validation on TSUBAME4.0

**Verdict: all 589 cases ran and were judged. Every check gives the verdict of Calliope's reference except
three, which fail where the reference passed: Euler1D `collide_U2_N512_est_chain` (a check of Calliope's
Stage 1 version), MHD1D `suite/C2_N4096` against s13d (last line 1.27e-11 > 1e-11 after 40,760 steps) and
Smooth_wave `2d_k102_a1_N32_default_p400` (growth of a parasitic mode seeded by rounding, γ/ω0 0.116 > 0.105).
The failures of Calliope's own runs are reproduced, and the 14 aborted runs are those of Calliope. Four long
Turb3D runs leave the s13d trajectory (listed, not judged).**

Date: 2026-11-14 to 2026-11-16. Machine: TSUBAME4.0 (H100 SXM5, 4 per node). Gallope `feature/isoth`, built
with `nvhpc/26.1_cuda13.1`, `openmpi/5.0.2-nvhpc`, `hdf5-parallel/1.14.3`, `netcdf-parallel/4.9.2`,
`MODEL = MHD_COMP_ISOTH`. Calliope: `develop` 1d924a9 for the inputs, the initial conditions and the analysis
scripts (its `src/model/MHD_COMP_ISOTH` is the source of the s13d executable, 07716cc). The Gallope model is a
port of that source.

## What was run

The 589 cases of Calliope's `tests/MHD_COMP_ISOTH`, the runs of the paper, ran with
`bin/gallope-MHD_COMP_ISOTH`:

| folder | cases | Calliope executable of the reference runs |
|---|---|---|
| `Euler1D` | 102 | `calliope_isoth_es` (Stage 1, model 2126274) |
| `MHD1D` | 215 | s13d (final, model 07716cc) |
| `Smooth_wave` | 167 | `calliope_isoth_s4` (Stage 4) |
| `OT2D` | 15 | s13d |
| `Turb3D` | 90 | s4, s4b, s6 (Stages 4–8, 20 cases), s13d (Stages 13–16, 70 cases) |

Every case was judged by the analysis scripts of Calliope with Calliope's criteria, unchanged. Two kinds of
Calliope's checks compare two Calliope versions and have no Gallope counterpart: D2, D3, D5, O3, O4 and CV2
of `REPORT_stage14`/`REPORT_stage16`, and the pairs of runs with `fftw_estimate`. They are replaced by a
direct comparison of every Gallope run with the s13d output of the same case (`compare_s13d_all.py`), where
s13d ran the case. Where Calliope ran an older version (`Euler1D`, `Smooth_wave`, `Turb3D` Stages 4–8),
the Gallope results are judged by the physics criteria only and differences from the Calliope values also
contain the changes between the Calliope versions.

What differs from a Calliope run is listed in `README.md` ("What differs from a Calliope run"). The
differences that change a check:
- MPI decomposition: Gallope needs P_fft to divide nlx and nly, so a 1D or 2D case (nly = 2) runs on 1 or
  2 GPUs. Calliope's decomposition tests on 2 to 256 ranks become 1-GPU against 2-GPU runs
  (`ngpu_map.md`).
- `pruned = T` (the 3/2 rule) does not exist in Gallope; the two `pruned` cases run the 2/3 rule.

## How to reproduce

The commands are in `README.md` (Setup, Running, Analysis). In short, on TSUBAME4.0:

```bash
$PY $MC setup --calliope $CALLIOPE_ROOT --runroot $RUNROOT
$PY $MC ic    --calliope $CALLIOPE_ROOT --runroot $RUNROOT --python $PY
$PY $MC check --runroot $RUNROOT
cd $RUNROOT
$PY $MC list --runroot . --folder <folder> --status todo --max-ngpu 1 > list.txt    # split into node_q lists
CASEFILE=list.txt RESOURCE=node_q H_RT=4:30:00 PYTHON=$PY $GALLOPE/tests/MHD_COMP_ISOTH/job-tsubame.sh
FOLDERS="Euler1D MHD1D Smooth_wave OT2D Turb3D" CALLIOPE_ROOT=$CALLIOPE_ROOT PYTHON=$PY REF=$REF \
    $GALLOPE/tests/MHD_COMP_ISOTH/judge-tsubame.sh
```

- Inputs: Calliope's `calliope.in` of each case, translated by the job (`calliope_run.in`, `README.md`).
  Initial conditions: Calliope's `stages/regen_inputs.py` with numpy 2.0.2 on TSUBAME; the md5 of all 487
  generated files equals Calliope's `generated.md5` (the other 102 cases have built-in initial conditions).
- Job scripts: `job-tsubame.sh` (runs), `judge-tsubame.sh` (judgments on `cpu_16`).
- `REF`: the s13d outputs relayed from ps1 without `out3d/`, NetCDF files and executables
  (`/gs/bs/tgh-26IDS/ua07691/s13d`).

## Direct comparison with s13d (`compare_s13d.py`)

Three outputs decide the verdict: t of the final restart (relative 1e-12), the seven restart fields
(relative L2 ≤ 1e-11) and the out1d lines (Calliope's Stage 14 measure `analyze_stage3.maxdiff`, ≤ 2e-11
on every line, ≤ 1e-11 on the last line of `suite` and `eta_suite`). The rows of `cfl.dat`,
`conservation.dat`, `budget.dat` and the records of `out2d/plane_{y0,z0}.dat` must pair by t over the
restart pieces of both runs; their differences, max|a − b|/max|b| per column, are listed and not judged.
Columns below 1e-13 in both codes are taken as rounding level and not listed.

Every listed column above 1e-6 of MHD1D and OT2D was checked; none comes from a difference of the solutions:
- `conservation.dat` 4, (M − M0)/M0: M0 is M at the start of each piece. Where Calliope split a run and
  Gallope did not (2 to 6 pieces), the column restarts from 0 in one code only: 1.0 to 4.6
  (`stage16/conv/RDpc_N8192`, 6 pieces). Within one piece the column is 1e-13 to 1e-11 and differs by the
  rounding of M (≤ 1e-14, e.g. `suite/TDS_w4_N2048`: 1.4e-13, difference 1.0e-14); M itself agrees to 1e-14.
- `budget.dat` 3, 7, 9 (dM/dt and dE/dt from the change over one step, and their error err_t): the
  rounding of M or E divided by dt, 5e-11 to 4e-10 for dt = 2e-5 to 2e-4. Where the column is itself of
  that size, the difference is of the order of the column. Gallope's values are again the smaller
  (`suite/TDS_w4_N2048` dM/dt: 2.2e-11 against 3.0e-10).
- `budget.dat` 24, 25, 29, 30 (R_svv, R_eta and their ratios to the dissipation): residuals of identities,
  R_svv ≤ 3.7e-14 and R_eta ≤ 1.7e-16; the ratios exceed 1e-13 (up to 1.5e-12) because they are divided by
  the dissipation.
- `conservation.dat` 7 in OT2D (P_z, zero by symmetry): ≤ 2.8e-14 in Gallope, ≤ 9.0e-13 in Calliope.
- `cfl.dat` 4 and 13–16 (dt_visc, λ_svv of the power method and its iterations): the power method stops at
  a relative change of 1e-3, and its vector is not in the restart files, so a run that restarts starts a
  new vector. After a seam of Calliope λ_svv differs by up to 2.6e-2 (`stage16/conv/RDpc_N4096`) and the
  iteration counts differ; in one piece up to 1.8e-3 (`suite/RDp_w0_N2048`) with equal counts. λ_svv
  enters only the bound dt_visc of dt; dt (column 2) is equal in every paired row of every case.

## Gates before the full runs

| gate | what | criterion | result |
|---|---|---|---|
| G0 regression (job 9386336) | the six existing models build; R1 runs of `tests/release` (MHD_INCOMP, MHD_INCOMP_SHEAR, RMHD; 1 and 4 GPUs) against v1.1 | restart files bitwise equal | 36 files, 0 differ |
| G1 FFT with a size-2 dimension (job 9386134) | cuFFTMp round trip and mode mapping, shapes 2×8×8, 8×2×8, 8×8×2, 6×6×6, 2×2×8, 2×8×2, 8×2×2, 16×32×64 on 1 and 2 GPUs | about 1e-15 | 16/16 PASS, worst 2.46e-15 |
| G2 input | `make_cases.py check`: every translated namelist group and key exists in Gallope | all 589 | PASS; initial conditions 487 md5 equal |
| G3 MHD1D | all 215 cases against s13d (below) | D2/D3 tolerances | 214/215 (below) |
| G3 OT2D, Turb3D (job 9398327) | `dec_N256_*_n10` (10 steps), `stage15/pre` (10 steps) against s13d | relative 1e-12 | 8/8 PASS: restart L2 ≤ 3.2e-16, planes ≤ 1.1e-12, t identical |

## Results

### Euler1D (`REPORT_stage1`; `judge_stage1.py`)

| check | cases | criterion (Calliope's) | Calliope (es) | Gallope |
|---|---|---|---|---|
| shock | 20 | ρ*, shock speed, positivity, conservation per case; L1(ρ) order fitted for N ≥ 512 | 20 pass; 4 order fits pass | 20 pass; 4 order fits pass (0.921–0.985); ρ_min, L1(ρ), ρ* err and the orders equal to the digits shown |
| scan | 31 | as shock; 8 rows expected to fail (Sod C0/C1, U2 C0/C1/nosvv, U4 C1/linear/nosvv) | 23 pass, 8 FAIL | 23 pass, the same 8 FAIL |
| smooth | 19 | sound-wave amplitude and phase, simple-wave convergence (N → 2N), D_svv ≥ 0, R_svv/A_svv ≤ 1e-10 | 13 pass, 3 reference, 3 without a verdict (N = 512) | the same; Linf(s) equal to 4 digits |
| dtconv | 5 | mass error order ≥ 2.7; `U4_N2048_cflsvv1` vs `shock/U4_N2048` ≤ 1e-3 | pass | pass; the table equal to the digits shown (order 3.03; Δ(ρ*/ρ*_exact) 1.07e-7) |
| decomp | 9 | max\|Δρ\|, max\|Δu\| vs `x_np1` ≤ 1e-10 | ≤ 2.93e-14 | x, y: 0; z: 1.15e-14 (ρ), 1.67e-14 (u) |
| smoke `selftest_np{1,2,4}` | 3 | ζ bitwise equal to np1 | 0 | 0 on 1, 2, 4 GPUs |
| smoke `nan`, `sound64`, `collide_U2_N512{,_est}` | 4 | stop on a non-finite value; conservation | pass | pass |
| smoke `collide_U2_N512_est_chain` | 1 | out1d equal to the unchained run; 1001 rows of `conservation.dat`; FFTW_ESTIMATE in both | pass | **FAIL**: out1d diff 0, 1002 rows, no FFTW_ESTIMATE (below) |
| smoke `collide_U2_N512_C1_nosvv_fs` | 1 | reference (no verdict) | stopped by the CFL check, t_end 0.0158, ρ_min 1.38e-3 | the same t_end and ρ_min; exit label `dt_stop` |
| cflmax | 6 | reference rows | f2.0, f2.4 exit 1 (dt·λ_svv check), f3.0 stopped by the CFL check | the same exits, t_end and dt·λ_svv to the digits shown |
| fast | 3 | fast eigenmode, ω against √(2 + √2) and the amplitude at t = 10 (no threshold; `cf9bb60:fast_wave.py`) | +2.96e-7, 0.999961 (ideal); −3.99e-7, 0.998010 (default, default_linear) | the same to the digits shown |

- `est_chain` checks that a run split by a restart repeats the unsplit run under FFTW_ESTIMATE. Gallope has no
  FFTW plans, and its `conservation.dat` repeats the row of the seam (t = 1e-4): the first run writes it when
  it ends and the restart writes it again. Calliope's s13d does the same (added in 3b784b1 after Stage 1;
  `stage16/conv/I1_N8192` of s13d repeats its seam rows). Without the repeated row, all columns but
  (M − M0)/M0 are bitwise equal to the unchained run. The FAIL is version-specific.
- In the decomposition test, `y` is a direction of the c2c transforms in Gallope and `z` the r2c direction;
  in Calliope `y` is the r2c direction. Both codes differ from `x_np1` only in the r2c orientation. Within
  each orientation the 1-GPU and 2-GPU runs are bitwise equal.
- The conservation residuals dMdt_raw, dP and R_ζ are about ten times smaller than Calliope's (for example
  Sod_N1024: 1.21e-16, 3.16e-17, 1.03e-15 against 9.87e-16, 1.37e-15, 1.31e-14), all far below the
  criteria.
- The two `pruned` cases repeat `*_N1024_fs_nt6` to all digits; both pass. Calliope's comparison of the 3/2 and
  2/3 rules (`REPORT_stage1` section 9) is not reproduced.

### MHD1D (`REPORT_stage14`, `REPORT_stage16`)

| check | cases | criterion | Calliope (s13d) | Gallope |
|---|---|---|---|---|
| D1 complete | 207 | every case reaches t_end with all out1d lines (`complete_1d`) | 207 | 207 |
| CV1 complete | 8 (`stage16/conv`) | as D1 | 8 | 8 |
| D4 decomposition | 15 + 10 | max diff vs `_x_np1` ≤ 1e-13, equal line steps | 1.15e-14, 6.86e-15 | 8.58e-15 (`F1_N256_z_np4`), 5.08e-15 (`I1_N256_eta2_z_np4`) |
| s13d: every out1d line | 215 | max\|a − b\|/max(1, max\|b\|) of the 6 variables ≤ 2e-11 (D3/CV2 tolerance) | – | 215/215; worst 1.80e-11 (`stage16/conv/RDpc_N8192`) |
| s13d: last line | 182 (`suite`, `eta_suite`) | ≤ 1e-11 (D2 tolerance) | – | 181/182; **FAIL** `suite/C2_N4096` 1.27e-11 |
| s13d: restart, t | 215 | – | – | restart L2 ≤ 4.3e-13; t equal in all |

- The difference from s13d grows with the step count. It exceeds 1e-12 in 39 cases and 1e-11 in 3, all with
  N = 4096–8192 and 4e4–1e5 steps. `suite/C2_N4096` (40,760 steps) passes the line tolerance 2e-11 (1.40e-11)
  but not the last-line tolerance 1e-11. The D2 tolerance was set for two CPU builds of Calliope (its worst:
  3.49e-12 `suite`, 5.11e-12 `eta_suite`); a GPU build sums in another order. The case is kept as a FAIL.
- The tolerances of this comparison were changed after its first results. The first measure, every column
  relative to its own maximum, failed 67 of 215 cases in columns at the rounding level (u_x ~ 1e-6 in the
  TDS problems, div u). It was replaced by Calliope's own measure of Stage 14 (`analyze_stage3.maxdiff`).

### Smooth_wave (`REPORT_stage4`; `analyze_stage4.py alfven|cpdiag|planar`)

| check | cases | criterion | Calliope (s4) | Gallope |
|---|---|---|---|---|
| alfven | 63 | error, symmetry, permutation, decomposition, dt³ order; growth of a parasitic mode γ ≤ 1.25 × the physical maximum for N ≥ 32 with the default dissipation | 1 FAIL (`3d_k112_a1_N32_default_p400_np8`, γ/ω0 0.132) | **2 FAIL**: `3d_k112_a1_N32_default_p400_np8` (0.124) and `2d_k102_a1_N32_default_p400` (0.116 > 1.25 × 0.084) |
| alfven `mpi` | (3 of the 63) | diff from the 1-GPU run | 0 | 0 (2 and 1 GPU) |
| cpdiag | 16 | all runs finish (no verdict) | finished | finished |
| planar | 88 | conservation, symmetry, equivalence of 1D and 2D, convergence | PASS | PASS (every check) |

- The second FAIL is the growth of a mode seeded by rounding errors and driven by the shock resistivity
  (Calliope `REPORT_stage4` 5.4). Calliope measured γ/ω0 = 0.092 on another dominant mode, (2, 0, −4);
  Gallope's dominant mode is (5, 0, 1). Rounding errors of another code seed other modes. The criterion is not
  relaxed.
- cpdiag γ/γ_galerkin, Gallope/Calliope: 2D N16 etaonly 6.06/7.04, default 5.92/6.12; 3D N16 etaonly
  5.64/5.17, default 5.77/4.84; 3D N32 etaonly 1.00/1.00, default 1.30/1.39. For 2D N16 `noeta` and `svvonly`
  the fit has 2 records above 10 A(T) in Gallope and 3 in Calliope (its minimum), so Gallope gives no γ.
- `2d_k102_a1_N64_ideal_p400` stops with `dt_stop`, which `analyze_stage4.py` allows (line 304).

### OT2D (`REPORT_stage14`)

| check | cases | criterion | Calliope (s13d) | Gallope |
|---|---|---|---|---|
| O1 complete, budget | 7 (N = 128–512 with C_η = 0 and 2, N = 1024 with C_η = 0) | complete; budget residual R ≤ TOL_RES = 1e-10 | 7; R ≤ 5.1e-14 | 7; R ≤ 4.1e-14 (`ot_N1024_eta0_p256`) |
| O2 decomposition | 8 (`dec_N256_*`, `*_n10`) | planes after 10 steps vs `np16` ≤ 1e-12 | ≤ 7.8e-16 (4, 64, 256 ranks vs 16) | 0 (2 GPUs vs 1); planes and conservation.dat up to t = 1 also 0 |
| O3 convergence | C_η = 0 (N 128–1024), C_η = 2 (N 128–512) | verdict and failing items equal to the reference | both FAIL: 2 items (C_η = 0), 9 items (C_η = 2) | both FAIL with the same items |
| O4 t ≤ π | 7 | relative difference ≤ 1e-6 in conservation.dat, the planes ρ, u, B and the S4.3 scalars | (old vs s13d) cons ≤ 1.5e-13, planes ≤ 1.1e-12, scalars ≤ 1.3e-7 | (vs s13d) cons ≤ 2.0e-12, planes ≤ 7.4e-12, scalars ≤ 1.9e-11; ζ area fractions equal |
| s13d: restart, t, rows | 15 | as for MHD1D | – | 15/15; restart L2 ≤ 7.3e-12, t equal; planes over t ≤ 2π ≤ 1.0e-10 (listed) |

- O3 and O4 are Calliope's own checks between two versions; here `analyze_stage14.py ot` compares Gallope
  with s13d, whose runs are linked in place of Calliope's "old" runs (`judge/ot_s13d`). The reference of O3
  is therefore s13d analyzed by the same scripts.
- `analyze_stage4.py ot` ends with `ot: FAIL`: the convergence verdicts fail, in s13d with the same items
  (Calliope `REPORT_stage14` 6.2: plasmoids form at t ≥ 4 and the pointwise values do not converge in N).
- After t = π the planes depart from s13d as they departed between two Calliope versions (`REPORT_stage14`
  6.4). At t ≈ 2π, Gallope vs s13d: 2e-13, 1e-11, 7e-11 (N = 128, 256, 512, C_η = 0) and 9e-14, 8e-13, 6e-11
  (C_η = 2); Calliope's old vs new: 1.3e-13, 3.1e-12, 2.2e-10 and 4.2e-14, 1.3e-12, 6.7e-12. max|J| and the
  ρ range at t = 3, 5, 6 equal Calliope's to the digits printed (e.g. N = 512, C_η = 0, t = 5: max|J| 143.2,
  ρ 0.195–7.157).
- Calliope ran N = 512 in 2 pieces and N = 1024 in 3; Gallope ran each case in one piece (1 GPU).

### Turb3D (`REPORT_stage4`–`REPORT_stage16`)

All 90 cases ran on 1 GPU except the decomposition cases (`ngpu_map.md`). τ0 is the eddy turnover time of
Calliope's scripts.

| check | cases | criterion | Calliope | Gallope |
|---|---|---|---|---|
| Stage 4 turb: runs | 5 (`stage4/turb`) | complete; budget residual R ≤ TOL_RES = 1e-10 | turb1 64³, 128³ pass; `turb4_N64` FAIL (`dt_stop`); `turb4_N64_a`, `turb4b_N64` cut at dt ≈ 2e-8 | the same; `turb4_N64` `dt_stop` at t = 0.1016; R ≤ 1.7e-14 |
| Stage 4 turb: ε(t) peak | turb1 64³, 128³ | reference | 35.3 at 0.755 τ0; 36.5 at 1.058 τ0 | 35.331 at 0.755 τ0; 36.461 at 1.058 τ0 |
| Stage 4 turb: convergence | turb1 64³ vs 128³ | at t_end ρ_min, ρ_max, E_K, E_M ≤ 2%; ε peak ≤ 5%; its t ≤ 2%; share of the shock viscosity ≤ 10% | FAIL: 13, 2.8, 1.7, 1.1; 3.1; 29; 15% | FAIL: 13.0, 2.81, 1.65, 1.08; 3.10; 28.6; 15.2% |
| restart seams (`check_seams.py`) | 5 cases in 2–6 pieces | E budget summed across the seams | PASS, residual 0 | PASS, residual 0 |
| Stage 5 res | `turb4_N128_a`, `turb4b_N128` | "collapse at grid scale": the hole deepens one cell thick at 128³ | collapse | collapse; the widths of Calliope's table at every t listed, e.g. thin axis 0.90 and 1.14 cells at t = 0.11, 0.21 |
| Stage 5 ms | M_s = 1, 2, 3, 4 at 64³ | S: 3τ0 with dt ≥ dt0/10; H: ρ_min < 1e-3 and dt < dt0/10 | S, not S, not S, H | the same; dt < dt0/10 from 0.477 τ0 (M_s = 2) and 0.245 τ0 (3); H at 0.234 τ0 (4) |
| Stage 6 track: decomposition | `dec_N64_np{1,8,64}`, 400 steps | dt column bitwise equal to the reference; E_K, E_M, ρ_min ≤ 1e-10 of the column maximum | dt equal; ≤ 5.4e-15 (P = 1, 8 vs 64) | dt equal; ≤ 4.6e-16 (2 and 4 GPUs vs 1) |
| Stage 6 track: vs Stage 4 | `turb1_N{64,128}_t` to 3τ0 | d_track ≤ d_res/10 for 14 quantities | PASS; 5,017 and 12,052 steps; at most 0.77% of the limit | PASS; 5,017 and 12,052 steps; at most 1.6% of the limit (ρ_max, 128³, 3τ0) |
| Stage 6 turb | `turb1_N256` | as Stage 4 turb; convergence 128³ vs 256³ | 40,064 steps; ρ_min 0.0090; ε peak 38.5 at 1.121 τ0; convergence FAIL (ε 5.2%, t 5.6%, share 14%, ρ_min 12%, ρ_max 11%, E_K 0.50%, E_M 0.78%) | 39,780 steps (below); ρ_min 0.0090; ε peak 38.478 at 1.121 τ0; R 5.2e-15; convergence FAIL (ε 5.24%, t 5.63%, share 14.3%, ρ_min 12.0%, ρ_max 11.2%, E_K 0.506%, E_M 0.783%) |
| Stage 7 tail, hole | turb1 64³, 128³, 256³ | low-density tail of 128³ vs 256³; share of holes not explained by ideal MHD | σ_s(1τ0) 0.366, 0.387, 0.399; ρ_min 0.114, 0.093, 0.042; "intermediate"; 12 of 36 holes suspect | σ_s 0.3658, 0.3869, 0.3992; ρ_min 0.1137, 0.0926, 0.0419; tail 128³–256³ at each t as Calliope's table (F_V ratios to 3 digits; n_fine 65 against 64 at 3τ0); "intermediate"; 1 of 6 holes suspect at 256³ (below) |
| Stage 8 coefficients | 5 (`svv_eps0` and `C_shock` doubled and halved, 128³) | noise of the control ≤ 0.010; an axis is sensitive if \|R\| ≥ 0.1 | noise −3.2e-7; R −0.0830, +0.0218 (svv_eps0), +0.0688, −0.0785 (C_shock); both insensitive | noise 0; R equal to the 4 digits; both insensitive |
| Stage 13 nomass | `svv_mass = F` and control, 128³ | label from L of N and C over 0.6–1.2 τ0 | N unstable at 0.466 τ0, no label; d ln dt/d ln ρ_min 1.4996 | the same: N at 0.466099 τ0, no label; 1.500; L of C −2.3746 (Calliope −2.3747) |
| Stage 15 pre | 12 | dt halving: residual ratio ≈ 8; P = 8, 64 planes bitwise equal; ε_ref | 7.98–8.00; equal; ε_ref 0.1424333 (0.7555 τ0), 0.01471188 (0.8487 τ0) | 7.98–8.00; equal (2 GPUs vs 1); 0.1424333 (0.7555 τ0), 0.01471188 (0.8487 τ0) |
| Stage 15 a, b, c | 30 + 14 + 2 | r*_L, r*_H at 64³; r_C at 128³; 256³ runs | r* (M_s = 1) L ≤ 1/4, H 2; r_C = H(8); H(8)-256 and H(32)-256 stopped | the same; the 8 aborted runs stop at the t/τ0 of Calliope (3 digits) |
| Stage 15 judge | – | prior rule on κ and cost | B* = H(8); undetermined; cost ratio 3.01 > 2.44 | B* = H(8); undetermined (M_s = 1, flagged: ES-256 is the only reference; M_s = 0.5 as Calliope's to 3 digits); κ below; cost ratio 3.31 > 2.69, reported only |
| Stage 16 hf | 7 + 3 + 1 | as Stage 15 with the filter of ln ρ off | r*_HF = 2; r_C,F = 8; HF(8)-256 stopped at 0.577 τ0 | the same; the 6 aborted runs stop at the t/τ0 of Calliope (3 digits) |
| Stage 16 rejudge, judge1 | – | Stage 15 judge with ES-256 as the only reference; post-hoc | "ES superior" post hoc; cost rule FAIL | "ES superior" in variants A and B; M_s = 0.5 undetermined as Calliope's to 3 digits; judge1: κ of HF(8) equal to that of H(8) to 3 digits (Calliope: within 0.007), B* and verdict unchanged; cost ratio 3.31 > 2.57, reported only |
| s13d: restart, t, rows | 70 (Stages 13–16) | listed, not judged (below) | – | 52 of 56 compared cases within the tolerances of `compare_s13d.py`, including `stage15/pre/dec_*` (10 steps, judged) 4/4; 14 aborted cases not compared; 4 long runs differ (below) |

- Calliope cut 9 cases before their nominal t_end: `turb4_N64_a`, `turb4b_N64`, `turb4_N128_a`,
  `turb4b_N128` (dt fell to 2e-8–4e-7 in the holes), `turb2b_N64`, `turb3b_N64` (dt < dt0/10, cut by the user),
  `stage13d/nomass/e0.5_C2` (stopped by the budget) and the two `ESf` runs of Stage 15 (stopped by the
  cost). Gallope runs these cases to t of Calliope's last restart (file `t_end`; the nominal value is kept in
  `t_end.nominal`), so they end with `exit 0` where Calliope's report says incomplete or stopped. The
  verdicts are unaffected: each check reads the run up to its last row, as for Calliope's cut runs.
- The 14 `exit 1` cases are the runs that Calliope reported as aborted (Stage 15: R0-64, H(1/4), H(1/2),
  H(1) at 64³, H(2), H(4) at 128³, H(8), H(32) at 256³; Stage 16: HF(1/4), HF(1/2), HF(1) at 64³, HF(2),
  HF(4) at 128³, HF(8) at 256³). Each stops at the t/τ0 of Calliope to 3 digits. ρ_min at the stop of the
  8 Stage 15 runs equals the value of `REPORT_stage15` to 3 digits in 7 (e.g. R0-64: 2.98e-21 at 0.512 τ0;
  H(8)-256: 2.36e-8 at 0.614 τ0) and differs by 1.4% in H(2)-128 (3.59e-4 against 3.64e-4).
- The decomposition cases `stage6q/track/dec_N64_np*` end with `dt_changed`, as in Calliope: 'track'
  changes dt and the runs have no `t_end`.
- The decomposition check of `analyze_stage4.py turb` (`stage4/turb/dec_N64_np{1,8}`) is empty: these
  cases are not among the 90.
- The core-h and cost columns of `analyze_stage6.py`, `analyze_stage15.py` and `analyze_stage16.py` are
  Gallope's wall time multiplied by Calliope's `nproc`, which means nothing for a GPU run. The cost rule of
  Stage 15/16 is therefore reported, not judged. Gallope's solve time: 0.0263 s/step at 64³ and 0.0553 s/step
  at 128³ and 0.401 s/step at 256³, all on 1 GPU.
- `turb1_N256` took 39,780 steps against Calliope's 40,064 (−0.7%); 'track' changed dt 82 times. Calliope ran
  Stage 6 with an older executable, which has no s13d output. 'track' sets dt to 0.8 dt_cfl rounded down to
  3 significant digits, so a rounding-level difference of dt_cfl at a rounding point can move a change of dt
  by one step (seen below in `stage15/b/turb1_N128_ES`). Whether this is the cause at 256³ was not checked.
- Stage 7 analyses the holes of 256³ in every out3d record. Calliope's chain wrote a record at the end of each
  job (0.578–0.965 τ0, `REPORT_stage7` §3); Gallope ran `turb1_N256` in one job and has the records of the
  output schedule only, so the 7 extra times print "no record". The 256³ set therefore has 6 holes at
  0.5 and 1τ0 against Calliope's 36 at 9 times. The minimum points at these two times are Calliope's:
  (134, 55, 147), ρ 0.1523, "balance" at 0.5τ0; (69, 18, 6), ρ 0.0419, suspect at 1τ0. Suspect holes are
  1 of 6 (Calliope 12 of 36), below half in both, so the next step is the same (step 3). 64³ and 128³:
  2 of 2 and 1 of 2 suspect, as Calliope.
- κ of Stage 15 (M_s = 1, 128³, pairs E_B W1, E_B W2, E_m W1, E_m W2): ES 0.706, 0.742, 0.739, 0.445 and
  H(8) 0.455, 0.535, 0.474, 0.555 (Calliope 0.713, 0.743, 0.744, 0.436 and 0.477, 0.534, 0.496, 0.557).
  Stage 16, variant A: ES 0.715, 0.740, 0.753, 0.433 and H(8) 0.471, 0.533, 0.481, 0.557 (Calliope 0.713,
  0.742, 0.746, 0.429 and 0.479, 0.538, 0.500, 0.558); δ 0.047, 0.047, 0.047, 0.596 (Calliope 0.052, 0.047,
  0.047, 0.449). The 128³ runs equal s13d (`turb1_N128_H_r8`: restart 1.3e-13), so these differences of
  κ come from the reference ES-256, which is each code's own `turb1_N256`. The M_s = 0.5 judge uses
  (64, 128) only and its κ, δ and D equal Calliope's to 3 digits.
- Variant A of Stage 16 averages κ over the out3d records in each window. Calliope's chained 128³ runs wrote a
  record at the end of every job, so it has 8 records in W1 (0.818–1.156 τ0) and 4 in W2 (1.5–2.0 τ0);
  Gallope ran each case in one job and has 1 (1.0 τ0) and 2 (1.5, 2.0 τ0). The verdict ('ES better') and the
  counts (n 4, need 3, plus 3, minus 0) are the same; the mean k* of ES and B* are 27.51 and 21.73 against
  Calliope's 27.37 and 22.10.
- s13d ran Stages 13–16 (70 cases). Long turbulent runs leave the s13d trajectory at the rounding level
  and the difference grows, so the comparison of these runs is listed and not judged (plan of this port);
  only the 10-step decomposition cases are judged. For the 14 aborted cases both codes write the restart
  files at t = 0 only. Of the other 56, 52 agree within the tolerances of `compare_s13d.py` (restart 1e-11,
  t 1e-12); the 4 others are:
  - `stage15/a/turb1_N64_ESf` and `stage13d/nomass/e0.5_C2`: t and dt equal s13d in every row (70 and 95
    changes of dt); the restart fields differ by 6.4e-11 and 3.4e-10. Both runs end in a collapsing hole
    with dt 3.5e-11 and 2.8e-13, where Calliope cut them.
  - `stage15/b/turb1_N128_ES`: a decrease of dt (0.000895 → 0.000715) came one step later in Gallope
    (t = 0.85851 against 0.85761). From there the rows have different t, and the runs end at t differing
    by 1.1e-4 relative, after 20 and 17 changes of dt.
  - `stage15/b/turb1_N128_ESf`: a decrease at t = 2.04092 gave dt 0.000329 in Gallope and 0.000328 in s13d,
    near the stop where dt falls to 1e-11; end t differ by 2.1e-12, the restart fields by 1.6e-4.

## Figures

`plot_overlay.py` draws Calliope's figures of the paper with the Gallope results over them (`fig/`, not
committed). Euler1D 4/4, MHD1D 6/6, Smooth_wave 4/5 (`cp_diag`: Calliope's script cannot print γ = None of the
2D N16 `noeta`/`svvonly` runs), OT2D 4/4. The Gallope curves lie on Calliope's in every panel checked
(`plot_overlay.py` leaves out a Gallope curve identical to Calliope's: 9 of the 17 spectra of `ot_spec`).
The maps of `ot_rho` and `ot_J` are drawn for each code and as the difference (Gallope − Calliope)/max|Calliope|:
at most 3.8e-8 (ρ) and 4.5e-8 (J) at isolated points, and 0 to 5e-14 in four of the eight maps. The figure
data hold ρ, B_x and B_z of the plane in float32 (`figdata.py`), and J is computed from the stored B. A stored
value differs by one float32 unit (about 6e-8 relative) where the double values of the two codes fall on either
side of a float32 rounding boundary.

Turb3D 11/11. The peaks of ε(t), the spectra at 1τ0, the density PDFs (`turb_rho_pdf`) and the dt histories
of `turb_track` equal Calliope's to the digits printed; `turb_track` and `turb_rho_pdf` draw no Gallope curve
that differs from Calliope's. Differences:
- `turb_spec` has Gallope curves at 1τ0 only. With t_end = 3π, `turb1_N128_a` and `turb1_N256` stop at 2.99997
  and 2.99994 τ0 (the run ends at tt ≥ t_end − dt/2), and `analyze_stage4.py` takes a spectrum at 3τ0 only
  if 3τ0 ≤ t_last. Calliope's runs, made with older executables, stopped past 3τ0. `plot_overlay.py` passes
  `keys=('1.00',)` to the Gallope pass of this figure.
- `turb_slice`: the z = 0 plane of 256³ at 1τ0, max|Gallope − Calliope|/max|Calliope| 3.0e-4 and 2.2e-4 in the
  two panels; ζ_max 3.17549e-2 against 3.17520e-2, area of ζ > 0.01 ζ_max 0.92336 against 0.92334. The
  Calliope run is Stage 6's older executable, which has no s13d output.
- `turb_coef`: m, L and R of the five Stage 8 runs equal to the 4 digits printed. The control repeats
  `turb1_N128_t` with max|Δ ln ρ_min| = 0, against 2.3e-3 in Calliope (an interpolation error at the kinks of
  ρ_min(t), `REPORT_stage8`).
- `turb_nomass`: the run without mass correction is equal (14,397 steps, d ln dt/d ln ρ_min 1.4991); the control
  takes 6,591 steps against 6,589.
- `turb_ms`: the step counts differ from Calliope's older executables (M_s = 2: 45,943 against 45,402;
  M_s = 3: 46,396 against 46,255; M_s = 4 at 64³: 51,054 against 50,733; at 128³: 16,045 against 15,885); dt0,
  t_last, min ρ_min, the onset of dt < dt0/10 and the verdicts are equal.
- `turb_K` and `turb_kappa` use each code's own ES-256 as reference, so κ differs as listed under Turb3D
  (e.g. ES, E_B W1: 0.7055 against 0.7127); the verdicts and B* = H(8) are the same. `turb_kappa` has fewer
  Gallope points because of the record counts above.

## Cost

TSUBAME4.0 charges P = N · c · (0.7 · max(T, 300 s) + 0.1 · h_rt) / 3600 points for a job of N resource units
of type c (node_f 1, node_h 0.5, node_q 0.25, cpu_16 0.06), wall time T and requested time h_rt. From `qacct`
of every job:

| jobs | count | GPU-h (allocated) | points |
|---|---|---|---|
| gates G0–G3 | 5 | 0.37 | 0.26 |
| throughput test (`bench_os`) | 1 | 0.09 | 0.03 |
| Euler1D, MHD1D, Smooth_wave | 14 | 26.47 | 6.04 |
| OT2D, Turb3D | 16 | 11.35 | 3.97 |
| judgments (`cpu_16`) | 4 | – | 0.09 |
| total | 40 | 38.28 | 10.39 |

The plan estimated 30–60 points with a stop at 60. The group balance (`t4-user-info group point`) was 207.01
points before the first job and 184.73 at 14:50 on 2026-11-16, with 16.11 points in deposit and no job in the
queue. Both values were the same at 14:19, although the judgment job 9414843 ran from 14:22 to 14:37, so they
are not the final charge.

The 1D jobs cost most (6.04 points). Two of them (`euler1d_f`, `euler1d_h`) packed 4 cases per GPU without
MPS, which ran 8.4 times slower than one case per GPU; `euler1d_h` alone took 0.78 points for 7,012 s. The
longest job is `turb1_N256` alone on one GPU (16,055 s, 1.38 points). Solve time per step on one H100:
0.0263 s at 64³, 0.0553 s at 128³, 0.401 s at 256³.

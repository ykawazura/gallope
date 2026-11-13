# EMHD validation on TSUBAME4.0

**Verdict: PASS.** All checks V0–V9 pass, except V4, which does not apply to EMHD. The Gallope
EMHD model reproduces Calliope EMHD to 9.4e-16 after 200 steps of a deterministic OT3 run.
Its results do not depend on the GPU count, and a restarted run is bitwise identical to an
uninterrupted one. The invariants W and H* drift at third order in dt. Helical whistler
modes have the analytic frequency to 2e-10 or better, on 1 and 4 GPUs. div b stays at round-off.
The energy budget closes, and fixed-power forcing injects the prescribed power.

Date: 2026-11-13. Machine: TSUBAME4.0 (H100 SXM5, 4 per node). Gallope `develop`, built with
`nvhpc/26.1_cuda13.1`, `openmpi/5.0.2-nvhpc`, `MODEL = EMHD`.

## How to reproduce

Run from the repository root on TSUBAME4.0:

```bash
CASES=$(./tests/EMHD/setup_runs.sh)        RESOURCE=node_f H_RT=0:30:00 NAME=emhd ./tests/run_cases-tsubame.sh
CASES=$(./tests/EMHD/setup_runs.sh linear) RESOURCE=node_f H_RT=0:20:00 NAME=elin ./tests/run_cases-tsubame.sh
# Calliope reference (V8): calliope_dev develop 1d924a9 (git archive) + tests/EMHD/calliope_OT3.patch
mkdir -p runs/c_emhd_ot3 && cp tests/EMHD/calliope.in runs/c_emhd_ot3/
CASES="runs/c_emhd_ot3:16" CALLIOPE_ROOT=<patched Calliope tree> RESOURCE=cpu_16 ./tests/run_calliope-tsubame.sh
```

Then check the results; this needs python with numpy and netCDF4:

```bash
python3 tests/EMHD/check_runs.py runs            # V0, V2-V6, V9, div b, energy budget
python3 tests/EMHD/check_runs.py runs linear     # V0, V7
python3 tests/compare_calliope.py --model EMHD --gallope runs/e_ot3_np1 \
    --calliope runs/c_emhd_ot3 --nlx 64 --nly 64 --nlz 32       # V8 (also e_ot3_np4)
```

- Input: every case is `tests/EMHD/gallope.in` modified by the `sed` lines in
  `tests/EMHD/setup_runs.sh`.
- Job script: `tests/run_cases-tsubame.sh`. It runs the cases in sequence inside one job; each
  case is `mpirun -n <np> bin/gallope-EMHD` launched through `batch/gpu-bind-tsubame.sh`.
- The Calliope job script is `tests/run_calliope-tsubame.sh` (16 MPI ranks, gcc + openmpi).
- Jobs: 9374328 (main set), 9374329 (linear set), 9374334 (Calliope).

The base case has the following setup:
- grid 64 × 64 × 32, de = 0.1, guide field b0 = z;
- hyper-resistivity η (k²/k²_max)^2 with η = 100, which is zero at k = 0;
- deterministic OT3 initial condition, bx = −sin(y + z), by = sin(2x + z), bz = 0;
- dt = 1e-3 held fixed, 200 steps, eSSPIFRK3.

The whistler CFL bound is about 1e-2 here, so dt is never reset. The random-IC cases use
b1 = 0.5 in |k_d| ≤ 3, b0 = 0 and init_seed = 1. The forced cases stir bx, by and bz at
|kx|, |ky| ≤ 1 and kz = 1, with stir_seed = 1.

## Results

| check | cases | criterion | result |
|---|---|---|---|
| V0 input read | `e_ot3_np1`, `e_frc_np1`, `e_ideal_dt2`, `e_rnd_ideal_dt2`, all four `e_lin_*` | η, η exponent, de and b0(t = 0) recorded in `gallope.out.nc` equal the input | as intended |
| V1 build | `MODEL = EMHD` | no error; no implicit array copy | implicit copies only of reduction scalars (`rms`, `*_sum`, `kw`, `bm`, `vx`, `vy`, `vz`) |
| V2 smoke | all 10 main-set cases | finite, reaches `t_end` | PASS |
| V3 decomposition | `e_ot3_np1` vs `e_ot3_np4` | series rel. diff < 1e-12 | 2.0e-13 in `wmag_dot_sum`, ≤ 2.4e-16 in the others; restart files identical |
| V4 redundancy | — | N/A: EMHD runs with P_m = P_s = 1 only | `e_pm2` (P_m = 2) aborts before the first step |
| V5 restart | `e_rs1` → `e_rs2` (100 + 100 steps) vs `e_ot3_np1` | restart bitwise equal | identical, tt = 0.20000000000000015 in both |
| V6 ideal invariants | `e_ideal_dt{1,2}` (OT3, b0 = z); `e_rnd_ideal_dt{1,2}` (random, b0 = 0); dt = 2e-3, 1e-3 | drift ratio ≈ 8 | W (OT3): 1.640e-8 / 2.049e-9 = 8.00; W (random): 1.534e-7 / 1.918e-8 = 8.00; H* (random): 3.982e-9 / 4.914e-10 = 8.10 |
| V7 linear waves | `e_lin_001_de{0,01}` (1 GPU), `e_lin_3m42_de{01,1}` (4 GPUs) | ω and amplitude rel. err < 1e-6; leakage < 1e-10 of the amplitude | ω: 3.3e-14, 3.1e-14, 1.6e-10, 5.1e-15; amplitude: 4.2e-11, 4.0e-11, 1.0e-7, 3.5e-13; leakage ≤ 3.9e-16; b0 unchanged exactly |
| V8 Calliope | `e_ot3_np1`, `e_ot3_np4` vs `c_emhd_ot3` | combined rel. L2 < 1e-10 | 9.4e-16 for both |
| V9 forcing | `e_frc_np1`, `e_frc_np4` (bx, by, bz) | p_ext ≈ ene_inj = 1; np1 vs np4 < 1e-10 | [1.0012, 1.0116]; np1 vs np4 1.9e-14 |
| div b | all 10 main-set cases and the linear set | divb_rel < 1e-12 | ≤ 4.8e-16 |
| energy budget | `e_ot3_np1` (decaying) | see below, < 1e-3 | 4.9e-6 |
| energy budget | `e_frc_np1` (driven) | \|dW/dt − (inj − dissip)\| / max\|dW/dt\| < 1e-2 | 6.4e-5 |

The V3 difference in `wmag_dot_sum` is round-off in the difference quotient (W − W_old)/dt.
W ≈ 1 differs by 2.2e-16 between the two runs, and |dW/dt| is only 4.5e-4.

The V7 frequencies follow from the implemented operator. Let d = 1 + de²k². For a helical
mode with k × e_s = i s|k| e_s, the fluctuation u = s|k| b and b* = d b are parallel to b.
The only surviving term is then b0 × u, and the mode is an exact solution of the nonlinear
equations:

  b_k(t) = b_k(0) exp(i s (k·b0)|k| t/d) exp(−η (k²/k²_max)^2 t/d).

The expected ω are 1 and −1/1.01 at (kx, ky, kz) = (0, 0, 1) with s = +1, de = 0 and
s = −1, de = 0.1. At (3, −4, 2) they are 2√29/1.29 = 8.3491 (s = +1, de = 0.1) and
−2√29/30 = −0.3590 (s = −1, de = 1). The (3, −4, 2) mode sits on the last of the four ky
ranks. The 1.0e-7 amplitude error at (3, −4, 2), de = 0.1, is the eSSPIFRK3 error of the
whistler term. b0 enters through b* × u, so that term is explicit, and the error is
N(ω dt)⁴/24 = 500 × (8.35e-3)⁴/24 = 1.0e-7.

The budget of a decaying run needs the integrator drift removed. The dissipation of OT3 is
at most 4.5e-4, and the eSSPIFRK3 drift of W at dt = 1e-3 is about 1.3e-8 per unit time
(V6). The check subtracts dW/dt of the same run without dissipation, `e_ideal_dt2`, and
divides the residual by the maximum dissipation. Without that subtraction, the residual is
3.1e-5 of the maximum dissipation.

## Differences from Calliope EMHD

- **CFL.** Gallope adds the whistler bound |b*|_max Kw/C to the velocity bound, with
  Kw = max k²/(1 + de²k²) over the dealiased modes and C = √3 for eSSPIFRK3. Calliope uses |u| only.
- **Initial conditions.**
  - Gallope adds `OT3`, and `single_mode` is a helical whistler mode.
  - In `random`, Gallope projects onto div b = 0 and scales the rms of |b| to b1.
    Calliope scales each component by b1(d), which breaks div b = 0, and its kmin test
    removes every mode with a zero component.
  - Every initial condition except `restart` sets b0 as the k = 0 mode.
- **Dissipation and forcing.**
  - The hyper-resistivity is zero at k = 0, so the guide field is not dissipated when
    eta_exp = 0.
  - Gallope applies `fix_power`; Calliope reads it but has no `normalize_force`.
- **Diagnostics.**
  - Gallope adds `hel_sum` (H*) and `divb_rel`.
  - It writes b0 as the real part of the k = 0 mode; Calliope writes its absolute value.
- **Parallelization.** EMHD runs with P_m = P_s = 1 only.

## Issues found during validation

- **Calliope EMHD has no OT3 initial condition.** `tests/EMHD/calliope_OT3.patch` adds the
  same field and b0, in the deployed reference copy only. Calliope itself is unchanged.
- **P_m = 2 is rejected by the common grid check.** With nm = 1, `src/grid.F90` stops with
  "nm has to divide evenly by the comm_m rank count (P_m)" before the EMHD-specific guard in
  `init_model_specific` is reached. `check_runs.py` accepts either message.
- **Calliope message about `&forced_fields`.** The Calliope run prints "Reading
  force_parameters failed" once. It comes from proc0 reading `&forced_fields`, which the
  undriven `calliope.in` does not have. `&force` itself is read, and the forcing is off.
- All inputs have one namelist item per line, since nvfortran stops reading a group at `;`.

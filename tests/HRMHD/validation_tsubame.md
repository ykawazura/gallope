# HRMHD validation on TSUBAME4.0

**Verdict: PASS.** All checks V0–V9 pass. The Gallope HRMHD model reproduces Calliope HRMHD to
1.6e-15 after 200 steps of a deterministic OT3 run. Its results do not depend on the GPU count
or the redundant decomposition, and a restarted run is bitwise identical to an uninterrupted
one. The invariants W and X drift at third order in dt. Both linear branches, KAW and ICW,
have the analytic frequency to 1e-11 or better. The energy budget closes. Fixed-power
forcing, including the Elsasser variant, injects the prescribed power.

Date: 2026-11-13. Machine: TSUBAME4.0 (H100 SXM5, 4 per node). Gallope `develop`, built with
`nvhpc/26.1_cuda13.1`, `openmpi/5.0.2-nvhpc`, `MODEL = HRMHD`.

## How to reproduce

Run from the repository root on TSUBAME4.0:

```bash
CASES=$(./tests/HRMHD/setup_runs.sh)        RESOURCE=node_f H_RT=0:30:00 NAME=hrmhd ./tests/run_cases-tsubame.sh
CASES=$(./tests/HRMHD/setup_runs.sh linear) RESOURCE=node_f H_RT=0:30:00 NAME=hlin  ./tests/run_cases-tsubame.sh
# Calliope reference (V8): calliope_dev develop 1d924a9 (git archive) + tests/HRMHD/calliope_OT3.patch
mkdir -p runs/c_ot3_n200 && cp tests/HRMHD/calliope.in runs/c_ot3_n200/
CASES="runs/c_ot3_n200:16" RESOURCE=cpu_16 ./tests/run_calliope-tsubame.sh
```

Then check the results; this needs python with numpy and netCDF4:

```bash
python3 tests/HRMHD/check_runs.py runs            # V0, V2-V6, V9, energy budget
python3 tests/HRMHD/check_runs.py runs linear     # V0, V7, rho -> 0 limit
python3 tests/compare_calliope.py --model HRMHD --gallope runs/h_ot3_np1 \
    --calliope runs/c_ot3_n200 --nlx 64 --nly 64 --nlz 32          # V8 (also h_ot3_np4)
```

- Input: every case is `tests/HRMHD/gallope.in` modified by the `sed` lines in
  `tests/HRMHD/setup_runs.sh`.
- Job script: `tests/run_cases-tsubame.sh`. It runs the cases in sequence inside one job; each
  case is `mpirun -n <np> bin/gallope-HRMHD` launched through `batch/gpu-bind-tsubame.sh`.
- The Calliope job script is `tests/run_calliope-tsubame.sh` (16 MPI ranks, gcc + openmpi).
- Jobs: 9373981 (main set), 9373987 (linear and Elsasser), 9374019 (`h_rho0_ideal`),
  9373898 (Calliope).

The base case has the following setup:
- grid 64 × 64 × 32, ρ = 0.5, σ = 1;
- hyper-dissipation of order 4 on all four fields, coefficient 10 perpendicular and
  1 parallel;
- deterministic OT3 initial condition;
- dt = 1e-3 held fixed, 200 steps, eSSPIFRK3.

## Results

| check | cases | criterion | result |
|---|---|---|---|
| V0 input read | `h_ot3_np1`, `h_frc_els_np1`, `h_ideal_dt2`, `h_lin_kaw_401`, `r_ot3` | parameters recorded in `gallope.out.nc` equal the input | as intended |
| V1 build | `MODEL = HRMHD` | no error; no implicit array copy | implicit copies only of reduction scalars (`*_sum`, `max_vel_*`, `chk`) |
| V2 smoke | all 14 main-set cases | finite, reaches `t_end` | PASS |
| V3 decomposition | `h_ot3_np1` vs `h_ot3_np4` | series rel. diff < 1e-10 | 4.7e-16; restart files identical |
| V4 redundancy | `h_ot3_np2` (P_fft = 2) vs `h_ot3_pm2` (P_fft = 2, P_m = 2) | restart bitwise equal | identical; `check_redundant_consistency` 22 × PASS |
| V5 restart | `h_rs1` → `h_rs2` (100 + 100 steps) vs `h_ot3_np1` | restart bitwise equal | identical, tt = 0.20000000000000015 in both |
| V6 ideal invariants | `h_ideal_dt1` (dt = 2e-3), `h_ideal_dt2` (dt = 1e-3) | drift ratio ≈ 8 | W: 5.63e-9 / 7.03e-10 = 8.00; X: 6.15e-9 / 7.68e-10 = 8.00 |
| V7 linear waves | `h_lin_{kaw,icw}_401` (1 GPU), `h_lin_{kaw,icw}_3m42` (4 GPUs) | ω rel. err < 1e-3 | KAW 1.1e-12, 3.5e-11; ICW 1.0e-14, 7.9e-15; no leakage to other modes |
| ρ → 0 limit | `h_rho0` vs RMHD `r_ot3` | φ, ψ agree | max rel. diff 0 (bitwise); upa = bpa = 0 |
| V8 Calliope | `h_ot3_np1`, `h_ot3_np4` vs `c_ot3_n200` | combined rel. L2 < 1e-10 | 1.6e-15 for both |
| V9 forcing | `h_frc_np1`, `h_frc_np4` (φ, ψ, upa, bpa) | p_aw + p_compr ≈ ene_inj = 1 | [1.0011, 1.0214]; np1 vs np4 3.5e-16 |
| V9 Elsasser forcing | `h_frc_els_np1`, `h_frc_els_np4`, xhl_inj = 0.5 | p_aw ≈ 1, p_xhl ≈ 0.5, p_compr = 0 | p_aw [1.0011, 1.0124], p_xhl [0.5004, 0.5116], p_compr = 0; np1 vs np4 4.4e-16 |
| V9 zeroing (device kernel) | `h_frc_els_np1`, `h_frc_els_np4` | unforced upa, bpa = 0 on the 18 stirred modes | exactly 0 on 1 and 4 GPUs; the neighbouring kz = 3 modes are O(1e-3) |
| energy budget | `h_ot3_np1`, `h_rho0` (decaying) | see below, < 1e-3 | 1.7e-5, 3.6e-5 |
| energy budget | `h_frc_np1`, `h_frc_els_np1` (driven) | \|dW/dt − (inj − dissip)\| / max\|dW/dt\| < 1e-2 | 1.5e-4, 7.2e-5 |

The V7 eigenvectors follow from the implemented linear operator. With λ = ω/k∥, each mode
satisfies

  λ⁴ − (1 + σ² + ρ²k⊥²) λ² + σ² = 0,
  (φ, ψ, upa, bpa) = (1, λ, σ(λ² − 1)/(ρλ), (λ² − 1)/ρ).

The expected ω are 2.4142, 0.4142 at (kx, ky, kz) = (4, 0, 1) and 5.7016, 0.7016 at (3, −4, 2).
The (3, −4, 2) mode sits on the last of the four ky ranks.

The budget of a decaying run needs the integrator drift removed. Its hyper-dissipation
(1e-9 to 1e-8) is as small as the eSSPIFRK3 truncation drift at dt = 1e-3, which is
about 3.5e-9 per unit time (V6). The check therefore subtracts dW/dt of the same run without
dissipation, `h_ideal_dt2` or `h_rho0_ideal`, and divides the residual by the dissipation.
Without that subtraction, the residual is 19 % of max\|dW/dt\|, and it is the drift itself.

## Issues found during validation

- **nvfortran silently truncates a namelist group at `;`.** An earlier version of the input
  had lines such as `nupe_x = 1d+1 ; nupe_x_exp = 4`. nvfortran read `nupe_x` and left every
  later item of `&physical_parameters` at its default (exponents 1, coefficients 0) without an
  error. Calliope (gfortran, one item per line) read everything. V8 exposed it: the
  difference was 7.9e-4 and grew linearly in time. All inputs now have one item per line,
  and V0 checks the parameters recorded in the output.
- **Calliope HRMHD OT3 initial condition.** `src/model/HRMHD/fields.F90` transforms `upa_r`
  and `bpa_r` into `psi`, which overwrites ψ. `tests/HRMHD/calliope_OT3.patch` fixes this in
  the deployed reference copy only. Calliope itself is unchanged.

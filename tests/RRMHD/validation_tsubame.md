# RRMHD validation on TSUBAME4.0

**Verdict: PASS.** All checks V0–V8 pass. After 200 steps of a deterministic OT3 run with
shear, rotation and both background gradients, the Gallope RRMHD model reproduces Calliope
RRMHD to 2.1e-15. Its results do not depend on the GPU count or the redundant decomposition,
and a restarted run is bitwise identical to an uninterrupted one. The energy drifts at third
order in dt. All five linear cases follow expm(iRt) of the implemented operator, which agrees
with the dispersion relation derived from the Calliope equations. The energy budget closes,
and the Elsasser energies satisfy their identities to round-off.

Date: 2026-11-13. Machine: TSUBAME4.0 (H100 SXM5, 4 per node). Gallope `develop`, built with
`nvhpc/26.1_cuda13.1`, `openmpi/5.0.2-nvhpc`, `MODEL = RRMHD`.

## How to reproduce

Run from the repository root on TSUBAME4.0:

```bash
CASES=$(./tests/RRMHD/setup_runs.sh)        RESOURCE=node_f H_RT=0:30:00 NAME=rrmhd ./tests/run_cases-tsubame.sh
CASES=$(./tests/RRMHD/setup_runs.sh linear) RESOURCE=node_f H_RT=0:30:00 NAME=rrlin ./tests/run_cases-tsubame.sh
# Calliope reference (V8): calliope_dev develop 1d924a9 (git archive) + tests/RRMHD/calliope_OT3.patch
mkdir -p runs/rc_ot3_n200 && cp tests/RRMHD/calliope.in runs/rc_ot3_n200/
CASES="runs/rc_ot3_n200:16" CALLIOPE_ROOT=$HOME/work/calliope_develop_rrmhd RESOURCE=cpu_16 \
    ./tests/run_calliope-tsubame.sh
```

Then check the results; this needs python with numpy, scipy and netCDF4:

```bash
python3 tests/RRMHD/check_runs.py runs            # V0, V2-V6, budget, Elsasser, driven
python3 tests/RRMHD/check_runs.py runs linear     # V0, V7
python3 tests/compare_calliope.py --model RRMHD --gallope runs/rr_ot3_np1 \
    --calliope runs/rc_ot3_n200 --nlx 64 --nly 64 --nlz 32          # V8 (also rr_ot3_np4)
```

- Input: every case is `tests/RRMHD/gallope.in` modified by the `sed` lines in
  `tests/RRMHD/setup_runs.sh`. The Calliope input is `tests/RRMHD/calliope.in`.
- Job script: `tests/run_cases-tsubame.sh`. It runs the cases in sequence inside one job; each
  case is `mpirun -n <np> bin/gallope-RRMHD` launched through `batch/gpu-bind-tsubame.sh`.
- The Calliope job script is `tests/run_calliope-tsubame.sh` (16 MPI ranks, gcc + openmpi).
- Jobs: 9374150 (main set, linear set and the `input_example/RRMHD.in` smoke run, in one job),
  9374070 (Calliope).

The base case has the following setup:
- grid 64 × 64 × 32, β = 1, γ = 5/3, so A = 1 + 2/(βγ) = 2.2 and c = βγ/2 = 5/6;
- q = 1.5, κ_b = 0.2, κ_p = 0.3;
- hyper-dissipation of order 4 on all four fields, coefficient 10 perpendicular and
  1 parallel;
- deterministic OT3 initial condition;
- dt = 1e-3 held fixed (the CFL bound is 0.196), 200 steps, eSSPIFRK3.

## Results

| check | cases | criterion | result |
|---|---|---|---|
| V0 input read | all 14 cases, `tests/RRMHD/gallope.in` | parameters recorded in `gallope.out.nc` equal the input; the base input holds the intended values | as intended |
| V1 build | `MODEL = RRMHD` | no error; no implicit array copy | implicit copies only of reduction scalars (`*_sum`, `max_vel_*`, `chk`) |
| V2 smoke | all 9 main-set cases | finite, reaches `t_end` | PASS |
| V3 decomposition | `rr_ot3_np1` vs `rr_ot3_np4` | series and restart rel. diff < 1e-12 | series 4.3e-16; restart files identical |
| V4 redundancy | `rr_ot3_np2` (P_fft = 2) vs `rr_ot3_pm2` (P_fft = 2, P_m = 2) | restart bitwise equal | identical; `check_redundant_consistency` 22 × PASS |
| V5 restart | `rr_rs1` → `rr_rs2` (100 + 100 steps) vs `rr_ot3_np1` | restart bitwise equal | identical, tt = 0.20000000000000015 in both |
| V6 ideal energy | `rr_ideal_dt1` (dt = 2e-3), `rr_ideal_dt2` (dt = 1e-3), q = κ = 0 | drift ratio ≈ 8 | 3.827e-9 / 4.784e-10 = 8.00; p_aw_rot + p_compr_rot = 0 and p_grd = 0 exactly |
| V7 Alfvén and slow waves | `rr_lin_ap`, `rr_lin_am`, (4, 0, 1), q = κ = 0 | \|X − expm(iRt) X(0)\| / \|X\| < 1e-6; ω rel. err < 1e-6 | 2.8e-11; ω = ±(1, 1, 0.674199862, 0.674199862) per component, as predicted; E drift 7.5e-11 |
| V7 rotation and gradients | `rr_lin_rot_a` (4 GPUs), `rr_lin_rot_b` (1 GPU), eigenvectors at (3, −4, 2) | same | 1.9e-9, 4.3e-11; ω = −2.614574714, 1.009233792 in all four components; E drift 3.9e-9, 8.7e-11 |
| V7 unstable mode | `rr_lin_unst`, (3, −4, 1) | expm < 1e-6; \|dE/dt − Σp_*\| / max(\|p_*\|, \|dE/dt\|) < 1e-6 | 1.8e-10, growth rate 0.414429; budget 6.7e-8 (max\|p_*\| = 17.1) |
| V7, all five cases | | no leakage to other modes; omg = −Kφ; quartic roots = −eig(R) | leakage 0; omg error 0; roots differ by ≤ 1.8e-15 |
| V8 Calliope | `rr_ot3_np1`, `rr_ot3_np4` vs `rc_ot3_n200` | combined rel. L2 < 1e-10 | 2.1e-15 for both (φ 2.4e-15, ψ 1.5e-15, u∥ 2.0e-15, b∥ 1.6e-15); restart times equal; axis mapping 3.8e-17 / 3.7e-18 |
| energy budget | `rr_ot3_np1` (eSSPIFRK3), `rr_gear3` | \|dE/dt − (Σp_* − dissip)\| / max\|dE/dt\| < 1e-2 | 8.6e-7, 8.9e-5 |
| Elsasser identities | `rr_ot3_np1`, `rr_ot3_np4` | see below; rel. < 1e-12, series and every (kz, k⊥) bin | ≤ 4.4e-16 |
| `driven = T` | `rr_driven` | stops with a message, no final restart | as intended |
| example input | `ex_rrmhd`: `input_example/RRMHD.in` with nstep = 100 | finite | reaches t = 9.25; the CFL lowers dt from 0.1 to 0.05; E grows from 0.416 to 0.991 |

The Elsasser energies follow Calliope: they are taken from the current fields without a
factor 1/2,

  zpep2 = K|φ + ψ|², zpem2 = K|φ − ψ|², zpap2 = |u∥ + √A b∥|², zpam2 = |u∥ − √A b∥|²,

so zpep2 + zpem2 = 4(upe2 + bpe2) and zpap2 + zpam2 = 4(upa2 + bpa2), where
bpa2 = A|b∥|²/2. These identities hold to round-off, which the check verifies. Applied to
output written with an earlier convention (midpoint fields, factor 1/2), the same check
reports 0.5 and fails.

## Linear operator and eigen-relations (V7)

Gallope advances X = (φ, ψ, u∥, b∥) at one Fourier mode as dX/dt = iRX, with K = kx² + ky²
and the real matrix

```
R = [[ 0,                        kz,       2ky/K,  -2(κ_b + cκ_p)ky/K ],
     [ kz,                       0,        0,       0                 ],
     [ (2 - q)ky,               -cκ_p ky,  0,       kz                ],
     [ -(2κ_b + (c - 1)κ_p)ky/A, q ky/A,   kz/A,    0                 ]]
```

`tests/RRMHD/check_runs.py` builds R from `advance::get_ext_terms` and compares the final
state with expm(iRt) X(0). Separately, it derives the dispersion relation from the Calliope
equations. For perturbations ∝ exp(−iωt) it is

  ω⁴ + a₂ω² + a₁ω + a₀ = 0,
  a₂ = −[AKkz² + Kkz² + (4 − 2q)Aky² + 2c²κ_p²ky² + 6cκ_bκ_p ky² − 2cκ_p²ky² + 4κ_b²ky²
        − 2κ_bκ_p ky²]/(AK),
  a₁ = −8ky²kz(κ_b + cκ_p)/(AK),
  a₀ = kz²[Kkz² − 2q ky² − 2c²κ_p²ky² − 2cκ_bκ_p ky²]/(AK).

Its roots equal −eig(R) to 2e-15 in every linear case. Below, frequencies are the
eigenvalues w of R (X ∝ exp(iwt)), as measured per component from the restart files;
ω = −w.

- ky = 0, (4, 0, 1), q = κ = 0: the Alfvén waves have w = ±kz with ψ = ±φ, and the slow
  waves have w = ±kz/√A = ±0.674199862 with u∥ = ±√A b∥. `rr_lin_ap` starts from
  (1, 1, √A, 1) and `rr_lin_am` from (1, −1, −√A, 1); each is the sum of one Alfvén and one
  slow eigenvector.
- (3, −4, 2) with the base q and κ: the eigenvalues of R are −2.614574714, −0.603836354,
  1.009233792 and 2.209177276, all real. `rr_lin_rot_a` (4 GPUs; the mode sits on the last
  ky rank) starts from the eigenvector of −2.614574714,
  (1, −0.76494276, 2.52585664, −1.91954908), and `rr_lin_rot_b` from that of 1.009233792,
  (1, 1.98170138, 11.9493176, 6.03897688).
- (3, −4, 1) with the base q and κ: Kkz² = 25 is smaller than 2q ky² = 48, and the
  eigenvalues are −1.720035225, 1.290363041 and 0.214836092 ± 0.414429373i. `rr_lin_unst`
  starts from (1, 0.5, −0.3, 0.2), which contains the growing mode, and the run grows at
  the rate 0.414429.
- A real eigenvector rotates with one phase in all four components, so every source term
  p_* ∝ Im(x y*) vanishes and E is constant. eSSPIFRK3 then changes E only through its
  amplitude error, (w dt)⁴/12 per step. For 1000 steps at dt = 1e-3 this predicts 3.9e-9
  for `rr_lin_rot_a` and 8.6e-11 for `rr_lin_rot_b`. In `rr_lin_ap` and `rr_lin_am` the
  Alfvén part (8.3e-11) and the slow part (1.7e-11) weighted by their energies give 7.5e-11.
  All three agree with the measured drifts in the table.

The budget of the OT3 run needs no drift subtraction, unlike HRMHD. The linear
sources p_rot and p_grd dominate dE/dt (max|p_rot + p_grd| / max|dE/dt| = 1.00), and they
are much larger than the eSSPIFRK3 truncation drift.

## Issues found during validation

- **nvfortran silently truncates a namelist group at `;`.** The first version of
  `tests/RRMHD/gallope.in` had lines with two items separated by `;`, as found in HRMHD.
  All inputs now have one item per line. V0 checks that every parameter recorded in
  `gallope.out.nc` equals the input, and that the base input holds the intended values.
- **Calliope RRMHD OT3 initial condition.** `src/model/RRMHD/fields.F90` transforms `upa_r`
  and `bpa_r` into `psi`, which overwrites ψ; u∥ and b∥ never receive their transforms.
  `tests/RRMHD/calliope_OT3.patch` fixes this in the deployed reference copy only. Calliope
  itself is unchanged.
- **Calliope RRMHD `gear3` drops the gradient terms.** Its `gear3` branch computes the
  explicit terms inline without the κ_b and κ_p terms, which `get_ext_terms` (eSSPIFRK3)
  includes. Gallope uses the same explicit terms for both schemes; `rr_gear3` checks its
  budget.
- **Calliope RRMHD diagnostic dissipation.** It uses the perpendicular exponents for u∥ and
  b∥. The evolution uses the correct exponents, so the comparison is unaffected when all
  exponents are equal, as in `tests/RRMHD/calliope.in`.
- **Eigenmode budget test.** A budget check of the form |dE/dt − Σp_*| / max|dE/dt| is
  meaningless for a real eigenvector, since both sides vanish. For such cases the check
  uses the energy drift instead.
- **`cfl = 1` reaches the RK3 stability limit at the dealiasing cutoff.** The linear bound
  for the waves along z is dt ≤ cfl·dlz, the same form as dlz/vmax_KAW in HRMHD. After the
  2/3 rule the largest kz is about nlz/3, so kz·dt reaches 2π/3 ≈ 2.09 (1.96 for nlz = 32)
  on the Alfvén branch. This exceeds √3, the imaginary-axis limit of RK3. `cfl ≤ 0.82`
  keeps it below; `input_example/RRMHD.in` uses 0.8. The code is unchanged.
- **No time-step control when `nonlinear = F`.** As in Calliope, the CFL check runs only
  with the nonlinear terms, so linear runs keep the input dt. The linear set uses
  dt = 1e-3.
- `Reading force_parameters failed` is printed once by `src/force_common.F90`, because
  RRMHD inputs have no `&forced_fields` group. It is harmless.

## Intentional differences from Calliope RRMHD

- `gear3` keeps the κ_b and κ_p terms (see above).
- The diagnostic dissipation uses the exponent of each field, and every exponent defaults
  to 1 (Calliope has no default).
- The sources p_aw_rot, p_compr_rot, p_aw_grd and p_compr_grd are written as explicit real
  forms at the midpoint of the last step, so that dE/dt = Σp_* − dissipation holds to
  O(dt²).
- No `ntrans_*` bins; `dissip_aw_bin` and `dissip_compr_bin` are added.
- `driven = T` stops the run with `MPI_Abort` (Calliope ignores it).
- The CFL condition uses the Elsasser velocities max|u_⊥ ± b_⊥| per direction, and the
  time step also obeys the linear bound dt ≤ cfl·min(dlz, 1/ω_rot), with
  ω_rot = 2 + |q| + 2|κ_b + cκ_p| + c|κ_p| + |2κ_b + (c − 1)κ_p|. Calliope bounds only
  the E×B velocity. The linear bound is 0.196 in the base case.
- The OT3 initial condition keeps ψ (Calliope overwrites it; see above).
- The random initial condition is reproducible for `init_seed ≥ 0`.

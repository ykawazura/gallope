# Release checks for Gallope v1.1 on TSUBAME4.0

**Verdict for rc `9c01367`: ready to release.** Every check passes (R0-R5). The two
bugs found on `f529a53` (Issues 1 and 2) are fixed in `9c01367`, and the checks they
failed now pass: R4 (d) in all 20 cases and R3′ bit for bit.

Date: 2026-11-13 (`f529a53`) and 2026-11-14 (`9c01367`). Machine: TSUBAME4.0 (H100 SXM5, 4 per node). Release candidate
`9c01367` (develop), deployed with `git archive` and built with `nvhpc/26.1_cuda13.1`,
`openmpi/5.0.2-nvhpc`, `arch = tsubame`. The first round ran on `f529a53` and found
Issues 1 and 2; `9c01367` fixes both. R0, R3, R3′ and R4 were rerun on `9c01367`; R1, R2
and R5 are carried over from `f529a53` for the reason given in the next section.

## Changes from f529a53 to 9c01367 and the carried-over checks

```
$ git log --oneline f529a53..9c01367
9c01367 Fix [KRMHD]: OT2 initial vorticity
a93e71c Fix [ALL]: out3d records written back to back
$ git diff --stat f529a53 9c01367
 src/model/KRMHD/fields.F90 | 5 +++--
 src/mpiio_write_var.F90    | 5 ++---
 2 files changed, 5 insertions(+), 5 deletions(-)
```

- a93e71c removes the second `disp` advance in `src/mpiio_write_var.F90` (Issue 1).
  This file is included only in `mpiio_write_var`, which only `loop_io_3D` calls.
  `loop_io_3D` runs only when `write_intvl_3D > 0`: the calls at t = 0 and at the end
  test it (`src/gallope.F90:61, 110`), and the main-loop call (`:83`) never fires,
  because an interval of 0 is replaced by `-huge(1)` when the input is read
  (`src/diagnostics_common.F90:84`), so `floor(tt/intvl)` is −1 at every step.
  Restart files are written by `mpiio_write_one`, which does not change.
- 9c01367 sets `omg = -kprp2*phi` in KRMHD `init_OT2` (Issue 2). Only KRMHD runs with
  `init_type = 'OT2'` take this path.
- R1 and R2 (MHD_INCOMP and RMHD, `write_intvl_3D = 0`, restart comparisons) execute
  neither changed routine, so their `f529a53` results hold for `9c01367`.
- R5 has `write_intvl_3D = 0` in every case. Its KRMHD cases start from
  `init_type = 'random'` (`init_random` is unchanged), not from OT2. So the whole of R5,
  KRMHD included, holds for `9c01367`.
- R3 and R3′ (KRMHD, OT2 in most cases) and R4 (3D output) take the changed paths and
  were rerun.

## How to reproduce

Deploy three trees next to each other and build them (`make arch=tsubame MODEL=<model>`):

| tree | source | models | installed in the rc `bin/` as |
|---|---|---|---|
| rc | `git archive 9c01367` (R0, R3, R3′, R4); `git archive f529a53` (R1, R2, R5) | MHD_INCOMP, RMHD, KRMHD, HRMHD, RRMHD, EMHD | `gallope-<MODEL>` |
| v1.0 | `git archive v1.0` + `tests/release/v1.0_nvhpc26.patch` | MHD_INCOMP, RMHD | `gallope-<MODEL>_v1.0` |
| v1.0 + 3D | the v1.0 tree + `tests/release/enable_io3d.patch` | MHD_INCOMP, RMHD | `gallope-<MODEL>_v1.0_io3d` |
| f05d2c1 | `git archive f05d2c1` | KRMHD | `gallope-KRMHD_f05d2c1` |

Then, from the root of the rc tree:

```bash
for s in r12 r3 r4; do
  CASES=$(./tests/release/setup_runs.sh $s) RESOURCE=node_f H_RT=0:50:00 NAME=rel$s \
    ./tests/run_cases-tsubame.sh
done
CASES=$(./tests/release/setup_runs.sh r5) RESOURCE=node_q H_RT=0:50:00 NAME=relr5 \
  ./tests/run_cases-tsubame.sh                      # 1 GPU only
python3 tests/release/check_runs.py runs/rel        # needs numpy and netCDF4
```

- Input: every case is one of `tests/release/{MHD_INCOMP,RMHD,KRMHD}.in` or
  `tests/{HRMHD,RRMHD,EMHD}/gallope.in`, modified by the `sed` lines in
  `tests/release/setup_runs.sh`. All have a 64 × 64 × 32 grid, dt = 1e-3 held fixed and
  200 steps (t = 0.2), except where a row says otherwise.
- Job script: `tests/run_cases-tsubame.sh`. It runs the given cases in sequence inside
  one job; each case is `mpirun -n <np> bin/gallope-<name>` through
  `batch/gpu-bind-tsubame.sh`. One `node_f` job per set works; the `9c01367` reruns
  split the r3 and r4 cases by GPU count instead (1 GPU on `node_q`, 2 on `node_h`, 4 on
  `node_f`), keeping `r3_kr_rs1` before `r3_kr_rs2` in the same job.
- Jobs on `9c01367`: r3 and r4 split by GPU count, 9380302 (1 GPU, `node_q`), 9380303
  (2 GPUs, `node_h`), 9380304 (4 GPUs, `node_f`). Jobs on `f529a53`: 9376480 (r12),
  9376667 (r5); 9376481 (r3) and 9376482 (r4) found Issues 1 and 2.

## Results

`check_runs.py` prints one line per check: PASS, FAIL or INFO (numbers without a
verdict). Relative differences are max|a − b| / max|a| over a field or a series.
"Bitwise" means the restart files (`restart/*.dat`) are byte-identical and the final
`tt` is equal.

### R0: build and inputs

| item | result |
|---|---|
| rc build `9c01367`, 6 models | all OK, 0 severe (MHD_INCOMP 64 s, RMHD 28 s, KRMHD 31 s, HRMHD 29 s, RRMHD 26 s, EMHD 23 s) |
| rc build `f529a53`, 6 models | all OK (MHD_INCOMP 53 s, RMHD 27 s, KRMHD 31 s, HRMHD 29 s, RRMHD 27 s, EMHD 23 s) |
| v1.0 build | fails with nvhpc 26.1 without `v1.0_nvhpc26.patch` (see Issues); OK with it |
| f05d2c1 KRMHD build | OK (42 s) |
| runs: r3 + r4 on `9c01367` (38 cases), r12 + r5 on `f529a53` (23 + 27 cases) | all exit with rc = 0 and print `Finished the main loop`; every time series is finite |
| V0: input echo | 88/88: every numeric input item stored in `gallope.out.nc` equals the input value (except dt in the R5 A cases, replaced by `dt_adjust` by design); no namelist read failure except `Reading force_parameters failed` in the undriven RRMHD cases (no force group in the input); dt never decreased |

### R1: undriven, v1.0 vs v1.1

| case | v1.0 vs v1.1, 1 GPU | v1.0 vs v1.1, 4 GPUs | v1.1 1 vs 4 GPUs |
|---|---|---|---|
| MHD_INCOMP OT3 (`r1_mi`) | bitwise | bitwise | bitwise; series max rel 2.8e-15 |
| MHD_INCOMP OT3, shear q = 1.5, t = 1 (`r1_ms`) | bitwise | bitwise | bitwise; series max rel 1.0e-15 |
| RMHD OT3 (`r1_rm`) | bitwise | bitwise | bitwise; series max rel 7.5e-16 |

The v1.0 vs v1.1 series agree exactly (rel 0). The 1- vs 4-GPU series differ only in
the last bits of the volume sums.

### R2: driven (fixed `stir_seed`), v1.0 vs v1.1

| case | v1.1 1 vs 4 GPUs | v1.0 1 vs 4 GPUs (info) | v1.0 vs v1.1, 1 GPU |
|---|---|---|---|
| MHD_INCOMP, kz ∈ [1, 2] (`r2_mi`) | restart max rel 6.7e-15, series 1.0e-14 | restart 0.24, series 0.38 | bitwise; series rel 0 |
| RMHD, kz ∈ [1, 2] (`r2_rm`) | restart 1.3e-16, series 1.7e-15 | restart 0.25, series 1.5 | bitwise; series rel 0 |
| RMHD, kz = 0 in the band (`r2_rk0`) | restart 2.4e-16, series 3.6e-14 | — | differ: restart 0.21, series 0.30, from t = 0.009 (info) |

- A driven run of v1.1 on 1 GPU reproduces v1.0 bit for bit when the forcing band
  excludes kz = 0 (both models).
- With kz = 0 in the band, v1.0 and v1.1 differ from the first output (t = 0.009). This
  is the intended v1.1 change in the kz = 0 forcing, not a regression.
- v1.0 on 4 GPUs differs from v1.0 on 1 GPU at O(1) (the forcing bug fixed in v1.1);
  v1.1 agrees to < 1e-12. A driven run is not bitwise across GPU counts, because the
  forcing normalization sums in a different order.

### R3: KRMHD (9c01367, jobs 9380302-9380304)

| check | result |
|---|---|
| driven OT2 + g1, 1 vs 4 GPUs | restart max rel 6.0e-16, series 6.9e-13 (`p_xhl_sum`; < 1e-12) |
| P_m = 1 vs 2 at P_fft = 2 | bitwise |
| P_m = 1 vs 4 at P_fft = 1 | bitwise |
| undriven random, P_m = 1 vs 2 at P_fft = 2 | bitwise |
| `check_redundant_consistency` (pm2, pm4, random pm2) | 22 PASS, 0 FAIL each |
| restart: 100 + 100 steps vs 200 steps (random) | bitwise, final tt equal |
| R3′: KRMHD with g = 0 vs RMHD, same input, undriven OT2 | phi, omg, psi byte-identical; the 12 common series agree exactly (rel 0); max\|g\| = 0 |
| R3′: the same, Elsasser-driven OT2 (phi, psi forced; g not) | phi, omg, psi byte-identical; the 12 common series agree exactly (rel 0); max\|g\| = 0 |

- R3′ is now bitwise, not only < 1e-12. Undriven OT2: upe2 = 0.5, 0.49999797 (t = 0.009),
  0.49885944 (t = 0.2) in both models; upe2 + bpe2 = 0.9999999868 at t = 0.2.
- On `f529a53`, R3′ failed (Issue 2): phi, omg, psi max rel 16.8 (undriven) and 4.6
  (Elsasser-driven). The KRMHD flow energy dropped from 0.5 to 1.8e-5 at t = 0.009.
- The driven 1 vs 4 GPU case also starts from OT2, so its numbers changed with the fix
  (series max rel 7.3e-13 on `f529a53`, 6.9e-13 now). The rows with `init_type =
  'random'` (P_m random, restart) do not use `init_OT2` and give the same verdicts as on
  `f529a53`.

### R4: 3D output, write_intvl_3D = 0.08 (9c01367, jobs 9380302-9380304)

| check | result |
|---|---|
| (a) out3d 1 vs 4 GPUs, 6 models | byte-identical for MHD_INCOMP (7 files), RMHD, HRMHD, RRMHD, EMHD; KRMHD (driven) max rel 2.2e-16, as its restart (6.0e-16) |
| (b) out3d P_m = 1 vs 2 at P_fft = 2 | byte-identical for RMHD, KRMHD, HRMHD, RRMHD |
| (c) last record == `restart/<field>.dat` | all 20 rc cases, every field |
| (d) file size | all 20 rc cases: every `out3d/<field>.dat` is 4 456 448 B = len(time.dat) (4) × 16·(nlz/2+1)·nly·nlx (1 114 112 B) |
| (d) all-zero records | only where the field is zero: uz, bz at t = 0 (MHD_INCOMP OT3) and upa, bpa at t = 0 (HRMHD, RRMHD OT3); none in RMHD, KRMHD, EMHD |
| (d) `time.dat` times | 0, 0.079, 0.159, 0.2 in every case, equal to the main-loop rule to 1e-14 |
| v1.0 + 3D calls vs v1.1, 1 GPU | v1.0 record 2i == v1.1 record i byte for byte (i = 0..3), v1.0 records 1, 3, 5 all zero; time.dat and restart identical (MHD_INCOMP, RMHD) |
| v1.0 + 3D calls, 4 GPUs (info) | wrong: the file is not a whole number of records, and the last record ≠ restart in all 4 ky bands |
| cost | timer_io_3D 0.06 s = 0.7 % of the run (MHD_INCOMP, 4 GPUs, out3d 25.5 MiB); 0.1 % (KRMHD, 4 GPUs, 8.5 MiB); 0.3 % (EMHD, 1 GPU, 12.8 MiB) |

- A 3D record is written at the last step before each multiple of write_intvl_3D
  (t = 0.079, 0.159), as for the other outputs, plus t = 0 and the end.
- How (d) decides that a field is zero: for each group whose energy is in the time
  series (`u2_sum` = [ux, uy, uz], `b2_sum` = [bx, by, bz], `upa2_sum` = [upa],
  `bpa2_sum` = [bpa]), it computes E_3D = Σ_k w_kz |f_k|² over the group's records
  (w = 1 at kz = 0 and the Nyquist plane, 2 elsewhere). E_3D / series must be one
  constant at all 4 record times (spread (max − min)/max < 1e-10), and E_3D must be
  exactly 0 where the series is 0. In all 10 cases with such groups (MHD_INCOMP, HRMHD,
  RRMHD), the ratio is 2 (0.909 for RRMHD `bpa2_sum`) with spread ≤ 5.6e-16, and the
  series is 0 at t = 0 for upa and bpa. An all-zero record of a field outside these
  groups fails (d).
- The check catches what it should. On the `f529a53` output, (d) fails in all 20 cases
  (file size). With one record zeroed by hand in `9c01367` output, (d) fails in each
  of three tests: ux at t = 0.159 (MHD_INCOMP; `u2_sum` spread 0.49), upa at t = 0.079
  (HRMHD; spread 1.0) and phi at t = 0.079 (RMHD; no energy group, so unexplained).
- On `f529a53`, (d) failed in all 20 cases (Issue 1): 7 records for 4 `time.dat`
  entries, records 1, 3, 5 all zero.

### R5: dt_adjust (f529a53, job 9376667)

Case A starts with dt = 2e-3 and a `dt_adjust` file holding 1.0d-3; case B runs with
dt = 1e-3 from the start. Both are undriven, 1 GPU, 200 steps; KRMHD starts from
`init_type = 'random'`. A must reproduce B bit for bit.

| model | eSSPIFRK3: A vs B | gear3: A vs B |
|---|---|---|
| MHD_INCOMP | bitwise; series rel 0 | bitwise; series rel 0 |
| RMHD | bitwise; series rel 0 | bitwise; series rel 0 |
| KRMHD | bitwise; series rel 0 | bitwise; series rel 0 |
| HRMHD | bitwise; series rel 0 | bitwise; series rel 0 |
| RRMHD | bitwise; series rel 0 | bitwise; series rel 0 |
| EMHD | bitwise; series rel 2.7e-43 | bitwise; series rel 2.7e-43 |

- In all 12 A cases, out.std prints `Manually dt is changed from 2.0000E-003 to
  1.0000E-003` once, and the `dt_adjust` file is deleted. Final tt is
  0.20000000000000015 in A and B.
- EMHD: the only series difference is `wmag_dot_sum` at t = 0, a round-off zero
  divided by dt (1.2e-46 in A, still at dt = 2e-3; 2.4e-46 in B; series scale 4.5e-4).
- Binaries without the dt_adjust fix, run as case A against the rc case B (info):
  v1.0 MHD_INCOMP restart max rel 9.3e-3, series 1.3e-2; v1.0 RMHD 1.8e-3, 1.4e-2;
  KRMHD f05d2c1 6.3e-3, 5.1e-2. All diverge from t = 0.009, the first output after
  t = 0: the host dt changes but the device keeps the old dt.
- No R5 case wrote a 3D record (every `out3d/*.dat` is empty), so `mpiio_write_var`
  never ran.

## Issues found during the checks

1. **3D output: an all-zero gap after every record** (all models, all
   decompositions). **Fixed in a93e71c**; R4 passes on `9c01367`.
   In `f529a53`, `src/mpiio_write_var.F90:26-29` advances `disp` by
   `sizes(1)*sizes(2)*sizes(3)*bytes` and then, for `complex_type`, by the same
   amount again. `bytes` is `MPI_TYPE_SIZE(MPI_DOUBLE_COMPLEX)` = 16, which already
   counts both parts; the 2DECOMP&FFT original multiplied by the size of the real
   type. A file with N records is therefore 2N − 1 records long. A reader that takes
   record i at offset i × (record size) reads zeros or the wrong time, and
   `misc/backup_output.py` (`merge_bindir`), which assumes size / len(time.dat) bytes
   per record, cuts mid-record when it trims overlapping chunks. The code has been
   unchanged since v1.0 but never ran: the `loop_io_3D` calls were commented out
   until bfa0e63. a93e71c deletes lines 27-29. `mpiio_write_var` is called only by
   `loop_io_3D`, so restart files are not affected.
2. **KRMHD `init_type = 'OT2'` starts with omg = 0.** **Fixed in 9c01367**; R3′ is
   bitwise on `9c01367`. In `f529a53`,
   `src/model/KRMHD/fields.F90:444-453` (`init_OT2`) sets phi and psi but not omg;
   the RMHD version sets `omg = -kprp2*phi` (`src/model/RMHD/fields.F90:344`). The
   first step rebuilds phi from omg, so the initial flow is lost. `init_single_mode`
   (`:160`) and `init_random` (`:287`) set omg. The shipped KRMHD inputs
   (`tests/KRMHD/gallope.in`, `input_example/KRMHD.in`) use `init_type = 'zero'` and
   are not affected.
3. **v1.0 does not build with nvhpc ≥ 26.1** (`NVFORTRAN-S-1295` at
   `src/model/MHD_INCOMP/diagnostics.F90:125` and `src/model/RMHD/diagnostics.F90:84`
   of v1.0: `allocate(src(...), source=(0.d0, 0.d0))` with a real `src`). Fixed in
   a0f1b7a; `tests/release/v1.0_nvhpc26.patch` applies the same change to a v1.0
   tree for the comparisons above. It only changes the zero used for the allocation.
4. **3D output is new in v1.1.** v1.0 + 3D calls on 1 GPU writes the same records as
   v1.1, with the zero gaps of Issue 1 between them. v1.0 itself never called `loop_io_3D` (`src/gallope.F90:61, 73, 83,
   110` commented out since 19999d4), and with the calls enabled it writes wrong
   files on more than one GPU (fixed in 2bbed36).

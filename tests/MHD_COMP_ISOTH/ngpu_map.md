# ngpu of the MHD_COMP_ISOTH cases on TSUBAME4.0

Written by `make_cases.py ngpu-map` (rule `ngpu_rule`); a file `ngpu` in a run directory overrides it.
A case runs on 1 GPU unless it is a member of a decomposition test (`DECOMP`): the reference of a family
runs on 1 GPU and the other members on 2, 4, ... GPUs in the order of their Calliope np (the j-th
distinct np gets 2^j), capped at 4 and at the largest power of two that divides nlx and nly.

Cases: 589; ngpu 1: 546 cases, 2: 41 cases, 4: 2 cases.

## Euler1D

| case | Calliope np | ngpu | reason |
|---|---:|---:|---|
| stage1_es/cflmax/U2_N512_f1.0 | 1 | 1 | single GPU |
| stage1_es/cflmax/U2_N512_f1.2 | 1 | 1 | single GPU |
| stage1_es/cflmax/U2_N512_f1.6 | 1 | 1 | single GPU |
| stage1_es/cflmax/U2_N512_f2.0 | 1 | 1 | single GPU |
| stage1_es/cflmax/U2_N512_f2.4 | 1 | 1 | single GPU |
| stage1_es/cflmax/U2_N512_f3.0 | 1 | 1 | single GPU |
| stage1_es/decomp/U2_N512_x_np1 | 1 | 1 | reference of decomposition test stage1_es/decomp/U2_N512, judge_stage1.judge_batch_decomp |
| stage1_es/decomp/U2_N512_x_np2 | 2 | 2 | decomposition test stage1_es/decomp/U2_N512: np 2 is no. 2 of np 1/2/4 -> 2^2, capped at 2 (the ranks divide nlx = 512 and nly = 2) |
| stage1_es/decomp/U2_N512_x_np4 | 4 | 2 | decomposition test stage1_es/decomp/U2_N512: np 4 is no. 3 of np 1/2/4 -> 2^3, capped at 2 (the ranks divide nlx = 512 and nly = 2) |
| stage1_es/decomp/U2_N512_y_np1 | 1 | 2 | decomposition test stage1_es/decomp/U2_N512: np 1 is no. 1 of np 1/2/4 -> 2^1 |
| stage1_es/decomp/U2_N512_y_np2 | 2 | 2 | decomposition test stage1_es/decomp/U2_N512: np 2 is no. 2 of np 1/2/4 -> 2^2, capped at 2 (the ranks divide nlx = 2 and nly = 512) |
| stage1_es/decomp/U2_N512_y_np4 | 4 | 2 | decomposition test stage1_es/decomp/U2_N512: np 4 is no. 3 of np 1/2/4 -> 2^3, capped at 2 (the ranks divide nlx = 2 and nly = 512) |
| stage1_es/decomp/U2_N512_z_np1 | 1 | 2 | decomposition test stage1_es/decomp/U2_N512: np 1 is no. 1 of np 1/2/4 -> 2^1 |
| stage1_es/decomp/U2_N512_z_np2 | 2 | 2 | decomposition test stage1_es/decomp/U2_N512: np 2 is no. 2 of np 1/2/4 -> 2^2, capped at 2 (the ranks divide nlx = 2 and nly = 2) |
| stage1_es/decomp/U2_N512_z_np4 | 4 | 2 | decomposition test stage1_es/decomp/U2_N512: np 4 is no. 3 of np 1/2/4 -> 2^3, capped at 2 (the ranks divide nlx = 2 and nly = 2) |
| stage1_es/dtconv/U2_N512_dt1 | 1 | 1 | single GPU |
| stage1_es/dtconv/U2_N512_dt2 | 1 | 1 | single GPU |
| stage1_es/dtconv/U2_N512_dt4 | 1 | 1 | single GPU |
| stage1_es/dtconv/U2_N512_dt8 | 1 | 1 | single GPU |
| stage1_es/dtconv/U4_N2048_cflsvv1 | 4 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage1_es/fast/default | 1 | 1 | single GPU |
| stage1_es/fast/default_linear | 1 | 1 | single GPU |
| stage1_es/fast/ideal | 1 | 1 | single GPU |
| stage1_es/scan/Sod_N1024_C0 | 1 | 1 | single GPU |
| stage1_es/scan/Sod_N1024_C1 | 1 | 1 | single GPU |
| stage1_es/scan/Sod_N1024_C4 | 1 | 1 | single GPU |
| stage1_es/scan/Sod_N1024_eps0.25 | 1 | 1 | single GPU |
| stage1_es/scan/Sod_N1024_eps1 | 1 | 1 | single GPU |
| stage1_es/scan/Sod_N1024_fs | 1 | 1 | single GPU |
| stage1_es/scan/Sod_N1024_fs_mc | 1 | 1 | single GPU |
| stage1_es/scan/Sod_N1024_fs_nt6 | 1 | 1 | single GPU |
| stage1_es/scan/Sod_N1024_linear | 1 | 1 | single GPU |
| stage1_es/scan/Sod_N1024_mc | 1 | 1 | single GPU |
| stage1_es/scan/Sod_N1024_nosvv | 1 | 1 | single GPU |
| stage1_es/scan/Sod_N1024_pruned | 1 | 1 | single GPU |
| stage1_es/scan/Sod_N1024_svvQ0 | 1 | 1 | single GPU |
| stage1_es/scan/U2_N1024_C0 | 1 | 1 | single GPU |
| stage1_es/scan/U2_N1024_C1 | 1 | 1 | single GPU |
| stage1_es/scan/U2_N1024_C4 | 1 | 1 | single GPU |
| stage1_es/scan/U2_N1024_eps0.25 | 1 | 1 | single GPU |
| stage1_es/scan/U2_N1024_eps1 | 1 | 1 | single GPU |
| stage1_es/scan/U2_N1024_fs | 1 | 1 | single GPU |
| stage1_es/scan/U2_N1024_fs_mc | 1 | 1 | single GPU |
| stage1_es/scan/U2_N1024_fs_nt6 | 1 | 1 | single GPU |
| stage1_es/scan/U2_N1024_linear | 1 | 1 | single GPU |
| stage1_es/scan/U2_N1024_mc | 1 | 1 | single GPU |
| stage1_es/scan/U2_N1024_nosvv | 1 | 1 | single GPU |
| stage1_es/scan/U2_N1024_pruned | 1 | 1 | single GPU |
| stage1_es/scan/U2_N1024_svvQ0 | 1 | 1 | single GPU |
| stage1_es/scan/U4_N1024_C1 | 1 | 1 | single GPU |
| stage1_es/scan/U4_N1024_C4 | 1 | 1 | single GPU |
| stage1_es/scan/U4_N1024_linear | 1 | 1 | single GPU |
| stage1_es/scan/U4_N1024_nosvv | 1 | 1 | single GPU |
| stage1_es/scan/U4_N1024_svvQ0 | 1 | 1 | single GPU |
| stage1_es/shock/Sod_N1024 | 1 | 1 | single GPU |
| stage1_es/shock/Sod_N2048 | 1 | 1 | single GPU |
| stage1_es/shock/Sod_N256 | 1 | 1 | single GPU |
| stage1_es/shock/Sod_N4096 | 1 | 1 | single GPU |
| stage1_es/shock/Sod_N512 | 1 | 1 | single GPU |
| stage1_es/shock/U1_N1024 | 1 | 1 | single GPU |
| stage1_es/shock/U1_N2048 | 1 | 1 | single GPU |
| stage1_es/shock/U1_N256 | 1 | 1 | single GPU |
| stage1_es/shock/U1_N4096 | 1 | 1 | single GPU |
| stage1_es/shock/U1_N512 | 1 | 1 | single GPU |
| stage1_es/shock/U2_N1024 | 1 | 1 | single GPU |
| stage1_es/shock/U2_N2048 | 1 | 1 | single GPU |
| stage1_es/shock/U2_N256 | 1 | 1 | single GPU |
| stage1_es/shock/U2_N4096 | 2 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage1_es/shock/U2_N512 | 1 | 1 | single GPU |
| stage1_es/shock/U4_N1024 | 1 | 1 | single GPU |
| stage1_es/shock/U4_N2048 | 2 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage1_es/shock/U4_N256 | 1 | 1 | single GPU |
| stage1_es/shock/U4_N4096 | 4 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage1_es/shock/U4_N512 | 1 | 1 | single GPU |
| stage1_es/smoke/collide_U2_N512 | 1 | 1 | single GPU |
| stage1_es/smoke/collide_U2_N512_C1_nosvv_fs | 1 | 1 | single GPU |
| stage1_es/smoke/collide_U2_N512_est | 1 | 1 | single GPU |
| stage1_es/smoke/collide_U2_N512_est_chain | 1 | 1 | single GPU |
| stage1_es/smoke/nan | 2 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage1_es/smoke/selftest_np1 | 1 | 1 | reference of decomposition test stage1_es/smoke/selftest, judge_stage1 (shock_selftest.dat) |
| stage1_es/smoke/selftest_np2 | 2 | 2 | decomposition test stage1_es/smoke/selftest: np 2 is no. 1 of np 2/4 -> 2^1 |
| stage1_es/smoke/selftest_np4 | 4 | 4 | decomposition test stage1_es/smoke/selftest: np 4 is no. 2 of np 2/4 -> 2^2 |
| stage1_es/smoke/sound64 | 1 | 1 | single GPU |
| stage1_es/smooth/simple_N128_C0_nosvv | 1 | 1 | single GPU |
| stage1_es/smooth/simple_N128_C0_svv | 1 | 1 | single GPU |
| stage1_es/smooth/simple_N128_C2_svv | 1 | 1 | single GPU |
| stage1_es/smooth/simple_N256_C0_nosvv | 1 | 1 | single GPU |
| stage1_es/smooth/simple_N256_C0_svv | 1 | 1 | single GPU |
| stage1_es/smooth/simple_N256_C2_svv | 1 | 1 | single GPU |
| stage1_es/smooth/simple_N32_C0_nosvv | 1 | 1 | single GPU |
| stage1_es/smooth/simple_N32_C0_svv | 1 | 1 | single GPU |
| stage1_es/smooth/simple_N32_C2_svv | 1 | 1 | single GPU |
| stage1_es/smooth/simple_N512_C0_nosvv | 1 | 1 | single GPU |
| stage1_es/smooth/simple_N512_C0_svv | 1 | 1 | single GPU |
| stage1_es/smooth/simple_N512_C2_svv | 1 | 1 | single GPU |
| stage1_es/smooth/simple_N64_C0_nosvv | 1 | 1 | single GPU |
| stage1_es/smooth/simple_N64_C0_svv | 1 | 1 | single GPU |
| stage1_es/smooth/simple_N64_C2_svv | 1 | 1 | single GPU |
| stage1_es/smooth/sound_N32_C0 | 1 | 1 | single GPU |
| stage1_es/smooth/sound_N32_C2 | 1 | 1 | single GPU |
| stage1_es/smooth/sound_N64_C0 | 1 | 1 | single GPU |
| stage1_es/smooth/sound_N64_C2 | 1 | 1 | single GPU |

## MHD1D

| case | Calliope np | ngpu | reason |
|---|---:|---:|---|
| stage14/mhd1d/decomp/C3_N256_x_np1 | 1 | 1 | reference of decomposition test stage14/mhd1d/decomp/C3_N256, analyze_stage14 (TOL_DEC1D) |
| stage14/mhd1d/decomp/C3_N256_x_np2 | 2 | 2 | decomposition test stage14/mhd1d/decomp/C3_N256: np 2 is no. 1 of np 2/4 -> 2^1 |
| stage14/mhd1d/decomp/C3_N256_x_np4 | 4 | 2 | decomposition test stage14/mhd1d/decomp/C3_N256: np 4 is no. 2 of np 2/4 -> 2^2, capped at 2 (the ranks divide nlx = 256 and nly = 2) |
| stage14/mhd1d/decomp/C3_N256_y_np2 | 2 | 2 | decomposition test stage14/mhd1d/decomp/C3_N256: np 2 is no. 1 of np 2/4 -> 2^1 |
| stage14/mhd1d/decomp/C3_N256_z_np4 | 4 | 2 | decomposition test stage14/mhd1d/decomp/C3_N256: np 4 is no. 2 of np 2/4 -> 2^2, capped at 2 (the ranks divide nlx = 2 and nly = 2) |
| stage14/mhd1d/decomp/F1_N256_x_np1 | 1 | 1 | reference of decomposition test stage14/mhd1d/decomp/F1_N256, analyze_stage14 (TOL_DEC1D) |
| stage14/mhd1d/decomp/F1_N256_x_np2 | 2 | 2 | decomposition test stage14/mhd1d/decomp/F1_N256: np 2 is no. 1 of np 2/4 -> 2^1 |
| stage14/mhd1d/decomp/F1_N256_x_np4 | 4 | 2 | decomposition test stage14/mhd1d/decomp/F1_N256: np 4 is no. 2 of np 2/4 -> 2^2, capped at 2 (the ranks divide nlx = 256 and nly = 2) |
| stage14/mhd1d/decomp/F1_N256_y_np2 | 2 | 2 | decomposition test stage14/mhd1d/decomp/F1_N256: np 2 is no. 1 of np 2/4 -> 2^1 |
| stage14/mhd1d/decomp/F1_N256_z_np4 | 4 | 2 | decomposition test stage14/mhd1d/decomp/F1_N256: np 4 is no. 2 of np 2/4 -> 2^2, capped at 2 (the ranks divide nlx = 2 and nly = 2) |
| stage14/mhd1d/decomp/RDh_w4_N256_x_np1 | 1 | 1 | reference of decomposition test stage14/mhd1d/decomp/RDh_w4_N256, analyze_stage14 (TOL_DEC1D) |
| stage14/mhd1d/decomp/RDh_w4_N256_x_np2 | 2 | 2 | decomposition test stage14/mhd1d/decomp/RDh_w4_N256: np 2 is no. 1 of np 2/4 -> 2^1 |
| stage14/mhd1d/decomp/RDh_w4_N256_x_np4 | 4 | 2 | decomposition test stage14/mhd1d/decomp/RDh_w4_N256: np 4 is no. 2 of np 2/4 -> 2^2, capped at 2 (the ranks divide nlx = 256 and nly = 2) |
| stage14/mhd1d/decomp/RDh_w4_N256_y_np2 | 2 | 2 | decomposition test stage14/mhd1d/decomp/RDh_w4_N256: np 2 is no. 1 of np 2/4 -> 2^1 |
| stage14/mhd1d/decomp/RDh_w4_N256_z_np4 | 4 | 2 | decomposition test stage14/mhd1d/decomp/RDh_w4_N256: np 4 is no. 2 of np 2/4 -> 2^2, capped at 2 (the ranks divide nlx = 2 and nly = 2) |
| stage14/mhd1d/eta_decomp/C3_N256_eta2_x_np1 | 1 | 1 | reference of decomposition test stage14/mhd1d/eta_decomp/C3_N256_eta2, analyze_stage14 (TOL_DEC1D) |
| stage14/mhd1d/eta_decomp/C3_N256_eta2_x_np2 | 2 | 2 | decomposition test stage14/mhd1d/eta_decomp/C3_N256_eta2: np 2 is no. 1 of np 2/4 -> 2^1 |
| stage14/mhd1d/eta_decomp/C3_N256_eta2_x_np4 | 4 | 2 | decomposition test stage14/mhd1d/eta_decomp/C3_N256_eta2: np 4 is no. 2 of np 2/4 -> 2^2, capped at 2 (the ranks divide nlx = 256 and nly = 2) |
| stage14/mhd1d/eta_decomp/C3_N256_eta2_y_np2 | 2 | 2 | decomposition test stage14/mhd1d/eta_decomp/C3_N256_eta2: np 2 is no. 1 of np 2/4 -> 2^1 |
| stage14/mhd1d/eta_decomp/C3_N256_eta2_z_np4 | 4 | 2 | decomposition test stage14/mhd1d/eta_decomp/C3_N256_eta2: np 4 is no. 2 of np 2/4 -> 2^2, capped at 2 (the ranks divide nlx = 2 and nly = 2) |
| stage14/mhd1d/eta_decomp/I1_N256_eta2_x_np1 | 1 | 1 | reference of decomposition test stage14/mhd1d/eta_decomp/I1_N256_eta2, analyze_stage14 (TOL_DEC1D) |
| stage14/mhd1d/eta_decomp/I1_N256_eta2_x_np2 | 2 | 2 | decomposition test stage14/mhd1d/eta_decomp/I1_N256_eta2: np 2 is no. 1 of np 2/4 -> 2^1 |
| stage14/mhd1d/eta_decomp/I1_N256_eta2_x_np4 | 4 | 2 | decomposition test stage14/mhd1d/eta_decomp/I1_N256_eta2: np 4 is no. 2 of np 2/4 -> 2^2, capped at 2 (the ranks divide nlx = 256 and nly = 2) |
| stage14/mhd1d/eta_decomp/I1_N256_eta2_y_np2 | 2 | 2 | decomposition test stage14/mhd1d/eta_decomp/I1_N256_eta2: np 2 is no. 1 of np 2/4 -> 2^1 |
| stage14/mhd1d/eta_decomp/I1_N256_eta2_z_np4 | 4 | 2 | decomposition test stage14/mhd1d/eta_decomp/I1_N256_eta2: np 4 is no. 2 of np 2/4 -> 2^2, capped at 2 (the ranks divide nlx = 2 and nly = 2) |
| stage14/mhd1d/eta_suite/C1_N1024_eta2 | 1 | 1 | single GPU |
| stage14/mhd1d/eta_suite/C1_N2048_eta2 | 2 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/eta_suite/C1_N256_eta2 | 1 | 1 | single GPU |
| stage14/mhd1d/eta_suite/C1_N4096_eta2 | 4 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/eta_suite/C1_N512_eta2 | 1 | 1 | single GPU |
| stage14/mhd1d/eta_suite/C2_N1024_eta2 | 1 | 1 | single GPU |
| stage14/mhd1d/eta_suite/C2_N2048_eta2 | 2 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/eta_suite/C2_N256_eta2 | 1 | 1 | single GPU |
| stage14/mhd1d/eta_suite/C2_N4096_eta2 | 4 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/eta_suite/C2_N512_eta2 | 1 | 1 | single GPU |
| stage14/mhd1d/eta_suite/C3_N1024_eta2 | 1 | 1 | single GPU |
| stage14/mhd1d/eta_suite/C3_N2048_eta2 | 2 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/eta_suite/C3_N256_eta2 | 1 | 1 | single GPU |
| stage14/mhd1d/eta_suite/C3_N4096_eta2 | 4 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/eta_suite/C3_N512_eta2 | 1 | 1 | single GPU |
| stage14/mhd1d/eta_suite/F1_N1024_eta2 | 1 | 1 | single GPU |
| stage14/mhd1d/eta_suite/F1_N2048_eta2 | 2 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/eta_suite/F1_N256_eta2 | 1 | 1 | single GPU |
| stage14/mhd1d/eta_suite/F1_N512_eta2 | 1 | 1 | single GPU |
| stage14/mhd1d/eta_suite/F2_N1024_eta2 | 1 | 1 | single GPU |
| stage14/mhd1d/eta_suite/F2_N2048_eta2 | 2 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/eta_suite/F2_N256_eta2 | 1 | 1 | single GPU |
| stage14/mhd1d/eta_suite/F2_N512_eta2 | 1 | 1 | single GPU |
| stage14/mhd1d/eta_suite/I1_N1024_eta2 | 1 | 1 | single GPU |
| stage14/mhd1d/eta_suite/I1_N2048_eta2 | 2 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/eta_suite/I1_N256_eta2 | 1 | 1 | single GPU |
| stage14/mhd1d/eta_suite/I1_N512_eta2 | 1 | 1 | single GPU |
| stage14/mhd1d/eta_suite/P1b_N1024_eta2 | 1 | 1 | single GPU |
| stage14/mhd1d/eta_suite/P1b_N2048_eta2 | 2 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/eta_suite/P1b_N256_eta2 | 1 | 1 | single GPU |
| stage14/mhd1d/eta_suite/P1b_N512_eta2 | 1 | 1 | single GPU |
| stage14/mhd1d/eta_suite/P2_N1024_eta2 | 1 | 1 | single GPU |
| stage14/mhd1d/eta_suite/P2_N2048_eta2 | 2 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/eta_suite/P2_N256_eta2 | 1 | 1 | single GPU |
| stage14/mhd1d/eta_suite/P2_N512_eta2 | 1 | 1 | single GPU |
| stage14/mhd1d/eta_suite/RDh_w4_N1024_eta2 | 1 | 1 | single GPU |
| stage14/mhd1d/eta_suite/RDh_w4_N256_eta2 | 1 | 1 | single GPU |
| stage14/mhd1d/eta_suite/RDh_w4_N512_eta2 | 1 | 1 | single GPU |
| stage14/mhd1d/eta_suite/RDpc_N1024_eta2 | 1 | 1 | single GPU |
| stage14/mhd1d/eta_suite/RDpc_N256_eta2 | 1 | 1 | single GPU |
| stage14/mhd1d/eta_suite/RDpc_N512_eta2 | 1 | 1 | single GPU |
| stage14/mhd1d/eta_suite/S1_N1024_eta2 | 1 | 1 | single GPU |
| stage14/mhd1d/eta_suite/S1_N2048_eta2 | 2 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/eta_suite/S1_N256_eta2 | 1 | 1 | single GPU |
| stage14/mhd1d/eta_suite/S1_N512_eta2 | 1 | 1 | single GPU |
| stage14/mhd1d/eta_suite/S2_N1024_eta2 | 1 | 1 | single GPU |
| stage14/mhd1d/eta_suite/S2_N2048_eta2 | 2 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/eta_suite/S2_N256_eta2 | 1 | 1 | single GPU |
| stage14/mhd1d/eta_suite/S2_N512_eta2 | 1 | 1 | single GPU |
| stage14/mhd1d/eta_suite/TDS_w1_N256_eta2 | 1 | 1 | single GPU |
| stage14/mhd1d/eta_suite/TDS_w1_N512_eta2 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/C1_N1024 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/C1_N2048 | 2 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/suite/C1_N256 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/C1_N4096 | 4 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/suite/C1_N512 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/C2_N1024 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/C2_N2048 | 2 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/suite/C2_N256 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/C2_N4096 | 4 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/suite/C2_N512 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/C3_N1024 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/C3_N2048 | 2 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/suite/C3_N256 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/C3_N4096 | 4 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/suite/C3_N512 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/F1_N1024 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/F1_N2048 | 2 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/suite/F1_N256 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/F1_N512 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/F2_N1024 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/F2_N2048 | 2 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/suite/F2_N256 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/F2_N512 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/I1_N1024 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/I1_N2048 | 2 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/suite/I1_N256 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/I1_N512 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/P1a_N1024 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/P1a_N2048 | 2 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/suite/P1a_N256 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/P1a_N512 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/P1b_N1024 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/P1b_N2048 | 2 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/suite/P1b_N256 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/P1b_N512 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/P1c_N1024 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/P1c_N2048 | 2 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/suite/P1c_N256 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/P1c_N512 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/P2_N1024 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/P2_N2048 | 2 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/suite/P2_N256 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/P2_N512 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/P3_N1024 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/P3_N2048 | 2 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/suite/P3_N256 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/P3_N512 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/RDh_w0_N1024 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/RDh_w0_N2048 | 2 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/suite/RDh_w0_N256 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/RDh_w0_N512 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/RDh_w4_N1024 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/RDh_w4_N2048 | 2 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/suite/RDh_w4_N256 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/RDh_w4_N512 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/RDh_w8_N1024 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/RDh_w8_N2048 | 2 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/suite/RDh_w8_N256 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/RDh_w8_N512 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/RDp_w0_N1024 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/RDp_w0_N2048 | 2 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/suite/RDp_w0_N256 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/RDp_w0_N512 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/RDp_w4_N1024 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/RDp_w4_N2048 | 2 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/suite/RDp_w4_N256 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/RDp_w4_N512 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/RDp_w8_N1024 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/RDp_w8_N2048 | 2 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/suite/RDp_w8_N256 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/RDp_w8_N512 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/RDpc_N1024 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/RDpc_N2048 | 2 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/suite/RDpc_N256 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/RDpc_N512 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/S1_N1024 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/S1_N2048 | 2 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/suite/S1_N256 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/S1_N512 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/S2_N1024 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/S2_N2048 | 2 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/suite/S2_N256 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/S2_N512 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/TDA_w0_N1024 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/TDA_w0_N2048 | 2 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/suite/TDA_w0_N256 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/TDA_w0_N512 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/TDA_w0pre_N1024 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/TDA_w0pre_N2048 | 2 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/suite/TDA_w0pre_N256 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/TDA_w0pre_N512 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/TDA_w1_N1024 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/TDA_w1_N2048 | 2 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/suite/TDA_w1_N256 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/TDA_w1_N512 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/TDA_w1pre_N1024 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/TDA_w1pre_N2048 | 2 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/suite/TDA_w1pre_N256 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/TDA_w1pre_N512 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/TDA_w4_N1024 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/TDA_w4_N2048 | 2 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/suite/TDA_w4_N256 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/TDA_w4_N512 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/TDA_w4pre_N1024 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/TDA_w4pre_N2048 | 2 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/suite/TDA_w4pre_N256 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/TDA_w4pre_N512 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/TDS_w0_N1024 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/TDS_w0_N2048 | 2 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/suite/TDS_w0_N256 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/TDS_w0_N512 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/TDS_w0pre_N1024 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/TDS_w0pre_N2048 | 2 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/suite/TDS_w0pre_N256 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/TDS_w0pre_N512 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/TDS_w1_N1024 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/TDS_w1_N2048 | 2 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/suite/TDS_w1_N256 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/TDS_w1_N512 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/TDS_w1pre_N1024 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/TDS_w1pre_N2048 | 2 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/suite/TDS_w1pre_N256 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/TDS_w1pre_N512 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/TDS_w4_N1024 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/TDS_w4_N2048 | 2 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/suite/TDS_w4_N256 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/TDS_w4_N512 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/TDS_w4pre_N1024 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/TDS_w4pre_N2048 | 2 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/mhd1d/suite/TDS_w4pre_N256 | 1 | 1 | single GPU |
| stage14/mhd1d/suite/TDS_w4pre_N512 | 1 | 1 | single GPU |
| stage16/conv/I1_N4096 | 4 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage16/conv/I1_N4096_eta2 | 4 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage16/conv/I1_N8192 | 8 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage16/conv/I1_N8192_eta2 | 8 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage16/conv/RDpc_N4096 | 4 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage16/conv/RDpc_N4096_eta2 | 4 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage16/conv/RDpc_N8192 | 8 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage16/conv/RDpc_N8192_eta2 | 8 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |

## Smooth_wave

| case | Calliope np | ngpu | reason |
|---|---:|---:|---|
| stage4/alfven/2d_k101_a01_N16_default_p1600 | 1 | 1 | single GPU |
| stage4/alfven/2d_k101_a01_N16_default_p400 | 1 | 1 | single GPU |
| stage4/alfven/2d_k101_a01_N16_default_p800 | 1 | 1 | single GPU |
| stage4/alfven/2d_k101_a01_N16_ideal_p1600 | 1 | 1 | single GPU |
| stage4/alfven/2d_k101_a01_N16_ideal_p400 | 1 | 1 | single GPU |
| stage4/alfven/2d_k101_a01_N16_ideal_p800 | 1 | 1 | single GPU |
| stage4/alfven/2d_k101_a01_N32_default_p400 | 1 | 1 | single GPU |
| stage4/alfven/2d_k101_a01_N32_ideal_p400 | 1 | 1 | single GPU |
| stage4/alfven/2d_k101_a01_N64_default_p400 | 1 | 1 | single GPU |
| stage4/alfven/2d_k101_a01_N64_ideal_p400 | 1 | 1 | single GPU |
| stage4/alfven/2d_k101_a1_N16_default_p1600 | 1 | 1 | single GPU |
| stage4/alfven/2d_k101_a1_N16_default_p3200 | 1 | 1 | single GPU |
| stage4/alfven/2d_k101_a1_N16_default_p800 | 1 | 1 | single GPU |
| stage4/alfven/2d_k101_a1_N16_ideal_p1600 | 1 | 1 | single GPU |
| stage4/alfven/2d_k101_a1_N16_ideal_p3200 | 1 | 1 | single GPU |
| stage4/alfven/2d_k101_a1_N16_ideal_p800 | 1 | 1 | single GPU |
| stage4/alfven/2d_k101_a1_N32_default_p800 | 1 | 1 | reference of decomposition test stage4/alfven/2d_k101_a1_N32_default_p800, analyze_stage4 (TOL_MPI) |
| stage4/alfven/2d_k101_a1_N32_default_p800_np16 | 16 | 2 | decomposition test stage4/alfven/2d_k101_a1_N32_default_p800: np 16 is no. 2 of np 4/16 -> 2^2, capped at 2 (the ranks divide nlx = 32 and nly = 2) |
| stage4/alfven/2d_k101_a1_N32_default_p800_np4 | 4 | 2 | decomposition test stage4/alfven/2d_k101_a1_N32_default_p800: np 4 is no. 1 of np 4/16 -> 2^1 |
| stage4/alfven/2d_k101_a1_N32_ideal_p800 | 1 | 1 | single GPU |
| stage4/alfven/2d_k101_a1_N64_default_p800 | 1 | 1 | single GPU |
| stage4/alfven/2d_k101_a1_N64_ideal_p800 | 1 | 1 | single GPU |
| stage4/alfven/2d_k102_a01_N16_default_p200 | 1 | 1 | single GPU |
| stage4/alfven/2d_k102_a01_N16_default_p400 | 1 | 1 | single GPU |
| stage4/alfven/2d_k102_a01_N16_default_p800 | 1 | 1 | single GPU |
| stage4/alfven/2d_k102_a01_N16_ideal_p200 | 1 | 1 | single GPU |
| stage4/alfven/2d_k102_a01_N16_ideal_p400 | 1 | 1 | single GPU |
| stage4/alfven/2d_k102_a01_N16_ideal_p800 | 1 | 1 | single GPU |
| stage4/alfven/2d_k102_a01_N32_default_p200 | 1 | 1 | single GPU |
| stage4/alfven/2d_k102_a01_N32_ideal_p200 | 1 | 1 | single GPU |
| stage4/alfven/2d_k102_a01_N64_default_p200 | 1 | 1 | single GPU |
| stage4/alfven/2d_k102_a01_N64_ideal_p200 | 1 | 1 | single GPU |
| stage4/alfven/2d_k102_a1_N16_default_p1600 | 1 | 1 | single GPU |
| stage4/alfven/2d_k102_a1_N16_default_p400 | 1 | 1 | single GPU |
| stage4/alfven/2d_k102_a1_N16_default_p800 | 1 | 1 | single GPU |
| stage4/alfven/2d_k102_a1_N16_ideal_p1600 | 1 | 1 | single GPU |
| stage4/alfven/2d_k102_a1_N16_ideal_p400 | 1 | 1 | single GPU |
| stage4/alfven/2d_k102_a1_N16_ideal_p800 | 1 | 1 | single GPU |
| stage4/alfven/2d_k102_a1_N32_default_p400 | 1 | 1 | single GPU |
| stage4/alfven/2d_k102_a1_N32_ideal_p400 | 1 | 1 | single GPU |
| stage4/alfven/2d_k102_a1_N64_default_p400 | 1 | 1 | single GPU |
| stage4/alfven/2d_k102_a1_N64_ideal_p400 | 1 | 1 | single GPU |
| stage4/alfven/3d_k112_a01_N16_default_p200 | 1 | 1 | single GPU |
| stage4/alfven/3d_k112_a01_N16_default_p400 | 1 | 1 | single GPU |
| stage4/alfven/3d_k112_a01_N16_default_p800 | 1 | 1 | single GPU |
| stage4/alfven/3d_k112_a01_N16_ideal_p200 | 1 | 1 | single GPU |
| stage4/alfven/3d_k112_a01_N16_ideal_p400 | 1 | 1 | single GPU |
| stage4/alfven/3d_k112_a01_N16_ideal_p800 | 1 | 1 | single GPU |
| stage4/alfven/3d_k112_a01_N32_default_p200_np8 | 8 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/alfven/3d_k112_a01_N32_ideal_p200_np8 | 8 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/alfven/3d_k112_a1_N16_default_p1600 | 1 | 1 | single GPU |
| stage4/alfven/3d_k112_a1_N16_default_p400 | 1 | 1 | reference of decomposition test stage4/alfven/3d_k112_a1_N16_default_p400, analyze_stage4 (TOL_MPI) |
| stage4/alfven/3d_k112_a1_N16_default_p400_np8 | 8 | 2 | decomposition test stage4/alfven/3d_k112_a1_N16_default_p400: np 8 is no. 1 of np 8 -> 2^1 |
| stage4/alfven/3d_k112_a1_N16_default_p800 | 1 | 1 | single GPU |
| stage4/alfven/3d_k112_a1_N16_ideal_p1600 | 1 | 1 | single GPU |
| stage4/alfven/3d_k112_a1_N16_ideal_p400 | 1 | 1 | single GPU |
| stage4/alfven/3d_k112_a1_N16_ideal_p800 | 1 | 1 | single GPU |
| stage4/alfven/3d_k112_a1_N32_default_p400_np8 | 8 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/alfven/3d_k112_a1_N32_ideal_p400_np8 | 8 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/alfven/3d_k121_a1_N16_default_p400 | 1 | 1 | single GPU |
| stage4/alfven/3d_k121_a1_N16_ideal_p400 | 1 | 1 | single GPU |
| stage4/alfven/3d_k211_a1_N16_default_p400 | 1 | 1 | single GPU |
| stage4/alfven/3d_k211_a1_N16_ideal_p400 | 1 | 1 | single GPU |
| stage4/cpdiag/2d_k101_a1_N16_etaonly_p800 | 1 | 1 | single GPU |
| stage4/cpdiag/2d_k101_a1_N16_noeta_p800 | 1 | 1 | single GPU |
| stage4/cpdiag/2d_k101_a1_N16_nosvv_p800 | 1 | 1 | single GPU |
| stage4/cpdiag/2d_k101_a1_N16_svvonly_p800 | 1 | 1 | single GPU |
| stage4/cpdiag/2d_k102_a1_N16_etaonly_p400 | 1 | 1 | single GPU |
| stage4/cpdiag/2d_k102_a1_N16_noeta_p400 | 1 | 1 | single GPU |
| stage4/cpdiag/2d_k102_a1_N16_nosvv_p400 | 1 | 1 | single GPU |
| stage4/cpdiag/2d_k102_a1_N16_svvonly_p400 | 1 | 1 | single GPU |
| stage4/cpdiag/3d_k112_a1_N16_etaonly_p400 | 1 | 1 | single GPU |
| stage4/cpdiag/3d_k112_a1_N16_noeta_p400 | 1 | 1 | single GPU |
| stage4/cpdiag/3d_k112_a1_N16_nosvv_p400 | 1 | 1 | single GPU |
| stage4/cpdiag/3d_k112_a1_N16_svvonly_p400 | 1 | 1 | single GPU |
| stage4/cpdiag/3d_k112_a1_N32_etaonly_p400_np8 | 8 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/cpdiag/3d_k112_a1_N32_noeta_p400_np8 | 8 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/cpdiag/3d_k112_a1_N32_nosvv_p400_np8 | 8 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/cpdiag/3d_k112_a1_N32_svvonly_p400_np8 | 8 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/planar/ctl_C2_eta2_N256_1d | 1 | 1 | single GPU |
| stage4/planar/ctl_C2_eta2_N256_2d | 16 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/planar/ctl_F1_N256_1d | 1 | 1 | single GPU |
| stage4/planar/ctl_F1_N256_2d | 16 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/planar/ctl_RDh_w4_N256_1d | 1 | 1 | single GPU |
| stage4/planar/ctl_RDh_w4_N256_2d | 16 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/planar/ctl_TDS_w1_N256_1d | 1 | 1 | single GPU |
| stage4/planar/ctl_TDS_w1_N256_2d | 16 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/planar/eqv_C2_eta2_k11_N256_1d | 1 | 1 | single GPU |
| stage4/planar/eqv_C2_eta2_k11_N256_2d | 16 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/planar/eqv_F1_eta2_k11_N256_1d | 1 | 1 | single GPU |
| stage4/planar/eqv_F1_eta2_k11_N256_2d | 16 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/planar/eqv_F1_k11_N256_1d | 1 | 1 | single GPU |
| stage4/planar/eqv_F1_k11_N256_2d | 16 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/planar/eqv_RDh_w4_k11_N256_1d | 1 | 1 | single GPU |
| stage4/planar/eqv_RDh_w4_k11_N256_2d | 16 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/planar/eqv_TDS_w1_k11_N256_1d | 1 | 1 | single GPU |
| stage4/planar/eqv_TDS_w1_k11_N256_2d | 16 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/planar/phy_C2_eta2_k11_N256_1d | 1 | 1 | single GPU |
| stage4/planar/phy_C2_eta2_k11_N256_2d | 16 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/planar/phy_C2_eta2_k11_N512_1d | 1 | 1 | single GPU |
| stage4/planar/phy_C2_eta2_k11_N512_2d | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/planar/phy_C2_eta2_k12_N256_1d | 1 | 1 | single GPU |
| stage4/planar/phy_C2_eta2_k12_N256_1dh | 1 | 1 | single GPU |
| stage4/planar/phy_C2_eta2_k12_N256_2d | 16 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/planar/phy_C2_eta2_k12_N512_1d | 1 | 1 | single GPU |
| stage4/planar/phy_C2_eta2_k12_N512_1dh | 1 | 1 | single GPU |
| stage4/planar/phy_C2_eta2_k12_N512_2d | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/planar/phy_C2_k11_N256_1d | 1 | 1 | single GPU |
| stage4/planar/phy_C2_k11_N256_2d | 16 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/planar/phy_C2_k11_N512_1d | 1 | 1 | single GPU |
| stage4/planar/phy_C2_k11_N512_2d | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/planar/phy_C2_k12_N256_1d | 1 | 1 | single GPU |
| stage4/planar/phy_C2_k12_N256_1dh | 1 | 1 | single GPU |
| stage4/planar/phy_C2_k12_N256_2d | 16 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/planar/phy_C2_k12_N512_1d | 1 | 1 | single GPU |
| stage4/planar/phy_C2_k12_N512_1dh | 1 | 1 | single GPU |
| stage4/planar/phy_C2_k12_N512_2d | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/planar/phy_F1_eta2_k11_N256_1d | 1 | 1 | single GPU |
| stage4/planar/phy_F1_eta2_k11_N256_2d | 16 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/planar/phy_F1_eta2_k11_N512_1d | 1 | 1 | single GPU |
| stage4/planar/phy_F1_eta2_k11_N512_2d | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/planar/phy_F1_eta2_k12_N256_1d | 1 | 1 | single GPU |
| stage4/planar/phy_F1_eta2_k12_N256_1dh | 1 | 1 | single GPU |
| stage4/planar/phy_F1_eta2_k12_N256_2d | 16 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/planar/phy_F1_eta2_k12_N512_1d | 1 | 1 | single GPU |
| stage4/planar/phy_F1_eta2_k12_N512_1dh | 1 | 1 | single GPU |
| stage4/planar/phy_F1_eta2_k12_N512_2d | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/planar/phy_F1_k11_N256_1d | 1 | 1 | single GPU |
| stage4/planar/phy_F1_k11_N256_2d | 16 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/planar/phy_F1_k11_N512_1d | 1 | 1 | single GPU |
| stage4/planar/phy_F1_k11_N512_2d | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/planar/phy_F1_k12_N256_1d | 1 | 1 | single GPU |
| stage4/planar/phy_F1_k12_N256_1dh | 1 | 1 | single GPU |
| stage4/planar/phy_F1_k12_N256_2d | 16 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/planar/phy_F1_k12_N512_1d | 1 | 1 | single GPU |
| stage4/planar/phy_F1_k12_N512_1dh | 1 | 1 | single GPU |
| stage4/planar/phy_F1_k12_N512_2d | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/planar/phy_RDh_w4_k11_N256_1d | 1 | 1 | single GPU |
| stage4/planar/phy_RDh_w4_k11_N256_2d | 16 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/planar/phy_RDh_w4_k11_N512_1d | 1 | 1 | single GPU |
| stage4/planar/phy_RDh_w4_k11_N512_2d | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/planar/phy_RDh_w4_k12_N256_1d | 1 | 1 | single GPU |
| stage4/planar/phy_RDh_w4_k12_N256_1dh | 1 | 1 | single GPU |
| stage4/planar/phy_RDh_w4_k12_N256_2d | 16 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/planar/phy_RDh_w4_k12_N512_1d | 1 | 1 | single GPU |
| stage4/planar/phy_RDh_w4_k12_N512_1dh | 1 | 1 | single GPU |
| stage4/planar/phy_RDh_w4_k12_N512_2d | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/planar/phy_RDpc_k11_N256_1d | 1 | 1 | single GPU |
| stage4/planar/phy_RDpc_k11_N256_2d | 16 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/planar/phy_RDpc_k11_N512_1d | 1 | 1 | single GPU |
| stage4/planar/phy_RDpc_k11_N512_2d | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/planar/phy_RDpc_k12_N256_1d | 1 | 1 | single GPU |
| stage4/planar/phy_RDpc_k12_N256_1dh | 1 | 1 | single GPU |
| stage4/planar/phy_RDpc_k12_N256_2d | 16 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/planar/phy_RDpc_k12_N512_1d | 1 | 1 | single GPU |
| stage4/planar/phy_RDpc_k12_N512_1dh | 1 | 1 | single GPU |
| stage4/planar/phy_RDpc_k12_N512_2d | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/planar/phy_TDS_w1_k11_N256_1d | 1 | 1 | single GPU |
| stage4/planar/phy_TDS_w1_k11_N256_2d | 16 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/planar/phy_TDS_w1_k11_N512_1d | 1 | 1 | single GPU |
| stage4/planar/phy_TDS_w1_k11_N512_2d | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/planar/phy_TDS_w1_k12_N256_1d | 1 | 1 | single GPU |
| stage4/planar/phy_TDS_w1_k12_N256_1dh | 1 | 1 | single GPU |
| stage4/planar/phy_TDS_w1_k12_N256_2d | 16 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/planar/phy_TDS_w1_k12_N512_1d | 1 | 1 | single GPU |
| stage4/planar/phy_TDS_w1_k12_N512_1dh | 1 | 1 | single GPU |
| stage4/planar/phy_TDS_w1_k12_N512_2d | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |

## OT2D

| case | Calliope np | ngpu | reason |
|---|---:|---:|---|
| stage14/ot/dec_N256_np16 | 16 | 1 | reference of decomposition test stage14/ot/dec_N256, analyze_stage4 (ot_decomp) |
| stage14/ot/dec_N256_np16_n10 | 16 | 1 | reference of decomposition test stage14/ot/dec_N256_n10, analyze_stage4 (ot_decomp) |
| stage14/ot/dec_N256_np256 | 256 | 2 | decomposition test stage14/ot/dec_N256: np 256 is no. 3 of np 4/64/256 -> 2^3, capped at 2 (the ranks divide nlx = 256 and nly = 2) |
| stage14/ot/dec_N256_np256_n10 | 256 | 2 | decomposition test stage14/ot/dec_N256_n10: np 256 is no. 3 of np 4/64/256 -> 2^3, capped at 2 (the ranks divide nlx = 256 and nly = 2) |
| stage14/ot/dec_N256_np4 | 4 | 2 | decomposition test stage14/ot/dec_N256: np 4 is no. 1 of np 4/64/256 -> 2^1 |
| stage14/ot/dec_N256_np4_n10 | 4 | 2 | decomposition test stage14/ot/dec_N256_n10: np 4 is no. 1 of np 4/64/256 -> 2^1 |
| stage14/ot/dec_N256_np64 | 64 | 2 | decomposition test stage14/ot/dec_N256: np 64 is no. 2 of np 4/64/256 -> 2^2, capped at 2 (the ranks divide nlx = 256 and nly = 2) |
| stage14/ot/dec_N256_np64_n10 | 64 | 2 | decomposition test stage14/ot/dec_N256_n10: np 64 is no. 2 of np 4/64/256 -> 2^2, capped at 2 (the ranks divide nlx = 256 and nly = 2) |
| stage14/ot/ot_N1024_eta0_p256 | 256 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/ot/ot_N128_eta0 | 4 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/ot/ot_N128_eta2 | 4 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/ot/ot_N256_eta0 | 16 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/ot/ot_N256_eta2 | 16 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/ot/ot_N512_eta0 | 128 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage14/ot/ot_N512_eta2 | 128 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |

## Turb3D

| case | Calliope np | ngpu | reason |
|---|---:|---:|---|
| stage13d/nomass/e0.5_C2 | 128 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/a/turb05_N64_ES25 | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/a/turb05_N64_ESf | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/a/turb05_N64_H_r0.25 | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/a/turb05_N64_H_r0.5 | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/a/turb05_N64_H_r1 | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/a/turb05_N64_H_r16 | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/a/turb05_N64_H_r2 | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/a/turb05_N64_H_r4 | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/a/turb05_N64_H_r8 | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/a/turb05_N64_L_r0.25 | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/a/turb05_N64_L_r0.5 | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/a/turb05_N64_L_r1 | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/a/turb05_N64_L_r2 | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/a/turb05_N64_L_r4 | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/a/turb05_N64_R0 | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/a/turb1_N64_ES25 | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/a/turb1_N64_ESf | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/a/turb1_N64_H_r0.25 | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/a/turb1_N64_H_r0.5 | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/a/turb1_N64_H_r1 | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/a/turb1_N64_H_r16 | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/a/turb1_N64_H_r2 | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/a/turb1_N64_H_r4 | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/a/turb1_N64_H_r8 | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/a/turb1_N64_L_r0.25 | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/a/turb1_N64_L_r0.5 | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/a/turb1_N64_L_r1 | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/a/turb1_N64_L_r2 | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/a/turb1_N64_L_r4 | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/a/turb1_N64_R0 | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/b/turb05_N128_ES | 128 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/b/turb05_N128_H_r0.25 | 128 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/b/turb05_N128_H_r0.5 | 128 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/b/turb05_N128_H_r1 | 128 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/b/turb05_N128_L_r1 | 128 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/b/turb05_N128_L_r2 | 128 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/b/turb05_N128_R0 | 128 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/b/turb1_N128_ES | 128 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/b/turb1_N128_ESf | 128 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/b/turb1_N128_H_r2 | 128 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/b/turb1_N128_H_r4 | 128 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/b/turb1_N128_H_r8 | 128 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/b/turb1_N128_L_r1 | 128 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/b/turb1_N128_L_r2 | 128 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/c/turb1_N256_H_r32 | 256 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/c/turb1_N256_H_r8 | 256 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/pre/dec_turb1_N64_H_r1_np64 | 64 | 1 | reference of decomposition test stage15/pre/dec_turb1_N64_H_r1, analyze_stage15 (P0B_NP) |
| stage15/pre/dec_turb1_N64_H_r1_np8 | 8 | 2 | decomposition test stage15/pre/dec_turb1_N64_H_r1: np 8 is no. 1 of np 8 -> 2^1 |
| stage15/pre/dec_turb1_N64_L_r1_np64 | 64 | 1 | reference of decomposition test stage15/pre/dec_turb1_N64_L_r1, analyze_stage15 (P0B_NP) |
| stage15/pre/dec_turb1_N64_L_r1_np8 | 8 | 2 | decomposition test stage15/pre/dec_turb1_N64_L_r1: np 8 is no. 1 of np 8 -> 2^1 |
| stage15/pre/dt1_turb1_N64_ES | 16 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/pre/dt1_turb1_N64_H_r1 | 16 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/pre/dt1_turb1_N64_L_r1 | 16 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/pre/dt2_turb1_N64_ES | 16 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/pre/dt2_turb1_N64_H_r1 | 16 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/pre/dt2_turb1_N64_L_r1 | 16 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/pre/turb05_N64_ES | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage15/pre/turb1_N64_ES | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage16/hf128/turb1_N128_HF_r2 | 128 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage16/hf128/turb1_N128_HF_r4 | 128 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage16/hf128/turb1_N128_HF_r8 | 128 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage16/hf256/turb1_N256_HF_r8 | 256 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage16/hf64/turb1_N64_HF_r0.25 | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage16/hf64/turb1_N64_HF_r0.5 | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage16/hf64/turb1_N64_HF_r1 | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage16/hf64/turb1_N64_HF_r16 | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage16/hf64/turb1_N64_HF_r2 | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage16/hf64/turb1_N64_HF_r4 | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage16/hf64/turb1_N64_HF_r8 | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/turb/turb1_N128_a | 128 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/turb/turb1_N64 | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/turb/turb4_N64 | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/turb/turb4_N64_a | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage4/turb/turb4b_N64 | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage5/ms/turb2b_N64 | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage5/ms/turb3b_N64 | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage5/res/turb4_N128_a | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage5/res/turb4b_N128 | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage6q/track/dec_N64_np1 | 1 | 2 | decomposition test stage6q/track/dec_N64: np 1 is no. 1 of np 1/8 -> 2^1 |
| stage6q/track/dec_N64_np64 | 64 | 1 | reference of decomposition test stage6q/track/dec_N64, analyze_stage6 (TDEC_REF) |
| stage6q/track/dec_N64_np8 | 8 | 4 | decomposition test stage6q/track/dec_N64: np 8 is no. 2 of np 1/8 -> 2^2 |
| stage6q/track/turb1_N128_t | 128 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage6q/track/turb1_N64_t | 64 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage6q/turb/turb1_N256 | 256 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage8/coef/e0.25_C2 | 128 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage8/coef/e0.5_C1 | 128 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage8/coef/e0.5_C2 | 128 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage8/coef/e0.5_C4 | 128 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |
| stage8/coef/e1.0_C2 | 128 | 1 | single GPU (Calliope np > 1 is not a decomposition test) |

#!/bin/bash
#=======================================================================
# Create the v1.1 release-check run directories under runs/rel/ and
# print the CASES string for tests/run_cases-tsubame.sh. Run from the
# root of the release-candidate tree on TSUBAME4.0.
#
#   CASES=$(./tests/release/setup_runs.sh r12) RESOURCE=node_f H_RT=0:40:00 \
#     NAME=rel12 ./tests/run_cases-tsubame.sh
#
# Binaries (bin/gallope-<name>, the MODEL field of each case):
#   <MODEL>              release candidate (9c01367; r12 and r5 were run with
#                        f529a53, see validation_tsubame.md)
#   <MODEL>_v1.0         v1.0 (MHD_INCOMP, RMHD)
#   <MODEL>_v1.0_io3d    v1.0 + tests/release/enable_io3d.patch (MHD_INCOMP,
#                        RMHD): v1.0 never called loop_io_3D
#   KRMHD_f05d2c1        KRMHD at f05d2c1, before the dt_adjust fix of bab2b91
#
# Sets:
#   r12  R1  undriven, v1.0 vs v1.1, 1 and 4 GPUs:
#            r1_mi_* (MHD_INCOMP OT3), r1_ms_* (MHD_INCOMP OT3, shear q = 1.5,
#            dt = 5e-3, t = 1 > tremap = 2/3), r1_rm_* (RMHD OT3)
#        R2  driven (stir_seed fixed), v1.0 vs v1.1, 1 and 4 GPUs:
#            r2_mi_* (MHD_INCOMP, kz in [1, 2]), r2_rm_* (RMHD, kz in [1, 2]),
#            r2_rk0_* (RMHD, kmin(3) = 0: the kz = 0 plane is in the band)
#   r3   R3  KRMHD: r3_kf_* driven OT2 + g1 (np1, np4, np2 = P_fft 2,
#            pm2 = P_fft 2 x P_m 2, pm4 = P_m 4); r3_kr_* undriven random
#            (rs1 -> rs2 = 100 + 100 steps vs full; np2 vs pm2)
#        R3' KRMHD vs RMHD with one input: r3p_{k,r}_ot2 undriven OT2,
#            r3p_{k,r}_els Elsasser-driven without g
#   r4   R4  3D output, write_intvl_3D = 0.08 (records at t = 0, ~0.08, ~0.16
#            and the end, 0.2): np1 vs np4 for all six models; np2 vs pm2
#            (P_fft = 2, P_m = 1 vs 2) for RMHD, KRMHD, HRMHD, RRMHD; v1.0
#            (patched) np1 vs np4 for MHD_INCOMP and RMHD
#   r5   R5  dt_adjust on 1 GPU: A = dt 2e-3 in the input + file dt_adjust
#            (1.0d-3) in the run directory from the start; B = dt 1e-3.
#            Undriven; MHD_INCOMP, RMHD, KRMHD, HRMHD, EMHD, RRMHD with
#            eSSPIFRK3 and gear3 (_g3). A also with the v1.0 binaries and
#            KRMHD_f05d2c1 (eSSPIFRK3).
#=======================================================================
set -e
set_name=${1:?usage: setup_runs.sh r12|r3|r4|r5}
R=tests/release
cases=""

mk() {  # mk <case> <nranks> <binary> <base input> <sed expressions...>
  local d=runs/rel/$1 np=$2 bin=$3 base=$4; shift 4
  rm -rf "$d"; mkdir -p "$d"
  cp "$base" "$d/gallope.in"
  for e in "$@"; do sed -i -e "$e" "$d/gallope.in"; done
  cases="${cases:+$cases+}$d:$np:$bin"
}

# Namelist items one per line: nvfortran stops reading a group at ';'.
drive=( 's/^  driven = F/  driven = T/' )
undrive=( 's/^  driven    = T/  driven    = F/' )
shear=( 's/^  shear  = .false./  shear  = .true./' 's/^  q      = 0.d0/  q      = 1.5d0/' \
        's/^  dt    = 1.0d-3/  dt    = 5.0d-3/' )
kz0=( 's/^  kmin = 0.d0  0.d0  1.d0/  kmin = 0.d0  0.d0  0.d0/' )
pm2=( 's/^  P_m = 1/  P_m = 2/' )
pm4=( 's/^  P_m = 1/  P_m = 4/' )
nm2=( 's/^  nlz = 32/  nlz = 32\n  nm  = 2/' )
random=( "s/^  init_type = 'OT2'/  init_type = 'random'/" )
io3d=( 's/^  write_intvl_3D      = 0.0d0/  write_intvl_3D      = 0.08d0/' \
       's/^  write_intvl_2D      = 0.1d0/  write_intvl_2D      = 0.0d0/' )
dt2=( 's/^  dt    = 1.0d-3/  dt    = 2.0d-3/' )
gear3=( "s/^  time_step_scheme = 'eSSPIFRK3'/  time_step_scheme = 'gear3'/" )
# Elsasser forcing of phi and psi only (no g)
els=( 's/^  ene_inj_g = 1.0d0/  ene_inj_g = 0.0d0/' 's/^  nfields = 3/  nfields = 2/' \
      "s/^  field_names = .*/  field_names = 'zppe', 'zmpe'/" \
      's/^  amplitudes  = .*/  amplitudes  =  1.0, 1.0/' \
      's/^  frequencies = .*/  frequencies = (0.9, -0.6), (0.9, -0.6)/' )
# RMHD from the KRMHD input: delete what RMHD's namelists lack
krm2rm=( '/^  nm  = /d' '/^  write_hermite_flux /d' \
         '/^  \(mu_hyper_perp\|nexp_perp\|nu_hyper_m\|nexp_m\|beta_i\|tau\|Zcharge\|alpha_root\) /d' )

case $set_name in
r12)
  for v in v10 v11; do
    s=""; [ $v = v10 ] && s=_v1.0
    for np in 1 4; do
      mk r1_mi_${v}_np$np $np MHD_INCOMP$s $R/MHD_INCOMP.in
      mk r1_ms_${v}_np$np $np MHD_INCOMP$s $R/MHD_INCOMP.in "${shear[@]}"
      mk r1_rm_${v}_np$np $np RMHD$s       $R/RMHD.in
      mk r2_mi_${v}_np$np $np MHD_INCOMP$s $R/MHD_INCOMP.in "${drive[@]}"
      mk r2_rm_${v}_np$np $np RMHD$s       $R/RMHD.in       "${drive[@]}"
    done
  done
  mk r2_rk0_v11_np1 1 RMHD      $R/RMHD.in "${drive[@]}" "${kz0[@]}"
  mk r2_rk0_v11_np4 4 RMHD      $R/RMHD.in "${drive[@]}" "${kz0[@]}"
  mk r2_rk0_v10_np1 1 RMHD_v1.0 $R/RMHD.in "${drive[@]}" "${kz0[@]}"
  ;;
r3)
  K=$R/KRMHD.in
  mk r3_kf_np1 1 KRMHD $K
  mk r3_kf_np4 4 KRMHD $K
  mk r3_kf_np2 2 KRMHD $K
  mk r3_kf_pm2 4 KRMHD $K "${pm2[@]}"
  mk r3_kf_pm4 4 KRMHD $K "${pm4[@]}"
  mk r3_kr_rs1 1 KRMHD $K "${undrive[@]}" "${random[@]}" 's/^  nstep = 200/  nstep = 100/'
  mk r3_kr_rs2 1 KRMHD $K "${undrive[@]}" "${random[@]}" 's/^  nstep = 200/  nstep = 100/' \
     "s/^  init_type = 'random'/  init_type = 'restart'/" \
     "s|^  restart_dir = './restart/'|  restart_dir = '../r3_kr_rs1/restart/'|"
  mk r3_kr_full 1 KRMHD $K "${undrive[@]}" "${random[@]}"
  mk r3_kr_np2  2 KRMHD $K "${undrive[@]}" "${random[@]}"
  mk r3_kr_pm2  4 KRMHD $K "${undrive[@]}" "${random[@]}" "${pm2[@]}"
  mk r3p_k_ot2 1 KRMHD $K "${undrive[@]}"
  mk r3p_r_ot2 1 RMHD  $K "${undrive[@]}" "${krm2rm[@]}" "${els[@]}"
  mk r3p_k_els 1 KRMHD $K "${els[@]}"
  mk r3p_r_els 1 RMHD  $K "${els[@]}" "${krm2rm[@]}"
  ;;
r4)
  for np in 1 4; do
    mk r4_mi_np$np $np MHD_INCOMP $R/MHD_INCOMP.in         "${io3d[@]}"
    mk r4_rm_np$np $np RMHD       $R/RMHD.in               "${io3d[@]}"
    mk r4_kf_np$np $np KRMHD      $R/KRMHD.in              "${io3d[@]}"
    mk r4_hr_np$np $np HRMHD      tests/HRMHD/gallope.in   "${io3d[@]}"
    mk r4_rr_np$np $np RRMHD      tests/RRMHD/gallope.in   "${io3d[@]}"
    mk r4_em_np$np $np EMHD       tests/EMHD/gallope.in    "${io3d[@]}"
    mk r4_mi_v10_np$np $np MHD_INCOMP_v1.0_io3d $R/MHD_INCOMP.in "${io3d[@]}"
    mk r4_rm_v10_np$np $np RMHD_v1.0_io3d       $R/RMHD.in       "${io3d[@]}"
  done
  # P_fft = 2: P_m = 1 (np2) vs P_m = 2 (pm2); nm = 2 where the input has no nm
  mk r4_kf_np2 2 KRMHD $R/KRMHD.in "${io3d[@]}"
  mk r4_kf_pm2 4 KRMHD $R/KRMHD.in "${io3d[@]}" "${pm2[@]}"
  for m in rm:RMHD:$R/RMHD.in hr:HRMHD:tests/HRMHD/gallope.in rr:RRMHD:tests/RRMHD/gallope.in; do
    IFS=: read -r c bin base <<< "$m"
    mk r4_${c}_np2 2 $bin $base "${io3d[@]}" "${nm2[@]}"
    mk r4_${c}_pm2 4 $bin $base "${io3d[@]}" "${nm2[@]}" "${pm2[@]}"
  done
  ;;
r5)
  for m in mi:MHD_INCOMP:$R/MHD_INCOMP.in rm:RMHD:$R/RMHD.in kr:KRMHD:$R/KRMHD.in \
           hr:HRMHD:tests/HRMHD/gallope.in em:EMHD:tests/EMHD/gallope.in \
           rr:RRMHD:tests/RRMHD/gallope.in; do
    IFS=: read -r c bin base <<< "$m"
    x=(); [ $c = kr ] && x=( "${undrive[@]}" "${random[@]}" )
    mk r5_${c}_B    1 $bin $base "${x[@]}"
    mk r5_${c}_A    1 $bin $base "${x[@]}" "${dt2[@]}"
    mk r5_${c}_B_g3 1 $bin $base "${x[@]}" "${gear3[@]}"
    mk r5_${c}_A_g3 1 $bin $base "${x[@]}" "${gear3[@]}" "${dt2[@]}"
  done
  mk r5_mi_A_v10 1 MHD_INCOMP_v1.0 $R/MHD_INCOMP.in "${dt2[@]}"
  mk r5_rm_A_v10 1 RMHD_v1.0       $R/RMHD.in       "${dt2[@]}"
  mk r5_kr_A_f05 1 KRMHD_f05d2c1   $R/KRMHD.in      "${undrive[@]}" "${random[@]}" "${dt2[@]}"
  for d in runs/rel/r5_*_A*; do
    echo 1.0d-3 > $d/dt_adjust
  done
  ;;
*) echo "unknown set: $set_name" >&2; exit 1 ;;
esac

echo "$cases"

#!/bin/bash
#=======================================================================
# Create the EMHD validation run directories under runs/ from
# tests/EMHD/gallope.in, then print the CASES string for
# tests/run_cases-tsubame.sh. Run from the repository root.
# Do not run it while a job uses runs/e_*: it recreates the directories.
#
#   ./tests/EMHD/setup_runs.sh
#   CASES=$(./tests/EMHD/setup_runs.sh) RESOURCE=node_f H_RT=0:30:00 \
#     NAME=emhd ./tests/run_cases-tsubame.sh
#
#   e_ot3_np1/np4      V2 smoke, V3 decomposition (1 vs 4 GPUs), V8 vs Calliope
#   e_rs1 -> e_rs2     V5 restart: 100 + 100 steps vs e_ot3_np1
#   e_ideal_dt1/2      V6 ideal W at dt and dt/2 (b0 = z); e_ideal_dt2 is the
#                      dissipation-free twin of e_ot3_np1 for its energy budget
#   e_rnd_ideal_dt1/2  V6 ideal W and H* at dt and dt/2 (random IC, b0 = 0,
#                      where H* is invariant)
#   e_frc_np1/np4      V9 forcing of bx, by, bz with fix_power, 1 vs 4 GPUs
#   e_pm2              P_m = 2 must abort (EMHD supports P_m = P_s = 1 only)
#
#   ./tests/EMHD/setup_runs.sh linear
#   e_lin_001_de0      V7 helical mode (kx,ky,kz)=(0,0,1), s=+1, de=0,   1 GPU
#   e_lin_001_de01     V7 helical mode (0,0,1),            s=-1, de=0.1, 1 GPU
#   e_lin_3m42_de01    V7 helical mode (3,-4,2), s=+1, de=0.1, 4 GPUs (ky on rank 3)
#   e_lin_3m42_de1     V7 helical mode (3,-4,2), s=-1, de=1,   4 GPUs
#=======================================================================
set -e
base=tests/EMHD/gallope.in
set_name=${1:-main}
cases=""

mk() {  # mk <dir> <nranks> <sed expressions...>
  local d=runs/$1 np=$2; shift 2
  rm -rf "$d"; mkdir -p "$d"
  cp "$base" "$d/gallope.in"
  for e in "$@"; do sed -i -e "$e" "$d/gallope.in"; done
  cases="${cases:+$cases+}$d:$np:${model:-EMHD}"
}

# Namelist items one per line: nvfortran stops reading a group at ';'.
ideal=( 's/^  eta     = .*/  eta     = 0.d0/' )

if [ "$set_name" = linear ]; then
  # A helical mode b_k = mode_amp e_s, k x e_s = i s |k| e_s, on the guide field
  # b0 = z is an exact solution of the implemented equations (b1 x u1 = 0):
  #   b_k(t) = b_k(0) exp(i s (k.b0)|k| t/(1 + de^2 k^2))
  #                   exp(-eta (k^2/k2_max)^eta_exp t/(1 + de^2 k^2)).
  # nonlinear = T, since b0 enters through b* x u.
  lin=( "s/^  init_type = 'OT3'/  init_type = 'single_mode'/" \
        's/write_intvl_2D      = 0.1d0/write_intvl_2D      = 0.0d0/' )
  m3m42=( 's/mode_kx  = 0/mode_kx  = 3/' 's/mode_ky  = 0/mode_ky  = -4/' 's/mode_kz  = 1/mode_kz  = 2/' \
          's/nstep = 200/nstep = 500/' )
  mk e_lin_001_de0  1 "${lin[@]}" 's/nstep = 200/nstep = 1000/' 's/^  de      = .*/  de      = 0.d0/'
  mk e_lin_001_de01 1 "${lin[@]}" 's/nstep = 200/nstep = 1000/' 's/mode_hel = 1/mode_hel = -1/'
  mk e_lin_3m42_de01 4 "${lin[@]}" "${m3m42[@]}"
  mk e_lin_3m42_de1  4 "${lin[@]}" "${m3m42[@]}" 's/mode_hel = 1/mode_hel = -1/' \
                       's/^  de      = .*/  de      = 1.d0/'
  echo "$cases"
  exit 0
fi

rnd=( "s/^  init_type = 'OT3'/  init_type = 'random'/" 's/^  b0 = .*/  b0 = 0.d0  0.d0  0.d0/' )

mk e_ot3_np1 1
mk e_ot3_np4 4
mk e_rs1 1 's/nstep = 200/nstep = 100/'
mk e_rs2 1 's/nstep = 200/nstep = 100/' "s/^  init_type = 'OT3'/  init_type = 'restart'/" \
           "s|restart_dir = './restart/'|restart_dir = '../e_rs1/restart/'|"
mk e_ideal_dt1 1 "${ideal[@]}" 's/dt    = 1.0d-3/dt    = 2.0d-3/' 's/nstep = 200/nstep = 100/'
mk e_ideal_dt2 1 "${ideal[@]}"
mk e_rnd_ideal_dt1 1 "${ideal[@]}" "${rnd[@]}" 's/dt    = 1.0d-3/dt    = 2.0d-3/' 's/nstep = 200/nstep = 100/'
mk e_rnd_ideal_dt2 1 "${ideal[@]}" "${rnd[@]}"
mk e_frc_np1 1 's/driven = F/driven = T/'
mk e_frc_np4 4 's/driven = F/driven = T/'
mk e_pm2 2 's/P_m = 1/P_m = 2/'

echo "$cases"

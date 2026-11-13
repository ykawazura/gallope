#!/bin/bash
#=======================================================================
# Create the RRMHD validation run directories under runs/ from
# tests/RRMHD/gallope.in, then print the CASES string for
# tests/run_cases-tsubame.sh. Run from the repository root (GNU sed).
#
#   ./tests/RRMHD/setup_runs.sh
#   CASES=$(./tests/RRMHD/setup_runs.sh) RESOURCE=node_f H_RT=0:30:00 \
#     NAME=rrmhd ./tests/run_cases-tsubame.sh
#
# The base case has q = 1.5, kappa_b = 0.2, kappa_p = 0.3 (beta = 1,
# gamma = 5/3), hyperdissipation, OT3.
#
#   rr_ot3_np1/np4   V2 smoke, V3 decomposition (1 vs 4 GPUs), budget
#   rr_ot3_np2       V4 reference (P_fft=2, P_m=1, nm=2)
#   rr_ot3_pm2       V4 redundancy (P_fft=2, P_m=2, nm=2)
#   rr_rs1 -> rr_rs2 V5 restart: 100 + 100 steps vs rr_ot3_np1
#   rr_ideal_dt1/2   V6 energy conservation (q = kappa = 0, no
#                    dissipation) at dt and dt/2
#   rr_gear3         budget with gear3 (kappa terms in the gear3 path)
#   rr_driven        driven = T must stop the run
#
#   ./tests/RRMHD/setup_runs.sh linear     (nonlinear = F, no dissipation)
#   rr_lin_ap, rr_lin_am  V7 (kx,ky,kz) = (4,0,1), q = kappa = 0:
#                    Alfven phi = +-psi (omega = +-kz) and slow
#                    upa = +-sqrt(A) bpa (omega = +-kz/sqrt(A)), A = 2.2
#   rr_lin_rot_a     V7 (3,-4,2), base q, kappa: eigenmode, 4 GPUs
#   rr_lin_rot_b     V7 (3,-4,2), base q, kappa: another eigenmode, 1 GPU
#   rr_lin_unst      V7 (3,-4,1), base q, kappa: generic initial vector
#                    containing an unstable pair (K kz^2 < 2 q ky^2)
#=======================================================================
set -e
base=tests/RRMHD/gallope.in
set_name=${1:-main}
cases=""

mk() {  # mk <dir> <nranks> <sed expressions...>
  local d=runs/$1 np=$2; shift 2
  if [ -d "$d" ]; then rm -r "$d"; fi
  mkdir -p "$d"
  cp "$base" "$d/gallope.in"
  for e in "$@"; do sed -i -e "$e" "$d/gallope.in"; done
  cases="${cases:+$cases+}$d:$np:${model:-RRMHD}"
}

# one namelist item per line: nvfortran stops reading a group at the first ';'
ideal=( 's/^  \(nu\|eta\)\(pe\|pa\)_\([xz]\) *= *.*/  \1\2_\3 = 0.d0/' )
norot=( 's/^  q       = 1.5d0/  q       = 0.0d0/' \
        's/^  kappa_b = 0.2d0/  kappa_b = 0.0d0/' \
        's/^  kappa_p = 0.3d0/  kappa_p = 0.0d0/' )

if [ "$set_name" = linear ]; then
  # Gallope advances d/dt (phi, psi, upa, bpa) = i R (phi, psi, upa, bpa) with
  # the real matrix R of advance::get_ext_terms; an eigenvector of R with real
  # eigenvalue w is real and every component rotates as exp(i w t).
  # tests/RRMHD/check_runs.py recomputes R and compares the final state with
  # expm(i R t) applied to the amplitudes below.
  lin=( "${ideal[@]}" 's/nonlinear = T/nonlinear = F/' \
        "s/^  init_type = 'OT3'/  init_type = 'single_mode'/" \
        's/write_intvl_2D      = 0.1d0/write_intvl_2D      = 0.0d0/' \
        's/nstep = 200/nstep = 1000/' )
  amps() { a=( "s/amp_phi = .*/amp_phi = $1/" "s/amp_psi = .*/amp_psi = $2/" \
               "s/amp_upa = .*/amp_upa = $3/" "s/amp_bpa = .*/amp_bpa = $4/" ); }
  m401=( 's/mode_kx = 1/mode_kx = 4/' )
  m3m42=( 's/mode_kx = 1/mode_kx = 3/' 's/mode_ky = 0/mode_ky = -4/' 's/mode_kz = 1/mode_kz = 2/' )
  m3m41=( 's/mode_kx = 1/mode_kx = 3/' 's/mode_ky = 0/mode_ky = -4/' )

  # sqrt(A) = sqrt(1 + 2/(beta gamma)) = sqrt(2.2)
  amps 1.d0  1.d0  1.483239697419133d0 1.d0
  mk rr_lin_ap 1 "${lin[@]}" "${norot[@]}" "${m401[@]}" "${a[@]}"
  amps 1.d0 -1.d0 -1.483239697419133d0 1.d0
  mk rr_lin_am 1 "${lin[@]}" "${norot[@]}" "${m401[@]}" "${a[@]}"

  # eigenvectors of R at (3,-4,2), q = 1.5, kappa_b = 0.2, kappa_p = 0.3:
  #   w = -2.614574714393141 and w = 1.009233791672201
  amps 1.d0 -7.649427606677560d-01 2.525856643936327d+00 -1.919549076375046d+00
  mk rr_lin_rot_a 4 "${lin[@]}" "${m3m42[@]}" "${a[@]}"
  amps 1.d0 1.981701382279518d+00 1.194931763733975d+01 6.038976882374196d+00
  mk rr_lin_rot_b 1 "${lin[@]}" "${m3m42[@]}" "${a[@]}"

  # (3,-4,1): w = 0.2148 +- 0.4144 i among the roots
  amps 1.d0 5.d-1 -3.d-1 2.d-1
  mk rr_lin_unst 1 "${lin[@]}" "${m3m41[@]}" "${a[@]}"
  echo "$cases"
  exit 0
fi

mk rr_ot3_np1 1
mk rr_ot3_np4 4
mk rr_ot3_np2 2 's/^  nlz = 32/  nlz = 32\n  nm  = 2/'
mk rr_ot3_pm2 4 's/^  nlz = 32/  nlz = 32\n  nm  = 2/' 's/P_m = 1/P_m = 2/'
mk rr_rs1 1 's/nstep = 200/nstep = 100/'
mk rr_rs2 1 's/nstep = 200/nstep = 100/' "s/^  init_type = 'OT3'/  init_type = 'restart'/" \
            "s|restart_dir = './restart/'|restart_dir = '../rr_rs1/restart/'|"
mk rr_ideal_dt1 1 "${ideal[@]}" "${norot[@]}" 's/dt    = 1.0d-3/dt    = 2.0d-3/' 's/nstep = 200/nstep = 100/'
mk rr_ideal_dt2 1 "${ideal[@]}" "${norot[@]}"
mk rr_gear3 1 "s/time_step_scheme = 'eSSPIFRK3'/time_step_scheme = 'gear3'/"
mk rr_driven 1 's/driven = F/driven = T/'

echo "$cases"

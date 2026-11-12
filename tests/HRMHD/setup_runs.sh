#!/bin/bash
#=======================================================================
# Create the HRMHD validation run directories under runs/ from
# tests/HRMHD/gallope.in, then print the CASES string for
# tests/run_cases-tsubame.sh. Run from the repository root.
#
#   ./tests/HRMHD/setup_runs.sh
#   CASES=$(./tests/HRMHD/setup_runs.sh) RESOURCE=node_f H_RT=0:30:00 \
#     NAME=hrmhd ./tests/run_cases-tsubame.sh
#
#   h_ot3_np1/np4   V2 smoke, V3 decomposition (1 vs 4 GPUs)
#   h_ot3_np2       V4 reference (P_fft=2, P_m=1)
#   h_ot3_pm2       V4 redundancy (P_fft=2, P_m=2)
#   h_rs1 -> h_rs2  V5 restart: 100 + 100 steps vs h_ot3_np1
#   h_ideal_dt1/2   V6 ideal invariants W, X at dt and dt/2
#   h_frc_np1/np4   V9 forcing (phi, psi, upa, bpa), 1 vs 4 GPUs
#   h_rho0          rho = 0 (KAW/ICW split off, RMHD + passive upa/bpa)
#   h_rho0_ideal    h_rho0 without dissipation: integrator drift for its budget
#   h_frc_els_np1/np4  V9 Elsasser forcing (zppe, zmpe), xhl_inj = 0.5
#
#   ./tests/HRMHD/setup_runs.sh linear
#   h_lin_{kaw,icw}_401   V7 eigenmode (kx,ky,kz)=(4,0,1), 1 GPU
#   h_lin_{kaw,icw}_3m42  V7 eigenmode (3,-4,2), 4 GPUs (ky on rank 3)
#   r_ot3                 rho -> 0 limit: RMHD with the input of h_rho0
#=======================================================================
set -e
base=tests/HRMHD/gallope.in
set_name=${1:-main}
cases=""

mk() {  # mk <dir> <nranks> <sed expressions...>
  local d=runs/$1 np=$2; shift 2
  rm -rf "$d"; mkdir -p "$d"
  cp "$base" "$d/gallope.in"
  for e in "$@"; do sed -i -e "$e" "$d/gallope.in"; done
  cases="${cases:+$cases+}$d:$np:${model:-HRMHD}"
}

# Namelist items one per line: nvfortran stops reading a group at ';'.
ideal=( 's/^  \(nu\|eta\)\(pe\|pa\)_\([xz]\) *= .*/  \1\2_\3 = 0.d0/' )

if [ "$set_name" = linear ]; then
  # Eigenvectors of the implemented linear operator (rho = 0.5, sgm = 1):
  #   lambda^4 - (1 + sgm^2 + rho^2 kprp^2) lambda^2 + sgm^2 = 0,
  #   (phi, psi, upa, bpa) = (1, lambda, sgm(lambda^2-1)/(rho lambda), (lambda^2-1)/rho),
  #   each Fourier coefficient rotating as exp(i kz lambda t).
  lin=( "${ideal[@]}" 's/nonlinear = T/nonlinear = F/' \
        "s/^  init_type = 'OT3'/  init_type = 'single_mode'/" \
        's/write_intvl_2D      = 0.1d0/write_intvl_2D      = 0.0d0/' )
  amps() { a=( "s/amp_phi = .*/amp_phi = $1/" "s/amp_psi = .*/amp_psi = $2/" \
               "s/amp_upa = .*/amp_upa = $3/" "s/amp_bpa = .*/amp_bpa = $4/" ); }
  m401=( 's/mode_kx = 1/mode_kx = 4/' 's/nstep = 200/nstep = 1000/' )
  m3m42=( 's/mode_kx = 1/mode_kx = 3/' 's/mode_ky = 0/mode_ky = -4/' 's/mode_kz = 1/mode_kz = 2/' \
          's/nstep = 200/nstep = 500/' )
  amps 1.d0 2.414213562373095d0 4.000000000000000d0 9.656854249492380d0
  mk h_lin_kaw_401 1 "${lin[@]}" "${m401[@]}" "${a[@]}"
  amps 1.d0 4.142135623730948d-1 -4.000000000000004d0 -1.656854249492381d0
  mk h_lin_icw_401 1 "${lin[@]}" "${m401[@]}" "${a[@]}"
  amps 1.d0 2.850781059358212d0 5.000000000000000d0 1.425390529679106d1
  mk h_lin_kaw_3m42 4 "${lin[@]}" "${m3m42[@]}" "${a[@]}"
  amps 1.d0 3.507810593582123d-1 -4.999999999999997d0 -1.753905296791061d0
  mk h_lin_icw_3m42 4 "${lin[@]}" "${m3m42[@]}" "${a[@]}"

  # RMHD with the input of h_rho0, less what its namelists lack
  # (rho, sgm, nupa_*, etapa_*, mode_*, amp_*, forcing of upa/bpa)
  model=RMHD mk r_ot3 1 '/^  rho = /d' '/^  sgm = /d' '/^  \(nu\|eta\)pa_[xz]/d' \
                        '/^  \(mode\|amp\)_/d' 's/nfields = 4/nfields = 2/' \
                        "s/^  field_names = .*/  field_names = 'phi', 'psi'/" \
                        's/^  amplitudes  = .*/  amplitudes  =  1.0, 1.0/' \
                        's/^  frequencies = .*/  frequencies = (0.9, -0.6), (0.9, -0.6)/'
  echo "$cases"
  exit 0
fi

mk h_ot3_np1 1
mk h_ot3_np4 4
mk h_ot3_np2 2 's/^  nlz = 32/  nlz = 32\n  nm  = 2/'
mk h_ot3_pm2 4 's/^  nlz = 32/  nlz = 32\n  nm  = 2/' 's/P_m = 1/P_m = 2/'
mk h_rs1 1 's/nstep = 200/nstep = 100/'
mk h_rs2 1 's/nstep = 200/nstep = 100/' "s/^  init_type = 'OT3'/  init_type = 'restart'/" \
           "s|restart_dir = './restart/'|restart_dir = '../h_rs1/restart/'|"
mk h_ideal_dt1 1 "${ideal[@]}" 's/dt    = 1.0d-3/dt    = 2.0d-3/' 's/nstep = 200/nstep = 100/'
mk h_ideal_dt2 1 "${ideal[@]}"
mk h_frc_np1 1 's/driven = F/driven = T/'
mk h_frc_np4 4 's/driven = F/driven = T/'
mk h_rho0 1 's/rho = 0.5d0/rho = 0.0d0/'
mk h_rho0_ideal 1 's/rho = 0.5d0/rho = 0.0d0/' "${ideal[@]}"
els=( 's/driven = F/driven = T/' 's/^  elsasser = .false./  elsasser = .true./' \
      's/^  ene_inj = 1.0d0/  ene_inj = 1.0d0\n  xhl_inj = 0.5d0/' 's/nfields = 4/nfields = 2/' \
      "s/^  field_names = .*/  field_names = 'zppe', 'zmpe'/" \
      's/^  amplitudes  = .*/  amplitudes  =  1.0, 1.0/' \
      's/^  frequencies = .*/  frequencies = (0.9, -0.6), (0.9, -0.6)/' )
mk h_frc_els_np1 1 "${els[@]}"
mk h_frc_els_np4 4 "${els[@]}"

echo "$cases"

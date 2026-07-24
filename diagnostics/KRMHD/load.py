# -*- coding: utf-8 -*-
# KRMHD diagnostics loader.
# Reads <runname>.out.nc and rms.dat, exposing the Alfvenic (EM) diagnostics
# inherited from RMHD plus the compressive Hermite-moment (g) diagnostics:
#   W_free(t)          -- Meyrand et al. 2019 eq.6 free energy (alpha-weighted)
#   W_m(t, m)          -- per-Hermite-moment free energy (k-integrated, exact)
#   g2_bin(t,kz,m,kp)  -- binned 1/2 |g_m|^2 spectrum
#   Gamma_m(t,kz,m,kp) -- Hermite flux (eq.9), only if write_hermite_flux was on
#   Gamma_m_kint(t,m)  -- k-integrated Hermite flux, only if present
import warnings
warnings.filterwarnings('ignore')

import numpy as np
from scipy.io import netcdf

####### time index for final cut #########
input_dir = '../../'
runname = 'gallope'
restart_num = ''
final_idx     = -1
final_fld_idx = -1

if final_idx != -1:
  print ('\n!!! CAUTION: final_idx = %d !!!\n' % final_idx)

####### movie or final cut #########
ismovie = False

####### ignore these points #########
ignored_points     = [0]
ignored_points_fld = [0]


####### load netcdf file #########
ncfile = netcdf.netcdf_file(input_dir+runname+'.out.nc'+restart_num, 'r')

# Load parameters

# Load coordinate
tt  = np.copy(ncfile.variables['tt' ][:]); tt  = np.delete(tt , ignored_points, axis = 0)
xx  = np.copy(ncfile.variables['xx' ][:])
yy  = np.copy(ncfile.variables['yy' ][:])
zz  = np.copy(ncfile.variables['zz' ][:])
kx  = np.copy(ncfile.variables['kx' ][:])
ky  = np.copy(ncfile.variables['ky' ][:])
kz  = np.copy(ncfile.variables['kz' ][:])
kpbin  = np.copy(ncfile.variables['kpbin'][:])
mm  = np.copy(ncfile.variables['mm' ][:])   # Hermite index m = 0..nm-1

nt  = tt.size
nlx = xx.size
nly = yy.size
nlz = zz.size
nkx = kx.size
nky = ky.size
nkz = kz.size
nkpolar = kpbin.size
nm  = mm.size

if nkz <= 2:
  is2D = True
else:
  is2D = False

# Load total energies (Alfvenic / EM)
upe2_sum        = np.copy(ncfile.variables['upe2_sum'       ][:]); upe2_sum        = np.delete(upe2_sum       , ignored_points, axis = 0)
bpe2_sum        = np.copy(ncfile.variables['bpe2_sum'       ][:]); bpe2_sum        = np.delete(bpe2_sum       , ignored_points, axis = 0)
zppe2_sum       = np.copy(ncfile.variables['zppe2_sum'      ][:]); zppe2_sum       = np.delete(zppe2_sum      , ignored_points, axis = 0)
zmpe2_sum       = np.copy(ncfile.variables['zmpe2_sum'      ][:]); zmpe2_sum       = np.delete(zmpe2_sum      , ignored_points, axis = 0)
upe2dot_sum     = np.copy(ncfile.variables['upe2dot_sum'    ][:]); upe2dot_sum     = np.delete(upe2dot_sum    , ignored_points, axis = 0)
bpe2dot_sum     = np.copy(ncfile.variables['bpe2dot_sum'    ][:]); bpe2dot_sum     = np.delete(bpe2dot_sum    , ignored_points, axis = 0)
upe2dissip_sum  = np.copy(ncfile.variables['upe2dissip_sum' ][:]); upe2dissip_sum  = np.delete(upe2dissip_sum , ignored_points, axis = 0)
bpe2dissip_sum  = np.copy(ncfile.variables['bpe2dissip_sum' ][:]); bpe2dissip_sum  = np.delete(bpe2dissip_sum , ignored_points, axis = 0)
p_phi_sum       = np.copy(ncfile.variables['p_phi_sum'      ][:]); p_phi_sum       = np.delete(p_phi_sum      , ignored_points, axis = 0)
p_psi_sum       = np.copy(ncfile.variables['p_psi_sum'      ][:]); p_psi_sum       = np.delete(p_psi_sum      , ignored_points, axis = 0)
p_xhl_sum       = np.copy(ncfile.variables['p_xhl_sum'      ][:]); p_xhl_sum       = np.delete(p_xhl_sum      , ignored_points, axis = 0)

# Load compressive free energy (g). W_free already includes the alpha-weighted
# m=0 rung (eq.6) and the factor-of-2 negative-kz compensation (see io.F90).
W_free          = np.copy(ncfile.variables['W_free'         ][:]); W_free          = np.delete(W_free         , ignored_points, axis = 0)
W_m             = np.copy(ncfile.variables['W_m'            ][:]); W_m             = np.delete(W_m            , ignored_points, axis = 0)
# g free-energy power balance: injection P_g and dissipation D_g (dW_free/dt = P_g - D_g)
p_g_sum         = np.copy(ncfile.variables['p_g_sum'        ][:]); p_g_sum         = np.delete(p_g_sum        , ignored_points, axis = 0)
Dg_sum          = np.copy(ncfile.variables['Dg_sum'         ][:]); Dg_sum          = np.delete(Dg_sum         , ignored_points, axis = 0)

# Load binned spectra (Alfvenic / EM)
upe2_bin        = np.copy(ncfile.variables['upe2_bin'       ][:]); upe2_bin        = np.delete(upe2_bin       , ignored_points, axis = 0)
bpe2_bin        = np.copy(ncfile.variables['bpe2_bin'       ][:]); bpe2_bin        = np.delete(bpe2_bin       , ignored_points, axis = 0)
zppe2_bin       = np.copy(ncfile.variables['zppe2_bin'      ][:]); zppe2_bin       = np.delete(zppe2_bin      , ignored_points, axis = 0)
zmpe2_bin       = np.copy(ncfile.variables['zmpe2_bin'      ][:]); zmpe2_bin       = np.delete(zmpe2_bin      , ignored_points, axis = 0)

# Load binned g spectrum: g2_bin(tt, kz, mm, kpbin) = 1/2 |g_m|^2 (binned in kperp)
g2_bin          = np.copy(ncfile.variables['g2_bin'         ][:]); g2_bin          = np.delete(g2_bin         , ignored_points, axis = 0)

# Hermite flux is optional (only written when write_hermite_flux = .true.)
has_hermite_flux = ('Gamma_m_kint' in ncfile.variables)
if has_hermite_flux:
  Gamma_m_kint  = np.copy(ncfile.variables['Gamma_m_kint'   ][:]); Gamma_m_kint    = np.delete(Gamma_m_kint   , ignored_points, axis = 0)
  Gamma_m       = np.copy(ncfile.variables['Gamma_m'        ][:]); Gamma_m         = np.delete(Gamma_m        , ignored_points, axis = 0)
else:
  Gamma_m_kint  = None
  Gamma_m       = None

ncfile.close()

# The final diagnostic record is written twice (normal write_intvl cadence +
# the termination flush), producing a duplicate timestamp (dt=0) that makes
# np.gradient(., tt) blow up. Drop duplicate-time records from every
# time-series array (same treatment as rms.dat below), keeping the first
# occurrence so tt stays strictly increasing.
_, _keep = np.unique(tt, return_index=True)
if _keep.size != tt.size:
  _ntt = tt.size
  for _name in ['tt', 'upe2_sum', 'bpe2_sum', 'zppe2_sum', 'zmpe2_sum',
                'upe2dot_sum', 'bpe2dot_sum', 'upe2dissip_sum', 'bpe2dissip_sum',
                'p_phi_sum', 'p_psi_sum', 'p_xhl_sum', 'W_free', 'W_m',
                'p_g_sum', 'Dg_sum', 'upe2_bin', 'bpe2_bin', 'zppe2_bin',
                'zmpe2_bin', 'g2_bin', 'Gamma_m_kint', 'Gamma_m']:
    _arr = globals().get(_name)
    if _arr is not None and getattr(_arr, 'shape', (0,))[0] == _ntt:
      globals()[_name] = _arr[_keep]
  nt = tt.size

tlab  = r'$(v_\mathrm{A}/L_\|)t$'
zlab  = r'$z/L_\|$'
kzlab = r'$L_\| k_z$'
mlab  = r'$m$'

# Load rms
import os
filename = input_dir+'rms.dat'+restart_num
if os.path.isfile(filename):
  data = np.loadtxt(filename).T

  # Reduce multiple points at the same time because of multistep update (RK or Gear)
  _, indices = np.unique(data[0], return_index=True)
  data = np.asarray([data.T[i] for i in indices]).T
  tt_rms, ux_rms, uy_rms, bx_rms, by_rms = data

  # Pick up where tt_rms = tt
  indices = [np.abs(tt_rms - t).argmin() for t in tt]
  tt_rms, ux_rms, uy_rms, bx_rms, by_rms = np.asarray([data.T[i] for i in indices]).T
else:
  tt_rms, ux_rms, uy_rms, bx_rms, by_rms = [np.array([0])] * 5

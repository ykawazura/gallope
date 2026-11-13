# -*- coding: utf-8 -*-
import warnings
warnings.filterwarnings('ignore')

import numpy as np
from scipy.io import netcdf

####### time index for final cut #########
input_dir = '../../'
runname = 'gallope'
restart_num = ''
final_idx      = -1
final_fld_idx  = -1

if final_idx != -1:
  print ('\n!!! CAUTION: final_idx = %d !!!\n' % final_idx)

####### movie or final cut #########
ismovie = False

####### ignore these points #########
ignored_points      = [0]
ignored_points_fld  = [ ]


####### load netcdf file #########
ncfile = netcdf.netcdf_file(input_dir+runname+'.out.nc'+restart_num, 'r')   

# Load parameters
de = np.copy(ncfile.variables['de'].data)

# Load coordinate
tt  = np.copy(ncfile.variables['tt' ][:]); tt  = np.delete(tt , ignored_points, axis = 0)
xx  = np.copy(ncfile.variables['xx' ][:])
yy  = np.copy(ncfile.variables['yy' ][:])
zz  = np.copy(ncfile.variables['zz' ][:])
kx  = np.copy(ncfile.variables['kx' ][:])
ky  = np.copy(ncfile.variables['ky' ][:])
kz  = np.copy(ncfile.variables['kz' ][:])
kpbin  = np.copy(ncfile.variables['kpbin'][:])

nt  = tt.size
nlx = xx.size
nly = yy.size
nlz = zz.size
nkx = kx.size
nky = ky.size
nkz = kz.size
nkpolar = kpbin.size

if nkz <= 2:
  is2D = True
else:
  is2D = False

# Load total energies
wmag_sum         = np.copy(ncfile.variables['wmag_sum'        ][:]); wmag_sum         = np.delete(wmag_sum        , ignored_points, axis = 0)
wmag_dot_sum     = np.copy(ncfile.variables['wmag_dot_sum'    ][:]); wmag_dot_sum     = np.delete(wmag_dot_sum    , ignored_points, axis = 0)
wmag_dissip_sum  = np.copy(ncfile.variables['wmag_dissip_sum' ][:]); wmag_dissip_sum  = np.delete(wmag_dissip_sum , ignored_points, axis = 0)
p_ext_sum        = np.copy(ncfile.variables['p_ext_sum'       ][:]); p_ext_sum        = np.delete(p_ext_sum       , ignored_points, axis = 0)

# mean magnetic field
b0            = np.copy(ncfile.variables['b0'           ][:]); b0            = np.delete(b0           , ignored_points, axis = 0)

# Load binned spectra
b2_bin        = np.copy(ncfile.variables['b2_bin'       ][:]); b2_bin        = np.delete(b2_bin       , ignored_points, axis = 0)
bx2_bin       = np.copy(ncfile.variables['bx2_bin'      ][:]); bx2_bin       = np.delete(bx2_bin      , ignored_points, axis = 0)
by2_bin       = np.copy(ncfile.variables['by2_bin'      ][:]); by2_bin       = np.delete(by2_bin      , ignored_points, axis = 0)
bz2_bin       = np.copy(ncfile.variables['bz2_bin'      ][:]); bz2_bin       = np.delete(bz2_bin      , ignored_points, axis = 0)

ncfile.close()

tlab  = r'$\Omega_\rmi t$'
xlab  = r'$x/d_\rmi$'
ylab  = r'$y/d_\rmi$'
zlab  = r'$z/d_\rmi$'
kxlab = r'$k_x d_\rmi$'
kylab = r'$k_y d_\rmi$'
kzlab = r'$k_z d_\rmi$'
# the spectral bins are 3D isotropic shells in |k|, not in k_perp
kplab = r'$k d_\rmi$'

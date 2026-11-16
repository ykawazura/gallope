# -*- coding: utf-8 -*-
from load import *
from fft import *
import sys
sys.path.append('../')
from plots import *

print('\nplotting kspectrum\n')
outdir = './fig_kspectrum/'

# The bins are 3D isotropic shells that already contain the -kz modes
# (Hermitian doubling), so no folding of +-kz is needed here.
# Only the 2/3-dealiased range of the shells is fully populated.
if nlz == nkz:
  kp_end = np.argmin(np.abs(kpbin - kpbin.max()*2./3.))
else:
  kp_end = int(kpbin.size*2./3.)

# If shear is on, add a vertical line indicating the fastest-MRI mode
if shear_flg == 1:
  b00 = np.max(b0[0])
  k_mri = np.sqrt(15.)/4./b00

def add_vertical_line(xs, ys, ls, legends):
  xs.append([k_mri, k_mri])
  ys.append([np.min([u2_bin[final_idx, 1:kp_end].min(), b2_bin[final_idx, 1:kp_end].min()]),
             np.max([u2_bin[final_idx, 1:kp_end].max(), b2_bin[final_idx, 1:kp_end].max()])])
  ls.append('k:')
  legends.append(r'$k = k_\mr{MRI} = (\sqrt{15}/4)\varpi_0/v_\rmA$')

  return xs, ys, ls, legends

#--------------------------------------------------------#
#                      plot 1D spectra                   #
#--------------------------------------------------------#
ys = [
       u2_bin[final_idx, 1:kp_end],
       b2_bin[final_idx, 1:kp_end],
       kpbin[1:kp_end]**(-5./3.)/kpbin[1]**(-5./3.)*b2_bin[final_idx,1:kp_end][0],
       kpbin[1:kp_end]**(-3./2.)/kpbin[1]**(-3./2.)*b2_bin[final_idx,1:kp_end][0]
     ]
xs = [
      kpbin[1:kp_end],
      kpbin[1:kp_end],
      kpbin[1:kp_end],
      kpbin[1:kp_end]
     ]
ls = [
        '',
        '',
        'k--',
        'k-.',
     ]
legends = [
            r'$E_{u}$',
            r'$E_{B}$',
            r'-5/3',
            r'-3/2',
          ]

if shear_flg == 1: xs, ys, ls, legends = add_vertical_line(xs, ys, ls, legends)

plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt[final_idx], ylab='', term=True, save=outdir+'k_spectra.pdf')


# u by components
ys = [
       ux2_bin[final_idx, 1:kp_end],
       uy2_bin[final_idx, 1:kp_end],
       uz2_bin[final_idx, 1:kp_end],
       kpbin[1:kp_end]**(-5./3.)/kpbin[1]**(-5./3.)*u2_bin[final_idx,1:kp_end][0],
       kpbin[1:kp_end]**(-3./2.)/kpbin[1]**(-3./2.)*u2_bin[final_idx,1:kp_end][0],
     ]
xs = [
      kpbin[1:kp_end],
      kpbin[1:kp_end],
      kpbin[1:kp_end],
      kpbin[1:kp_end],
      kpbin[1:kp_end],
     ]
ls = [
        '',
        '',
        '',
        'k--',
        'k-.',
     ]
legends = [
            r'$E_{u_x}$',
            r'$E_{u_y}$',
            r'$E_{u_z}$',
            r'-5/3',
            r'-3/2',
          ]

if shear_flg == 1: xs, ys, ls, legends = add_vertical_line(xs, ys, ls, legends)

plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt[final_idx], ylab='', term=True, save=outdir+'k_spectra_u.pdf')


# B by components
ys = [
       bx2_bin[final_idx, 1:kp_end],
       by2_bin[final_idx, 1:kp_end],
       bz2_bin[final_idx, 1:kp_end],
       kpbin[1:kp_end]**(-5./3.)/kpbin[1]**(-5./3.)*b2_bin[final_idx,1:kp_end][0],
       kpbin[1:kp_end]**(-3./2.)/kpbin[1]**(-3./2.)*b2_bin[final_idx,1:kp_end][0],
     ]
xs = [
      kpbin[1:kp_end],
      kpbin[1:kp_end],
      kpbin[1:kp_end],
      kpbin[1:kp_end],
      kpbin[1:kp_end],
     ]
ls = [
        '',
        '',
        '',
        'k--',
        'k-.',
     ]
legends = [
            r'$E_{B_x}$',
            r'$E_{B_y}$',
            r'$E_{B_z}$',
            r'-5/3',
            r'-3/2',
          ]

if shear_flg == 1: xs, ys, ls, legends = add_vertical_line(xs, ys, ls, legends)

plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt[final_idx], ylab='', term=True, save=outdir+'k_spectra_b.pdf')


# Elsasser fields
ys = [
       zp2_bin[final_idx, 1:kp_end],
       zm2_bin[final_idx, 1:kp_end],
       kpbin[1:kp_end]**(-5./3.)/kpbin[1]**(-5./3.)*zp2_bin[final_idx,1:kp_end][0],
       kpbin[1:kp_end]**(-3./2.)/kpbin[1]**(-3./2.)*zp2_bin[final_idx,1:kp_end][0]
     ]
xs = [
      kpbin[1:kp_end],
      kpbin[1:kp_end],
      kpbin[1:kp_end],
      kpbin[1:kp_end]
     ]
ls = [
        '',
        '',
        'k--',
        'k-.',
     ]
legends = [
            r'$E_{Z^+}$',
            r'$E_{Z^-}$',
            r'-5/3',
            r'-3/2',
          ]

if shear_flg == 1: xs, ys, ls, legends = add_vertical_line(xs, ys, ls, legends)

plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt[final_idx], ylab='', term=True, save=outdir+'k_spectra_ELS.pdf')


# Thermodynamic fields
ys = [
       rho_bin[final_idx, 1:kp_end],
       kpbin[1:kp_end]**(-5./3.)/kpbin[1]**(-5./3.)*rho_bin[final_idx,1:kp_end][0],
       kpbin[1:kp_end]**(-3./2.)/kpbin[1]**(-3./2.)*rho_bin[final_idx,1:kp_end][0]
     ]
xs = [
      kpbin[1:kp_end],
      kpbin[1:kp_end],
      kpbin[1:kp_end]
     ]
ls = [
        '',
        'k--',
        'k-.',
     ]
legends = [
            r'$E_{\rho}$',
            r'-5/3',
            r'-3/2',
          ]

if shear_flg == 1: xs, ys, ls, legends = add_vertical_line(xs, ys, ls, legends)

plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt[final_idx], ylab='', term=True, save=outdir+'k_spectra_therm.pdf')

#------------------#
#   output ascii   #
#------------------#
np.savetxt(outdir + 'Ek.txt'  , np.column_stack((kpbin[:kp_end],
                                                  u2_bin[final_idx,:kp_end],
                                                 ux2_bin[final_idx,:kp_end],
                                                 uy2_bin[final_idx,:kp_end],
                                                 uz2_bin[final_idx,:kp_end],
                                                  b2_bin[final_idx,:kp_end],
                                                 bx2_bin[final_idx,:kp_end],
                                                 by2_bin[final_idx,:kp_end],
                                                 bz2_bin[final_idx,:kp_end],
                                                 rho_bin[final_idx,:kp_end],
                                                 )), fmt='%E')

#--------------------------------------------------------#
#                       plot movie                       #
#--------------------------------------------------------#
if ismovie:
  ys = []
  ys.append(u2_bin[:, 1:kp_end])
  ys.append(b2_bin[:, 1:kp_end])
  ys.append(np.tile(kpbin[1:kp_end], (nt, 1))**(-5./3.)/kpbin[1]*np.tile(u2_bin[:,1], (kp_end-1, 1)).T)
  ys.append(np.tile(kpbin[1:kp_end], (nt, 1))**(-3./2.)/kpbin[1]*np.tile(u2_bin[:,1], (kp_end-1, 1)).T)
  ys = np.transpose(ys, (1, 0, 2))

  xs = [
        kpbin[1:kp_end],
        kpbin[1:kp_end],
        kpbin[1:kp_end],
        kpbin[1:kp_end]
       ]
  ls = [
          '',
          '',
          'k--',
          'k-.',
       ]
  legends = [
              r'$E_{u}$',
              r'$E_{b}$',
              r'-5/3',
              r'-3/2',
            ]
  movie_log1d_many(tt, xs, ys, legends=legends, ls=ls, legendloc='lower left', xlab=kplab, save=outdir+'k_spectra.gif')


  # Elsasser fields
  ys = []
  ys.append(zp2_bin[:, 1:kp_end])
  ys.append(zm2_bin[:, 1:kp_end])
  ys.append(np.tile(kpbin[1:kp_end], (nt, 1))**(-5./3.)/kpbin[1]*np.tile(zp2_bin[:,1], (kp_end-1, 1)).T)
  ys.append(np.tile(kpbin[1:kp_end], (nt, 1))**(-3./2.)/kpbin[1]*np.tile(zp2_bin[:,1], (kp_end-1, 1)).T)
  ys = np.transpose(ys, (1, 0, 2))

  xs = [
        kpbin[1:kp_end],
        kpbin[1:kp_end],
        kpbin[1:kp_end],
        kpbin[1:kp_end]
       ]
  ls = [
          '',
          '',
          'k--',
          'k-.',
       ]
  legends = [
              r'$E_{Z^+}$',
              r'$E_{Z^-}$',
              r'-5/3',
              r'-3/2',
            ]
  movie_log1d_many(tt, xs, ys, legends=legends, ls=ls, legendloc='lower left', xlab=kplab, save=outdir+'k_spectra_ELS.gif')


  # Thermodynamic fields
  ys = []
  ys.append(rho_bin[:, 1:kp_end])
  ys.append(np.tile(kpbin[1:kp_end], (nt, 1))**(-5./3.)/kpbin[1]*np.tile(rho_bin[:,1], (kp_end-1, 1)).T)
  ys.append(np.tile(kpbin[1:kp_end], (nt, 1))**(-3./2.)/kpbin[1]*np.tile(rho_bin[:,1], (kp_end-1, 1)).T)
  ys = np.transpose(ys, (1, 0, 2))

  xs = [
        kpbin[1:kp_end],
        kpbin[1:kp_end],
        kpbin[1:kp_end]
       ]
  ls = [
          '',
          'k--',
          'k-.',
       ]
  legends = [
              r'$E_{\rho}$',
              r'-5/3',
              r'-3/2',
            ]
  movie_log1d_many(tt, xs, ys, legends=legends, ls=ls, legendloc='lower left', xlab=kplab, save=outdir+'k_spectra_therm.gif')

#--------------------------------------------------------#
#                       2D spectra                       #
#--------------------------------------------------------#
# out2d/time.dat has two columns: tt, tsc
tt_fld = np.transpose(np.loadtxt(input_dir+'out2d'+restart_num+'/time.dat'))[0]
nt_fld = tt_fld.size
if nt_fld == 1 : tt_fld = [tt_fld]
tt_fld = np.delete(tt_fld, ignored_points_fld, axis = 0)

# kx and ky are full and stored in FFT order, so they are shifted to ascending order.
# kz is the r2c half (0 <= kz <= kz_max) and is not shifted; the kz > 0 planes
# already carry the factor 2 for the -kz modes (diagnostics.F90).
nkx_mid = np.where(np.sign(kx) == -1)[0][0]
nky_mid = np.where(np.sign(ky) == -1)[0][0]
kx_ = np.concatenate([kx[nkx_mid:], kx[0:nkx_mid]])
ky_ = np.concatenate([ky[nky_mid:], ky[0:nky_mid]])

# record shapes (Fortran order): kxy (nkx, nky), kyz (nky, nkz), kxz (nkx, nkz)
u2_kxy  = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/u2_kxy_sum_kz.dat').reshape(nt_fld, nky, nkx)); u2_kxy  = np.transpose(u2_kxy , axes=(2, 1, 0));  u2_kxy  = np.delete(u2_kxy , ignored_points_fld, axis = 0)
u2_kxy  = np.concatenate([u2_kxy [:, nky_mid:, :], u2_kxy [:, 0:nky_mid, :]], axis=1)
u2_kxy  = np.concatenate([u2_kxy [:, :, nkx_mid:], u2_kxy [:, :, 0:nkx_mid]], axis=2)
u2_kyz  = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/u2_kyz_sum_kx.dat').reshape(nt_fld, nkz, nky)); u2_kyz  = np.transpose(u2_kyz , axes=(2, 1, 0));  u2_kyz  = np.delete(u2_kyz , ignored_points_fld, axis = 0)
u2_kyz  = np.concatenate([u2_kyz [:, :, nky_mid:], u2_kyz [:, :, 0:nky_mid]], axis=2)
u2_kxz  = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/u2_kxz_sum_ky.dat').reshape(nt_fld, nkz, nkx)); u2_kxz  = np.transpose(u2_kxz , axes=(2, 1, 0));  u2_kxz  = np.delete(u2_kxz , ignored_points_fld, axis = 0)
u2_kxz  = np.concatenate([u2_kxz [:, :, nkx_mid:], u2_kxz [:, :, 0:nkx_mid]], axis=2)

b2_kxy  = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/b2_kxy_sum_kz.dat').reshape(nt_fld, nky, nkx)); b2_kxy  = np.transpose(b2_kxy , axes=(2, 1, 0));  b2_kxy  = np.delete(b2_kxy , ignored_points_fld, axis = 0)
b2_kxy  = np.concatenate([b2_kxy [:, nky_mid:, :], b2_kxy [:, 0:nky_mid, :]], axis=1)
b2_kxy  = np.concatenate([b2_kxy [:, :, nkx_mid:], b2_kxy [:, :, 0:nkx_mid]], axis=2)
b2_kyz  = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/b2_kyz_sum_kx.dat').reshape(nt_fld, nkz, nky)); b2_kyz  = np.transpose(b2_kyz , axes=(2, 1, 0));  b2_kyz  = np.delete(b2_kyz , ignored_points_fld, axis = 0)
b2_kyz  = np.concatenate([b2_kyz [:, :, nky_mid:], b2_kyz [:, :, 0:nky_mid]], axis=2)
b2_kxz  = np.transpose(np.fromfile(input_dir+'out2d'+restart_num+'/b2_kxz_sum_ky.dat').reshape(nt_fld, nkz, nkx)); b2_kxz  = np.transpose(b2_kxz , axes=(2, 1, 0));  b2_kxz  = np.delete(b2_kxz , ignored_points_fld, axis = 0)
b2_kxz  = np.concatenate([b2_kxz [:, :, nkx_mid:], b2_kxz [:, :, 0:nkx_mid]], axis=2)

# Modes removed by the 2/3 dealias filter are exactly zero, so log10 gives -inf
# there and u.min() in plot_2d would break the color scale. Take the color range
# from the finite values only; the -inf cells are left blank.
def log10_finite(u):
  with np.errstate(divide='ignore'):
    lu = np.log10(u)
  lu_finite = lu[np.isfinite(lu)]
  if lu_finite.size == 0:
    return lu, None, None
  return lu, lu_finite.min(), lu_finite.max()

# The kx and ky edges are dropped as in Calliope. For kz, only the Nyquist plane
# kz_max is dropped so that the kz = 0 plane stays in the plot.
lu, umin, umax = log10_finite(u2_kxy[final_fld_idx,1:-1,1:-1]); plot_2d(lu, kx_[1:-1], ky_[1:-1], umin=umin, umax=umax, xlab=r'$k_x$', ylab=r'$k_y$', title=r'$E_{u}    (t = %.2E) $' % tt_fld[final_fld_idx], save=outdir+'u2_kxy.pdf' )
lu, umin, umax = log10_finite(b2_kxy[final_fld_idx,1:-1,1:-1]); plot_2d(lu, kx_[1:-1], ky_[1:-1], umin=umin, umax=umax, xlab=r'$k_x$', ylab=r'$k_y$', title=r'$E_{b}    (t = %.2E) $' % tt_fld[final_fld_idx], save=outdir+'b2_kxy.pdf' )
# kyz and kxz have a single kz = 0 row in 2D runs, so they are plotted in 3D only
if not is2D:
  lu, umin, umax = log10_finite(u2_kyz[final_fld_idx,:-1,1:-1]); plot_2d(lu, ky_[1:-1], kz[:-1], umin=umin, umax=umax, xlab=r'$k_y$', ylab=r'$k_z$', title=r'$E_{u}    (t = %.2E) $' % tt_fld[final_fld_idx], save=outdir+'u2_kyz.pdf' )
  lu, umin, umax = log10_finite(u2_kxz[final_fld_idx,:-1,1:-1]); plot_2d(lu, kx_[1:-1], kz[:-1], umin=umin, umax=umax, xlab=r'$k_x$', ylab=r'$k_z$', title=r'$E_{u}    (t = %.2E) $' % tt_fld[final_fld_idx], save=outdir+'u2_kxz.pdf' )
  lu, umin, umax = log10_finite(b2_kyz[final_fld_idx,:-1,1:-1]); plot_2d(lu, ky_[1:-1], kz[:-1], umin=umin, umax=umax, xlab=r'$k_y$', ylab=r'$k_z$', title=r'$E_{b}    (t = %.2E) $' % tt_fld[final_fld_idx], save=outdir+'b2_kyz.pdf' )
  lu, umin, umax = log10_finite(b2_kxz[final_fld_idx,:-1,1:-1]); plot_2d(lu, kx_[1:-1], kz[:-1], umin=umin, umax=umax, xlab=r'$k_x$', ylab=r'$k_z$', title=r'$E_{b}    (t = %.2E) $' % tt_fld[final_fld_idx], save=outdir+'b2_kxz.pdf' )

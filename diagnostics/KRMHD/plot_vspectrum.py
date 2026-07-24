# -*- coding: utf-8 -*-
# Hermite (v_parallel) spectrum of the compressive fluctuations g.
# NOTE: the Alfvenic / EM fields (phi, psi, u_perp, b_perp) carry NO Hermite
# index m, so an m-spectrum exists ONLY for g. There is nothing EM to plot here.
from load import *
from fft import *
import sys
sys.path.append('../')
from plots import *

print('\nplotting vspectrum (Hermite m spectrum of g)\n')
outdir = './fig_vspectrum/'
import os
if not os.path.isdir(outdir):
  os.makedirs(outdir)

# m >= 1 for the log-log phase-mixing spectrum (m=0 is the density rung).
m1 = mm[1:]

# E_g(m): k-integrated per-moment free energy (exact, W_m) and, as an internal
# consistency check, the same quantity reconstructed from the binned spectrum
# g2_bin summed over (kz, kperp). The two agree up to the top-bin dealias cut.
Eg_m_Wm = W_m[final_idx, 1:]
Eg_m_g2 = np.sum(g2_bin[final_idx], axis=(0, 2))[1:]   # (kz,mm,kperp) -> mm

ys = [
       Eg_m_Wm,
       Eg_m_g2,
       m1**(-1./2.)/m1[0]**(-1./2.)*Eg_m_Wm[0],
       m1**(-1.)/m1[0]**(-1.)*Eg_m_Wm[0],
     ]
xs = [
       m1,
       m1,
       m1,
       m1,
     ]
ls = [
        '',
        'x',
        'k--',
        'k--',
     ]
legends = [
            r'$E_g(m) = W_m$',
            r'$\sum_{k_z,k_\+}\tfrac12|g_m|^2$',
            r'$m^{-1/2}$',
            r'$m^{-1}$',
          ]
plot_log1d_many(xs, ys, xlab=mlab, legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt[final_idx], ylab='', term=True, save=outdir+'m_spectra_g.pdf')

# Hermite flux Gamma_m(m) (eq.9), k-integrated. Only if it was written.
# In the inertial range of stochastic echo Gamma_m ~ 0; for linear phase mixing
# it is roughly constant and non-zero.
if has_hermite_flux:
  ys = [
         Gamma_m_kint[final_idx, 1:],
       ]
  xs = [
         m1,
       ]
  ls = [
         '',
       ]
  legends = [
              r'$\Gamma_m$',
            ]
  plot_semilogx1d_many(xs, ys, xlab=mlab, legends=legends, ls=ls, legendloc='upper right', title=r'$t = %.2E $' % tt[final_idx], ylab='', term=True, save=outdir+'Gamma_m.pdf')
else:
  print('  (Gamma_m not present in nc: write_hermite_flux was off -- skipping flux plot)')

#------------------#
#   output ascii   #
#------------------#
if has_hermite_flux:
  np.savetxt(outdir + 'Em.txt', np.column_stack((mm,
                                                  W_m[final_idx, :],
                                                  np.sum(g2_bin[final_idx], axis=(0, 2)),
                                                  Gamma_m_kint[final_idx, :],
                                                )), fmt='%E')
else:
  np.savetxt(outdir + 'Em.txt', np.column_stack((mm,
                                                  W_m[final_idx, :],
                                                  np.sum(g2_bin[final_idx], axis=(0, 2)),
                                                )), fmt='%E')

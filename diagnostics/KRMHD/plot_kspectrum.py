# -*- coding: utf-8 -*-
from load import *
from fft import *
import sys
sys.path.append('../')
from plots import *

print('\nplotting kspectrum\n')
outdir = './fig_kspectrum/'
import os
if not os.path.isdir(outdir):
  os.makedirs(outdir)

if nlz == nkz:
  kp_end = np.argmin(np.abs(kpbin - kpbin.max()*2./3.))
  if not is2D:
    kz_end = nkz
else:
  kp_end = int(kpbin.size*2./3.)
  kz_end = int(nkz*2./3.)

#--------------------------------------------------------#
#                      plot 1D spectra                   #
#--------------------------------------------------------#
# kprp spectrum (Alfvenic / EM)
ys = [
       np.sum(upe2_bin[final_idx, :, 1:kp_end], axis=0),
       np.sum(bpe2_bin[final_idx, :, 1:kp_end], axis=0),
       kpbin[1:kp_end]**(-5./3.)/kpbin[1]**(-5./3.)*np.sum(bpe2_bin[final_idx,:,1:kp_end], axis=0)[0],
       kpbin[1:kp_end]**(-3./2.)/kpbin[1]**(-3./2.)*np.sum(bpe2_bin[final_idx,:,1:kp_end], axis=0)[0]
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
            r'$E_{u_\+}$',
            r'$E_{\delta B_\+}$',
            r'-5/3',
            r'-3/2',
          ]
plot_log1d_many(xs, ys, xlab='$k_\+ L_\+$', legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt[final_idx], ylab='', term=True, save=outdir+'kprp_spectra.pdf')

## Elsasser fields
ys = [
       np.sum(zppe2_bin[final_idx, :, 1:kp_end], axis=0),
       np.sum(zmpe2_bin[final_idx, :, 1:kp_end], axis=0),
       kpbin[1:kp_end]**(-5./3.)/kpbin[1]**(-5./3.)*np.sum(bpe2_bin[final_idx,:,1:kp_end], axis=0)[0],
       kpbin[1:kp_end]**(-3./2.)/kpbin[1]**(-3./2.)*np.sum(bpe2_bin[final_idx,:,1:kp_end], axis=0)[0]
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
            r'$E_{Z_\+^+}$',
            r'$E_{Z_\+^-}$',
            r'-5/3',
            r'-3/2',
          ]
plot_log1d_many(xs, ys, xlab='$k_\+ L_\+$', legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt[final_idx], ylab='', term=True, save=outdir+'kprp_spectra_ELS.pdf')

# kprp spectrum (compressive g): E_g(kprp) = sum_{kz,m} 1/2 |g_m|^2
Eg_kprp = np.sum(g2_bin[final_idx, :, :, 1:kp_end], axis=(0, 1))
ys = [
       Eg_kprp,
       kpbin[1:kp_end]**(-5./3.)/kpbin[1]**(-5./3.)*Eg_kprp[0],
       kpbin[1:kp_end]**(-3./2.)/kpbin[1]**(-3./2.)*Eg_kprp[0]
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
            r'$E_{g}$',
            r'-5/3',
            r'-3/2',
          ]
plot_log1d_many(xs, ys, xlab='$k_\+ L_\+$', legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt[final_idx], ylab='', term=True, save=outdir+'kprp_spectra_g.pdf')

# kz spectrum (Alfvenic / EM)
if not is2D:
  ys = [
         np.sum(upe2_bin[final_idx, 1:kz_end, :kp_end], axis=1),
         np.sum(bpe2_bin[final_idx, 1:kz_end, :kp_end], axis=1),
       ]
  xs = [
          kz[1:kz_end],
          kz[1:kz_end],
       ]
  ls = [
          '',
          '',
       ]
  legends = [
              r'$E_{u_\+}$',
              r'$E_{\delta B_\+}$',
            ]
  plot_log1d_many(xs, ys, xlab=kzlab, legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt[final_idx], ylab='', term=True, save=outdir+'kz_spectra.pdf')

## Elsasser fields
if not is2D:
  ys = [
         np.sum(zppe2_bin[final_idx, 1:kz_end, :kp_end], axis=1),
         np.sum(zmpe2_bin[final_idx, 1:kz_end, :kp_end], axis=1),
       ]
  xs = [
          kz[1:kz_end],
          kz[1:kz_end],
       ]
  ls = [
          '',
          '',
       ]
  legends = [
              r'$E_{Z_\+^+}$',
              r'$E_{Z_\+^-}$',
            ]
  plot_log1d_many(xs, ys, xlab=kzlab, legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt[final_idx], ylab='', term=True, save=outdir+'kz_spectra_ELS.pdf')

## kz spectrum (compressive g): E_g(kz) = sum_{kprp,m} 1/2 |g_m|^2
if not is2D:
  Eg_kz = np.sum(g2_bin[final_idx, 1:kz_end, :, :kp_end], axis=(1, 2))
  ys = [
         Eg_kz,
         kz[1:kz_end]**(-2.)/kz[1]**(-2.)*Eg_kz[0],
       ]
  xs = [
          kz[1:kz_end],
          kz[1:kz_end],
       ]
  ls = [
          '',
          'k--',
       ]
  legends = [
              r'$E_{g}$',
              r'-2',
            ]
  plot_log1d_many(xs, ys, xlab=kzlab, legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt[final_idx], ylab='', term=True, save=outdir+'kz_spectra_g.pdf')

#--------------------------------------------------------#
#                      plot 2D spectra                   #
#--------------------------------------------------------#
if not is2D:
  plot_log2d(upe2_bin[final_idx, 1:kz_end, 1:kp_end], kpbin[1:kp_end], kz[1:kz_end], xlab='$k_\+ L_\+$', ylab=kzlab,
      title=r'$E_{u_{\+}}$' + ' $(t = $ %.2E' % tt[final_idx] + '$)$', save=outdir + 'upe2.pdf')
  plot_log2d(bpe2_bin[final_idx, 1:kz_end, 1:kp_end], kpbin[1:kp_end], kz[1:kz_end], xlab='$k_\+ L_\+$', ylab=kzlab,
      title=r'$E_{\delta B_\+}$' + ' $(t = $ %.2E' % tt[final_idx] + '$)$', save=outdir + 'bpe2.pdf')
  plot_log2d(zppe2_bin[final_idx, 1:kz_end, 1:kp_end], kpbin[1:kp_end], kz[1:kz_end], xlab='$k_\+ L_\+$', ylab=kzlab,
      title=r'$E_{Z_{\+}^+}$' + ' $(t = $ %.2E' % tt[final_idx] + '$)$', save=outdir + 'zppe2.pdf')
  plot_log2d(zmpe2_bin[final_idx, 1:kz_end, 1:kp_end], kpbin[1:kp_end], kz[1:kz_end], xlab='$k_\+ L_\+$', ylab=kzlab,
      title=r'$E_{Z_{\+}^-}$' + ' $(t = $ %.2E' % tt[final_idx] + '$)$', save=outdir + 'zmpe2.pdf')
  # compressive g, summed over Hermite index m
  g2_2d = np.sum(g2_bin[final_idx, 1:kz_end, :, 1:kp_end], axis=1)
  plot_log2d(g2_2d, kpbin[1:kp_end], kz[1:kz_end], xlab='$k_\+ L_\+$', ylab=kzlab,
      title=r'$E_{g}$' + ' $(t = $ %.2E' % tt[final_idx] + '$)$', save=outdir + 'g2.pdf')

#------------------#
#   output ascii   #
#------------------#
np.savetxt(outdir + 'Ekprp.txt'  , np.column_stack((kpbin[:kp_end],
                                                     np.sum(upe2_bin [final_idx,:kz_end,:kp_end], axis=0),
                                                     np.sum(bpe2_bin [final_idx,:kz_end,:kp_end], axis=0),
                                                     np.sum(zppe2_bin[final_idx,:kz_end,:kp_end], axis=0),
                                                     np.sum(zmpe2_bin[final_idx,:kz_end,:kp_end], axis=0),
                                                     np.sum(g2_bin   [final_idx,:kz_end,:,:kp_end], axis=(0, 1)),
                                                   )), fmt='%E')
if not is2D:
  np.savetxt(outdir + 'Ekz.txt'  , np.column_stack((kz[:kz_end],
                                                         np.sum(upe2_bin [final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(bpe2_bin [final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(zppe2_bin[final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(zmpe2_bin[final_idx,:kz_end,:kp_end], axis=1),
                                                         np.sum(g2_bin   [final_idx,:kz_end,:,:kp_end], axis=(1, 2)),
                                                       )), fmt='%E')

del upe2_bin
del bpe2_bin

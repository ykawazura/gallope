# -*- coding: utf-8 -*-
# Time-averaged diagnostics over a user-specified interval.
# The interval is specified by tt VALUE (avg_start_t, avg_end_t), not index.
# Produces averaged EM power balance / energy, EM & g k-spectra, and the g
# Hermite m-spectrum, together with their standard deviations.
from load import *
from fft import *
import sys
sys.path.append('../')
from plots import *
from scipy.integrate import trapezoid as trapz  # scipy>=1.14 removed the legacy `trapz` alias
from scipy import interpolate
import os

#--------------------------------------------------------#
#             averaging window (by tt value)             #
#--------------------------------------------------------#
# Set avg_end_t = None to average up to the last record.
avg_start_t = 0.0
avg_end_t   = None

avg_start = int(np.argmin(np.abs(tt - avg_start_t)))
if avg_end_t is None:
  avg_end = tt.size - 1
else:
  avg_end = int(np.argmin(np.abs(tt - avg_end_t)))
# inclusive slice; guarantee at least two points for the trapezoidal rule
if avg_end <= avg_start:
  avg_end = min(avg_start + 1, tt.size - 1)
avg = slice(avg_start, avg_end + 1)


def time_average(x, y, axis=0): # x: 1D array, y: any-D array
  return trapz(y, x, axis=axis)/(x[-1] - x[0])


def std_dev(x, y): # x, y: 1D array
  y_intp = interpolate.interp1d(x, y)
  return np.std(y_intp(np.linspace(x[0], x[-1], x.size)), ddof=0)


##########################################################
#              average energy time evolution             #
##########################################################
print('\nplotting energy\n')
outdir = './fig_energy/'
if not os.path.isdir(outdir):
  os.makedirs(outdir)

# Alfvenic (EM) power balance quantities (dW/dt = P - D)
Waw      = upe2_sum + bpe2_sum
Waw_dot  = upe2dot_sum + bpe2dot_sum
Daw      = upe2dissip_sum + bpe2dissip_sum
Paw      = p_phi_sum + p_psi_sum   # NB: RMHD time_average used a nonexistent p_omg_sum

Waw_avg      = time_average(tt[avg], Waw     [avg])
Waw_dot_avg  = time_average(tt[avg], Waw_dot [avg])
Daw_avg      = time_average(tt[avg], Daw     [avg])
Paw_avg      = time_average(tt[avg], Paw     [avg])
upe2_sum_avg = time_average(tt[avg], upe2_sum[avg])
bpe2_sum_avg = time_average(tt[avg], bpe2_sum[avg])
W_free_avg   = time_average(tt[avg], W_free  [avg])

Waw_err      = std_dev(tt[avg], Waw     [avg])
Waw_dot_err  = std_dev(tt[avg], Waw_dot [avg])
Daw_err      = std_dev(tt[avg], Daw     [avg])
Paw_err      = std_dev(tt[avg], Paw     [avg])
upe2_sum_err = std_dev(tt[avg], upe2_sum[avg])
bpe2_sum_err = std_dev(tt[avg], bpe2_sum[avg])
W_free_err   = std_dev(tt[avg], W_free  [avg])

s =     'average over t     = [%.3E' % tt[avg_start] + ', %.3E' % tt[avg_end] + ']' + '\n'
s = s + 'average over index = [' + str(avg_start) + ', ' + str(avg_end) + ']' + '\n\n'
s = s + '  [Alfvenic / EM]' + '\n'
s = s + '  Waw              = %.3E \pm %.3E'  % (Waw_avg      , Waw_err      ) + '\n'
s = s + '    (upe2          = %.3E \pm %.3E)' % (upe2_sum_avg , upe2_sum_err ) + '\n'
s = s + '    (bpe2          = %.3E \pm %.3E)' % (bpe2_sum_avg , bpe2_sum_err ) + '\n'
s = s + '  Waw_dot          = %.3E \pm %.3E'  % (Waw_dot_avg  , Waw_dot_err  ) + '\n'
s = s + '  Paw              = %.3E \pm %.3E'  % (Paw_avg      , Paw_err      ) + '\n'
s = s + '  Daw              = %.3E \pm %.3E'  % (Daw_avg      , Daw_err      ) + '\n\n'
s = s + '  [compressive / g]' + '\n'
s = s + '  W_free           = %.3E \pm %.3E'  % (W_free_avg   , W_free_err   ) + '\n'
print (s)

f = open('time_average.txt', 'w')
f.write(s)
f.close()

# averaged EM power balance
ys = [
       upe2dot_sum + bpe2dot_sum,
       upe2dissip_sum + bpe2dissip_sum,
       -p_phi_sum,
       -p_psi_sum,
       upe2dot_sum + bpe2dot_sum + upe2dissip_sum + bpe2dissip_sum - p_phi_sum - p_psi_sum,
       np.full([tt[avg].size], Waw_dot_avg),
       np.full([tt[avg].size], Daw_avg),
       np.full([tt[avg].size], -Paw_avg),
     ]
xs = [
       tt, tt, tt, tt, tt,
       tt[avg], tt[avg], tt[avg],
     ]
ls = [
        '', '', '', '', 'k--', '', '', '',
     ]
legends = [
       r'$\rmd W/\rmd t$',
       r'$D_\mr{AW}$',
       r'$-P_{\Phi}$',
       r'$-P_{\Psi}$',
       r'balance',
       '', '', '',
     ]
plot_1d_many_average(xs, ys, tt[avg_start], tt[avg_end], xlab=tlab, legends=legends, ls=ls, legendloc='upper left', title='', ylab='', term=True, save=outdir + 'balance_all_avg.pdf')

# averaged energies (EM + compressive free energy)
ys = [
       upe2_sum,
       bpe2_sum,
       W_free,
       np.full([tt[avg].size], upe2_sum_avg),
       np.full([tt[avg].size], bpe2_sum_avg),
       np.full([tt[avg].size], W_free_avg),
     ]
ls = [
        '', '', '', '', '', '',
     ]
xs = [
       tt, tt, tt,
       tt[avg], tt[avg], tt[avg],
     ]
legends = [
       r'$W_{u_\+}$',
       r'$W_{\delta B_\+}$',
       r'$W_\mr{free}$',
       '', '', '',
     ]
plot_1d_many_average(xs, ys, tt[avg_start], tt[avg_end], xlab=tlab, legends=legends, ls=ls, legendloc='upper left', title='', ylab='', term=True, save=outdir + 'energy_all_avg.pdf')


##########################################################
#                    average kspectrum                   #
##########################################################
print('\nplotting kspectrum\n')
outdir = './fig_kspectrum/'
if not os.path.isdir(outdir):
  os.makedirs(outdir)

# dealias cut-offs (same convention as plot_kspectrum.py: no negative-kz folding,
# because the KRMHD nc stores the non-negative kz half with the factor-of-2 already
# applied by the Fortran diagnostics).
if nlz == nkz:
  kp_end = np.argmin(np.abs(kpbin - kpbin.max()*2./3.))
  if not is2D:
    kz_end = nkz
else:
  kp_end = int(kpbin.size*2./3.)
  kz_end = int(nkz*2./3.)

# time-average the binned spectra
upe2_bin_avg = time_average(tt[avg], upe2_bin[avg], axis=0)   # (nkz, nkpolar)
bpe2_bin_avg = time_average(tt[avg], bpe2_bin[avg], axis=0)
g2_bin_avg   = time_average(tt[avg], g2_bin  [avg], axis=0)   # (nkz, nm, nkpolar)

# kprp spectra (EM)
ys = [
       np.sum(upe2_bin_avg[:, 1:kp_end], axis=0),
       np.sum(bpe2_bin_avg[:, 1:kp_end], axis=0),
       kpbin[1:kp_end]**(-5./3.)/kpbin[1]**(-5./3.)*np.sum(bpe2_bin_avg[:, 1:kp_end], axis=0)[0],
     ]
xs = [
       kpbin[1:kp_end], kpbin[1:kp_end], kpbin[1:kp_end],
     ]
ls = [
        '', '', 'k--',
     ]
legends = [
            r'$E_{u_\+}$',
            r'$E_{\delta B_\+}$',
            r'-5/3',
          ]
plot_log1d_many(xs, ys, xlab='$k_\+ L_\+$', legends=legends, ls=ls, legendloc='lower left', ylab='', term=True, save=outdir+'kprp_spectra_avg.pdf')

# kprp spectrum (compressive g): E_g(kprp) = sum_{kz,m} 1/2 |g_m|^2
Eg_kprp_avg = np.sum(g2_bin_avg[:, :, 1:kp_end], axis=(0, 1))
ys = [
       Eg_kprp_avg,
       kpbin[1:kp_end]**(-5./3.)/kpbin[1]**(-5./3.)*Eg_kprp_avg[0],
       kpbin[1:kp_end]**(-3./2.)/kpbin[1]**(-3./2.)*Eg_kprp_avg[0],
     ]
xs = [
       kpbin[1:kp_end], kpbin[1:kp_end], kpbin[1:kp_end],
     ]
ls = [
        '', 'k--', 'k-.',
     ]
legends = [
            r'$E_{g}$',
            r'-5/3',
            r'-3/2',
          ]
plot_log1d_many(xs, ys, xlab='$k_\+ L_\+$', legends=legends, ls=ls, legendloc='lower left', ylab='', term=True, save=outdir+'kprp_spectra_g_avg.pdf')

# kz spectra (EM + g)
if not is2D:
  ys = [
         np.sum(upe2_bin_avg[1:kz_end, :kp_end], axis=1),
         np.sum(bpe2_bin_avg[1:kz_end, :kp_end], axis=1),
       ]
  xs = [
          kz[1:kz_end], kz[1:kz_end],
       ]
  ls = [
          '', '',
       ]
  legends = [
              r'$E_{u_\+}$',
              r'$E_{\delta B_\+}$',
            ]
  plot_log1d_many(xs, ys, xlab=kzlab, legends=legends, ls=ls, legendloc='lower left', ylab='', term=True, save=outdir+'kz_spectra_avg.pdf')

  Eg_kz_avg = np.sum(g2_bin_avg[1:kz_end, :, :kp_end], axis=(1, 2))
  ys = [
         Eg_kz_avg,
         kz[1:kz_end]**(-2.)/kz[1]**(-2.)*Eg_kz_avg[0],
       ]
  xs = [
          kz[1:kz_end], kz[1:kz_end],
       ]
  ls = [
          '', 'k--',
       ]
  legends = [
              r'$E_{g}$',
              r'-2',
            ]
  plot_log1d_many(xs, ys, xlab=kzlab, legends=legends, ls=ls, legendloc='lower left', ylab='', term=True, save=outdir+'kz_spectra_g_avg.pdf')

#--------------------------------------------------------#
#                      plot 2D spectra                   #
#--------------------------------------------------------#
if not is2D:
  plot_log2d(upe2_bin_avg[1:kz_end, 1:kp_end], kpbin[1:kp_end], kz[1:kz_end], xlab='$k_\+ L_\+$', ylab=kzlab,
      title=r'$E_{u_{\+}}$', save=outdir + 'upe2_avg.pdf')
  plot_log2d(bpe2_bin_avg[1:kz_end, 1:kp_end], kpbin[1:kp_end], kz[1:kz_end], xlab='$k_\+ L_\+$', ylab=kzlab,
      title=r'$E_{\delta B_{\+}}$', save=outdir + 'bpe2_avg.pdf')
  g2_2d_avg = np.sum(g2_bin_avg[1:kz_end, :, 1:kp_end], axis=1)
  plot_log2d(g2_2d_avg, kpbin[1:kp_end], kz[1:kz_end], xlab='$k_\+ L_\+$', ylab=kzlab,
      title=r'$E_{g}$', save=outdir + 'g2_avg.pdf')

#------------------#
#   output ascii   #
#------------------#
np.savetxt(outdir + 'Ekprp_avg.txt'  , np.column_stack((kpbin[:kp_end],
                                                     np.sum(upe2_bin_avg[:, :kp_end], axis=0),
                                                     np.sum(bpe2_bin_avg[:, :kp_end], axis=0),
                                                     np.sum(g2_bin_avg  [:, :, :kp_end], axis=(0, 1)),
                                                   )), fmt='%E')
if not is2D:
  np.savetxt(outdir + 'Ekz_avg.txt'  , np.column_stack((kz[:kz_end],
                                                     np.sum(upe2_bin_avg[:kz_end, :], axis=1),
                                                     np.sum(bpe2_bin_avg[:kz_end, :], axis=1),
                                                     np.sum(g2_bin_avg  [:kz_end, :, :], axis=(1, 2)),
                                                   )), fmt='%E')


##########################################################
#              average vspectrum (Hermite m)             #
##########################################################
print('\nplotting vspectrum (Hermite m spectrum of g)\n')
outdir = './fig_vspectrum/'
if not os.path.isdir(outdir):
  os.makedirs(outdir)

# k-integrated per-moment free energy, time-averaged
W_m_avg   = time_average(tt[avg], W_m[avg], axis=0)          # exact (nm,)
Eg_m_avg  = np.sum(g2_bin_avg, axis=(0, 2))                  # from binned spectrum (nm,)

m1 = mm[1:]
ys = [
       W_m_avg[1:],
       Eg_m_avg[1:],
       m1**(-1./2.)/m1[0]**(-1./2.)*W_m_avg[1],
       m1**(-1.)/m1[0]**(-1.)*W_m_avg[1],
     ]
xs = [
       m1, m1, m1, m1,
     ]
ls = [
        '', 'x', 'k--', 'k--',
     ]
legends = [
            r'$E_g(m) = W_m$',
            r'$\sum_{k_z,k_\+}\tfrac12|g_m|^2$',
            r'$m^{-1/2}$',
            r'$m^{-1}$',
          ]
plot_log1d_many(xs, ys, xlab=mlab, legends=legends, ls=ls, legendloc='lower left', ylab='', term=True, save=outdir+'m_spectra_g_avg.pdf')

np.savetxt(outdir + 'Em_avg.txt', np.column_stack((mm, W_m_avg, Eg_m_avg)), fmt='%E')


del upe2_bin
del bpe2_bin
del g2_bin

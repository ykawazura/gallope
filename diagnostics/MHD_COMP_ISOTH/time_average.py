# -*- coding: utf-8 -*-
from load import *
from fft import *
import sys
sys.path.append('../')
from plots import *
from scipy import interpolate
try:                                    # numpy >= 2.0
  from numpy import trapezoid
except ImportError:                     # numpy < 2.0
  from numpy import trapz as trapezoid

avg_start = 0
avg_end   = -1

def time_average(x, y, axis=0): # x: 1D array, y: any-D array
  return trapezoid(y, x, axis=axis)/(x[-1] - x[0])

def std_dev(x, y): # x, y: 1D array
  y_intp = interpolate.interp1d(x, y)
  return np.std(y_intp(np.linspace(x[0], x[-1], x.size)), ddof=0)

def cov(x, y1, y2): # x, y: 1D array
  y1_intp  = interpolate.interp1d(x, y1)
  y2_intp  = interpolate.interp1d(x, y2)
  y1_intp_ = y1_intp(np.linspace(x[0], x[-1], x.size))
  y2_intp_ = y2_intp(np.linspace(x[0], x[-1], x.size))

  return np.cov( np.stack((y1_intp_, y2_intp_), axis=0) )

if avg_end == -1:
  avg_end = np.argwhere(tt == tt[avg_end])[0][0]

##########################################################
#              average energy time evolution             #
##########################################################
print('\nplotting energy\n')
outdir = './fig_energy/'

W        = wkin_sum + wmag_sum + wrho_sum
W_dot    = wkin_dot_sum + wmag_dot_sum + wrho_dot_sum
D        = wkin_dissip_sum + wmag_dissip_sum + wrho_dissip_sum
P        = p_ext_sum + p_re_sum + p_ma_sum

W_avg        = time_average(tt[avg_start:avg_end], W       [avg_start:avg_end], axis=0)
W_dot_avg    = time_average(tt[avg_start:avg_end], W_dot   [avg_start:avg_end], axis=0)
D_avg        = time_average(tt[avg_start:avg_end], D       [avg_start:avg_end], axis=0)
P_avg        = time_average(tt[avg_start:avg_end], P       [avg_start:avg_end], axis=0)
wkin_sum_avg = time_average(tt[avg_start:avg_end], wkin_sum[avg_start:avg_end], axis=0)
wmag_sum_avg = time_average(tt[avg_start:avg_end], wmag_sum[avg_start:avg_end], axis=0)
wrho_sum_avg = time_average(tt[avg_start:avg_end], wrho_sum[avg_start:avg_end], axis=0)

W_err        = std_dev(tt[avg_start:avg_end], W       [avg_start:avg_end])
W_dot_err    = std_dev(tt[avg_start:avg_end], W_dot   [avg_start:avg_end])
D_err        = std_dev(tt[avg_start:avg_end], D       [avg_start:avg_end])
P_err        = std_dev(tt[avg_start:avg_end], P       [avg_start:avg_end])
wkin_sum_err = std_dev(tt[avg_start:avg_end], wkin_sum[avg_start:avg_end])
wmag_sum_err = std_dev(tt[avg_start:avg_end], wmag_sum[avg_start:avg_end])
wrho_sum_err = std_dev(tt[avg_start:avg_end], wrho_sum[avg_start:avg_end])

s =     'average over t     = [%.3E' % tt[avg_start] + ', %.3E' % tt[avg_end] + ']' + '\n'
s = s + 'average over index = [' + str(avg_start) + ', ' + str(avg_end) + ']' + '\n'
s = s + ' error of ratio is calculated by (a + da)/(b + db) ~ a/b*[1 + sqrt( (da/a)^2 + (db/b)^2 - 2cov/(a*b))]' + '\n\n'
s = s + r'  W     = %.3E \pm %.3E'  % (W_avg    , W_err    ) + '\n'
s = s + r'  W_dot = %.3E \pm %.3E'  % (W_dot_avg, W_dot_err) + '\n'
s = s + r'  P     = %.3E \pm %.3E'  % (P_avg    , P_err    ) + '\n'
s = s + r'  D     = %.3E \pm %.3E'  % (D_avg    , D_err    ) + '\n'
print (s)

f = open('time_average.txt','w')
f.write(s)
f.close()


ys = [
       wkin_dot_sum + wmag_dot_sum + wrho_dot_sum,
       wkin_dissip_sum + wmag_dissip_sum + wrho_dissip_sum,
       -p_ext_sum,
       -p_re_sum,
       -p_ma_sum,
       wkin_dot_sum + wmag_dot_sum + wrho_dot_sum + wkin_dissip_sum + wmag_dissip_sum + wrho_dissip_sum - p_ext_sum - p_re_sum - p_ma_sum,
       np.full([tt[avg_start:avg_end].size], W_dot_avg),
       np.full([tt[avg_start:avg_end].size], D_avg),
       np.full([tt[avg_start:avg_end].size], -P_avg),
     ]
xs = [
       tt,
       tt,
       tt,
       tt,
       tt,
       tt,
       tt[avg_start:avg_end],
       tt[avg_start:avg_end],
       tt[avg_start:avg_end],
     ]
ls = [
        '',
        '',
        '',
        '',
        '',
        'k--',
        '',
        '',
        '',
     ]
legends = [
       r'$\mathrm{d}W/\mathrm{d} t$',
       r'$D$',
       r'$P_\mr{ext}$',
       r'$P_\mr{Re}$',
       r'$P_\mr{M}$',
       r'balance',
       '',
       '',
       '',
     ]
plot_1d_many_average(xs, ys, tt[avg_start], tt[avg_end], xlab='$'+tlab+'$', legends=legends, ls=ls, legendloc='upper left', title='', ylab='', term=True, save=outdir + 'balance_all_avg.pdf')


ys = [
       wkin_sum,
       wmag_sum,
       wrho_sum,
       np.full([tt[avg_start:avg_end].size], wkin_sum_avg),
       np.full([tt[avg_start:avg_end].size], wmag_sum_avg),
       np.full([tt[avg_start:avg_end].size], wrho_sum_avg),
     ]
ls = [
        '',
        '',
        '',
        '',
        '',
        '',
     ]
xs = [
       tt,
       tt,
       tt,
       tt[avg_start:avg_end],
       tt[avg_start:avg_end],
       tt[avg_start:avg_end],
     ]
legends = [
       r'$W_\mr{kin}$',
       r'$W_\mr{mag}$',
       r'$W_\rho$',
       '',
       '',
       '',
     ]
plot_1d_many_average(xs, ys, tt[avg_start], tt[avg_end], xlab='$'+tlab+'$', legends=legends, ls=ls, legendloc='upper left', title='', ylab='', term=True, save=outdir + 'energy_all_avg.pdf')

##########################################################
#                    average kspectrum                   #
##########################################################
print('\nplotting kspectrum\n')
outdir = './fig_kspectrum/'

u2_bin_avg  = time_average(tt[avg_start:avg_end], u2_bin [avg_start:avg_end], axis=0)
ux2_bin_avg = time_average(tt[avg_start:avg_end], ux2_bin[avg_start:avg_end], axis=0)
uy2_bin_avg = time_average(tt[avg_start:avg_end], uy2_bin[avg_start:avg_end], axis=0)
uz2_bin_avg = time_average(tt[avg_start:avg_end], uz2_bin[avg_start:avg_end], axis=0)
b2_bin_avg  = time_average(tt[avg_start:avg_end], b2_bin [avg_start:avg_end], axis=0)
bx2_bin_avg = time_average(tt[avg_start:avg_end], bx2_bin[avg_start:avg_end], axis=0)
by2_bin_avg = time_average(tt[avg_start:avg_end], by2_bin[avg_start:avg_end], axis=0)
bz2_bin_avg = time_average(tt[avg_start:avg_end], bz2_bin[avg_start:avg_end], axis=0)
rho_bin_avg = time_average(tt[avg_start:avg_end], rho_bin[avg_start:avg_end], axis=0)
zp2_bin_avg = time_average(tt[avg_start:avg_end], zp2_bin[avg_start:avg_end], axis=0)
zm2_bin_avg = time_average(tt[avg_start:avg_end], zm2_bin[avg_start:avg_end], axis=0)

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
  ys.append([np.min([u2_bin_avg[1:kp_end].min(), b2_bin_avg[1:kp_end].min()]),
             np.max([u2_bin_avg[1:kp_end].max(), b2_bin_avg[1:kp_end].max()])])
  ls.append('k:')
  legends.append(r'$k = k_\mr{MRI} = (\sqrt{15}/4)\varpi_0/v_\rmA$')

  return xs, ys, ls, legends

#--------------------------------------------------------#
#                      plot 1D spectra                   #
#--------------------------------------------------------#
ys = [
       u2_bin_avg[1:kp_end],
       b2_bin_avg[1:kp_end],
       kpbin[1:kp_end]**(-5./3.)/kpbin[1]**(-5./3.)*b2_bin_avg[1:kp_end][0],
       kpbin[1:kp_end]**(-3./2.)/kpbin[1]**(-3./2.)*b2_bin_avg[1:kp_end][0]
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

plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', ylab='', term=True, save=outdir+'k_spectra_avg.pdf')


# u by components
ys = [
       ux2_bin_avg[1:kp_end],
       uy2_bin_avg[1:kp_end],
       uz2_bin_avg[1:kp_end],
       kpbin[1:kp_end]**(-5./3.)/kpbin[1]**(-5./3.)*u2_bin_avg[1:kp_end][0],
       kpbin[1:kp_end]**(-3./2.)/kpbin[1]**(-3./2.)*u2_bin_avg[1:kp_end][0],
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

plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', ylab='', term=True, save=outdir+'k_spectra_u_avg.pdf')


# B by components
ys = [
       bx2_bin_avg[1:kp_end],
       by2_bin_avg[1:kp_end],
       bz2_bin_avg[1:kp_end],
       kpbin[1:kp_end]**(-5./3.)/kpbin[1]**(-5./3.)*b2_bin_avg[1:kp_end][0],
       kpbin[1:kp_end]**(-3./2.)/kpbin[1]**(-3./2.)*b2_bin_avg[1:kp_end][0],
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

plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', ylab='', term=True, save=outdir+'k_spectra_b_avg.pdf')


# Elsasser fields
ys = [
       zp2_bin_avg[1:kp_end],
       zm2_bin_avg[1:kp_end],
       kpbin[1:kp_end]**(-5./3.)/kpbin[1]**(-5./3.)*zp2_bin_avg[1:kp_end][0],
       kpbin[1:kp_end]**(-3./2.)/kpbin[1]**(-3./2.)*zp2_bin_avg[1:kp_end][0],
     ]
xs = [
      kpbin[1:kp_end],
      kpbin[1:kp_end],
      kpbin[1:kp_end],
      kpbin[1:kp_end],
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

plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', ylab='', term=True, save=outdir+'k_spectra_ELS_avg.pdf')

#------------------#
#   output ascii   #
#------------------#
np.savetxt(outdir + 'Ek_avg.txt'  , np.column_stack((kpbin[:kp_end],
                                                      u2_bin_avg[:kp_end],
                                                     ux2_bin_avg[:kp_end],
                                                     uy2_bin_avg[:kp_end],
                                                     uz2_bin_avg[:kp_end],
                                                      b2_bin_avg[:kp_end],
                                                     bx2_bin_avg[:kp_end],
                                                     by2_bin_avg[:kp_end],
                                                     bz2_bin_avg[:kp_end],
                                                     rho_bin_avg[:kp_end],
                                                     )), fmt='%E')

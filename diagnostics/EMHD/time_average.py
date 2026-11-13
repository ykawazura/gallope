# -*- coding: utf-8 -*-
from load import *
from fft import *
import sys
sys.path.append('../')
from plots import *
from scipy.integrate import simps, trapz
from scipy import interpolate

avg_start = 0
avg_end   = -1

@jit
def time_average(x, y, axis=0): # x: 1D array, y: any-D array 
  return trapz(y, x, axis=axis)/(x[-1] - x[0])
  # return trapz(y, x, axis=axis)/(x[-1] - x[0])

@jit
def std_dev(x, y): # x, y: 1D array
  y_intp = interpolate.interp1d(x, y)
  return np.std(y_intp(np.linspace(x[0], x[-1], x.size)), ddof=0)

@jit
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

W        = wmag_sum
W_dot    = wmag_dot_sum
D        = wmag_dissip_sum
P        = p_ext_sum

W_avg        = time_average(tt[avg_start:avg_end], W       [avg_start:avg_end], axis=0)
W_dot_avg    = time_average(tt[avg_start:avg_end], W_dot   [avg_start:avg_end], axis=0)
D_avg        = time_average(tt[avg_start:avg_end], D       [avg_start:avg_end], axis=0)
P_avg        = time_average(tt[avg_start:avg_end], P       [avg_start:avg_end], axis=0)
wmag_sum_avg = time_average(tt[avg_start:avg_end], wmag_sum[avg_start:avg_end], axis=0)
             
W_err        = std_dev(tt[avg_start:avg_end], W       [avg_start:avg_end])
W_dot_err    = std_dev(tt[avg_start:avg_end], W_dot   [avg_start:avg_end])
D_err        = std_dev(tt[avg_start:avg_end], D       [avg_start:avg_end])
P_err        = std_dev(tt[avg_start:avg_end], P       [avg_start:avg_end])
wmag_sum_err = std_dev(tt[avg_start:avg_end], wmag_sum[avg_start:avg_end])

s =     'average over t     = [%.3E' % tt[avg_start] + ', %.3E' % tt[avg_end] + ']' + '\n'
s = s + 'average over index = [' + str(avg_start) + ', ' + str(avg_end) + ']' + '\n'
s = s + ' error of ratio is calculated by (a + da)/(b + db) ~ a/b*[1 + sqrt( (da/a)^2 + (db/b)^2 - 2cov/(a*b))]' + '\n\n' 
s = s + '  W     = %.3E \pm %.3E'  % (W_avg    , W_err    ) + '\n'
s = s + '  W_dot = %.3E \pm %.3E'  % (W_dot_avg, W_dot_err) + '\n'
s = s + '  P     = %.3E \pm %.3E'  % (P_avg    , P_err    ) + '\n'
s = s + '  D     = %.3E \pm %.3E'  % (D_avg    , D_err    ) + '\n'
print (s)

f = open('time_average.txt','w') 
f.write(s) 
f.close() 


ys = [ 
       wmag_dot_sum, 
       wmag_dissip_sum, 
       -p_ext_sum, 
       wmag_dot_sum + wmag_dissip_sum - p_ext_sum,
       np.full([tt[avg_start:avg_end].size], W_dot_avg),
       np.full([tt[avg_start:avg_end].size], D_avg),
       np.full([tt[avg_start:avg_end].size], -P_avg),
     ]
xs = [
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
        'k--', 
        '', 
        '', 
        '', 
     ]
legends = [ 
       r'$\mathrm{d}W/\mathrm{d} t$', 
       r'$D_B$', 
       r'$-P$', 
       r'balance', 
       '', 
       '', 
       '', 
     ]
plot_1d_many_average(xs, ys, tt[avg_start], tt[avg_end], xlab='$'+tlab+'$', legends=legends, ls=ls, legendloc='upper left', title='', ylab='', term=True, save=outdir + 'balance_all_avg.pdf')


ys = [ 
       wmag_sum, 
       np.full([tt[avg_start:avg_end].size], wmag_sum_avg),
     ]
ls = [ 
        '', 
        '', 
     ]
xs = [
       tt,
       tt[avg_start:avg_end],
     ]
legends = [ 
       r'$W_{B}$', 
       '',
     ]
plot_1d_many_average(xs, ys, tt[avg_start], tt[avg_end], xlab='$'+tlab+'$', legends=legends, ls=ls, legendloc='upper left', title='', ylab='', term=True, save=outdir + 'energy_all_avg.pdf')

##########################################################
#                    average kspectrum                   #
##########################################################
print('\nplotting kspectrum\n')
outdir = './fig_kspectrum/'

b2_bin_avg  = time_average(tt[avg_start:avg_end], b2_bin [avg_start:avg_end], axis=0)
bx2_bin_avg = time_average(tt[avg_start:avg_end], bx2_bin[avg_start:avg_end], axis=0)
by2_bin_avg = time_average(tt[avg_start:avg_end], by2_bin[avg_start:avg_end], axis=0)
bz2_bin_avg = time_average(tt[avg_start:avg_end], bz2_bin[avg_start:avg_end], axis=0)

# The bins are 3D isotropic shells that already contain the -kz modes
# (Hermitian doubling), so no folding of +-kz is needed here.
# Only the 2/3-dealiased range of the shells is fully populated.
if nlz == nkz:
  kp_end = np.argmin(np.abs(kpbin - kpbin.max()*2./3.))
else:
  kp_end = int(kpbin.size*2./3.)

#--------------------------------------------------------#
#                      plot 1D spectra                   #
#--------------------------------------------------------#
ys = [ 
       b2_bin_avg[1:kp_end], 
       kpbin[1:kp_end]**(-7./3.)/kpbin[1]**(-7./3.)*b2_bin_avg[1:kp_end][0],
     ]
xs = [ 
      kpbin[1:kp_end], 
      kpbin[1:kp_end]  
     ]
ls = [ 
        '', 
        'k-.', 
     ]
legends = [ 
            r'$E_{B}$',
            r'-7/3',
          ]

plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', ylab='', term=True, save=outdir+'k_spectra_avg.pdf')


# B by components
ys = [ 
       bx2_bin_avg[1:kp_end], 
       by2_bin_avg[1:kp_end], 
       bz2_bin_avg[1:kp_end], 
       kpbin[1:kp_end]**(-7./3.)/kpbin[1]**(-7./3.)*b2_bin_avg[1:kp_end][0],
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
        '', 
        'k-.', 
     ]
legends = [ 
            r'$E_{B_x}$', 
            r'$E_{B_y}$', 
            r'$E_{B_z}$', 
            r'-7/3',
          ]

plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', ylab='', term=True, save=outdir+'k_spectra_b_avg.pdf')

#------------------#
#   output ascii   #
#------------------#
np.savetxt(outdir + 'Ek_avg.txt'  , np.column_stack((kpbin[:kp_end], 
                                                  b2_bin_avg[:kp_end],
                                                 bx2_bin_avg[:kp_end],
                                                 by2_bin_avg[:kp_end],
                                                 bz2_bin_avg[:kp_end],
                                                 )), fmt='%E')

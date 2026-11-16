# -*- coding: utf-8 -*-

#####################################################
## main program for making plots from AstroGK data ##
#####################################################
import numpy as np
try:                                    # numpy >= 2.0
  from numpy import trapezoid
except ImportError:                     # numpy < 2.0
  from numpy import trapz as trapezoid

def time_average(x, y, axis=0): # x: 1D array, y: any-D array
  return trapezoid(y, x, axis=axis)/(x[-1] - x[0])

avg_start = 0
avg_end   = -1

##########################################################
#              average energy time evolution             #
##########################################################
print('\nplotting energy\n')
outdir = './fig_energy/'

# load data
# columns of energies.txt: tt, wkin, wmag, wrho, wkin_dot, wmag_dot, wrho_dot,
#                          wkin_dissip, wmag_dissip, wrho_dissip, p_ext, p_re, p_ma
time         = np.transpose(np.loadtxt(outdir + 'energies.txt'    ))[0]
wkin_dot_sum = np.transpose(np.loadtxt(outdir + 'energies.txt'    ))[4]
wmag_dot_sum = np.transpose(np.loadtxt(outdir + 'energies.txt'    ))[5]
wrho_dot_sum = np.transpose(np.loadtxt(outdir + 'energies.txt'    ))[6]

energy_dot = wkin_dot_sum + wmag_dot_sum + wrho_dot_sum

energy_dot_avg = time_average(time[avg_start:avg_end], energy_dot[avg_start:avg_end], axis=0)
print ('energy_dot     = %.3E' % energy_dot_avg + ' ! This must be zero when stationary')

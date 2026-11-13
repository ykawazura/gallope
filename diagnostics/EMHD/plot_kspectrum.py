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

def add_vertical_line(xs, ys, ls, legends):
  xs.append([1./de, 1./de])
  ys.append([b2_bin[final_idx, 1:kp_end].min(), 
             b2_bin[final_idx, 1:kp_end].max()])
  ls.append('k:')
  legends.append(r'$kd_\rme = 1$')

  return xs, ys, ls, legends

#--------------------------------------------------------#
#                      plot 1D spectra                   #
#--------------------------------------------------------#
ys = [ 
       b2_bin[final_idx, 1:kp_end], 
       kpbin[1:kp_end]**(-7./3.)/kpbin[1]**(-7./3.)*b2_bin[final_idx,1:kp_end][0],
     ]
xs = [ 
      kpbin[1:kp_end], 
      kpbin[1:kp_end], 
     ]
ls = [ 
        '', 
        'k--', 
     ]
legends = [ 
            r'$E$',
            r'-7/3',
          ]
if de > 0.0: xs, ys, ls, legends = add_vertical_line(xs, ys, ls, legends)

plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt[final_idx], ylab='', term=True, save=outdir+'k_spectra.pdf')


# B by components
ys = [ 
       bx2_bin[final_idx, 1:kp_end], 
       by2_bin[final_idx, 1:kp_end], 
       bz2_bin[final_idx, 1:kp_end], 
       kpbin[1:kp_end]**(-7./3.)/kpbin[1]**(-7./3.)*b2_bin[final_idx,1:kp_end][0],
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
        '', 
        'k--', 
     ]
legends = [ 
            r'$E_{B_x}$', 
            r'$E_{B_y}$', 
            r'$E_{B_z}$', 
            r'-7/3',
          ]

plot_log1d_many(xs, ys, xlab=kplab, legends=legends, ls=ls, legendloc='lower left', title=r'$t = %.2E $' % tt[final_idx], ylab='', term=True, save=outdir+'k_spectra_b.pdf')

#------------------#
#   output ascii   #
#------------------#
np.savetxt(outdir + 'Ek.txt'  , np.column_stack((kpbin    [:kp_end], 
                                                  b2_bin  [final_idx,:kp_end],
                                                 bx2_bin  [final_idx,:kp_end],
                                                 by2_bin  [final_idx,:kp_end],
                                                 bz2_bin  [final_idx,:kp_end],
                                                 )), fmt='%E')

#--------------------------------------------------------#
#                       plot movie                       #
#--------------------------------------------------------#
if ismovie:
  ys = []
  ys.append(b2_bin[:, 1:kp_end])
  ys.append(np.tile(kpbin[1:kp_end], (nt, 1))**(-7./3.)/kpbin[1]*np.tile(b2_bin[:,1], (kp_end-1, 1)).T)
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
              r'$E_{b}$', 
              r'-7/3',
            ]
  movie_log1d_many(tt, xs, ys, legends=legends, ls=ls, legendloc='lower left', xlab=kplab, save=outdir+'Ek.gif')

# -*- coding: utf-8 -*-
import warnings
warnings.filterwarnings('ignore')
#-------------------------------------------------------------#
#                      Matplotlib setting                     #
#-------------------------------------------------------------#
import pylab
import matplotlib as mpl
# mpl.use('Agg')
from matplotlib import animation
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import *
import matplotlib.ticker as ticker
from parula import parula_map
from scipy.io import loadmat

interpolation        = 'nearest' # or 'nearest'
default_colormap     = parula_map # or parula_map or 'viridis' for sequential colormaps and RdBu for diverging colormap
default_colormap_log = 'inferno' # or parula_map or 'viridis'

# setup some plot defaults
linewidth = 2
fontsize  = 30
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', serif='Times')
plt.rc('font', size=fontsize)
plt.rc('axes', linewidth=linewidth)
plt.rc('axes', labelsize=fontsize)
plt.rc('legend', fontsize=fontsize)
plt.rc('xtick', labelsize=fontsize)
plt.rc('xtick', top=True)
plt.rc('xtick.major', width=linewidth)
plt.rc('xtick.major', size=17)
plt.rc('xtick.minor', width=linewidth)
plt.rc('xtick.minor', visible=True)
plt.rc('xtick.minor', size=8)
plt.rc('ytick', labelsize=fontsize)
plt.rc('ytick', right=True)
plt.rc('ytick.major', width=linewidth)
plt.rc('ytick.major', size=17)
plt.rc('ytick.minor', width=linewidth)
plt.rc('ytick.minor', visible=True)
plt.rc('ytick.minor', size=8)
# plt.rc('xtick', labelbottom='off')
plt.rc('xtick', direction='in')
# plt.rc('ytick', labelleft='off')
plt.rc('ytick', direction='in')
rcParams.update({'figure.autolayout': True})
from latex_preamble import preamble
rcParams['text.latex.preamble'] = preamble  
tab10 = plt.get_cmap('tab10').colors
#-------------------------------------------------------------#

import numpy as np

names   = [r'$\calT_{uu}$', r'$\calT_{BB}$', r'$\calT_{uB}$', r'$\calT_{Bu}$']

# load grid
data = loadmat('../grid')
kp   = data['kpbin_log'][0]

# load u
data = loadmat('../trans')
trans_uu = data['nltrans_uu']
trans_bb = data['nltrans_bb']
trans_ub = data['nltrans_ub']
trans_bu = data['nltrans_bu']

k0idx = 3

transes = [trans_uu, trans_bb, trans_ub, trans_bu]

cuts = [10, 30, 70]

fig, axes = plt.subplots(6, 2, figsize=(15.2, 20))
for i, trans in enumerate(transes):
  col = i%2
  row = int(i/2)

  umax = np.max(abs(trans))

  # 2D map
  ax = plt.subplot2grid((6, 2), (3*row, col), rowspan=2)
  X, Y = np.meshgrid(kp[k0idx:], kp[k0idx:])
  im = ax.pcolormesh(X, Y, trans.T[k0idx:, k0idx:], cmap='RdBu_r', vmin = -umax, vmax = umax)

  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  cb = plt.colorbar(im, cax=cax)
  cb.formatter.set_powerlimits((0, 0))
  cb.ax.yaxis.set_major_locator(MaxNLocator(5,prune='upper'))
  cb.ax.yaxis.set_offset_position('left')  
  cb.locator = ticker.MaxNLocator(nbins=9)
  cb.update_ticks()

  ax.set_title(names[i])

  ax.minorticks_on()

  ax.set_xscale('log')
  ax.set_yscale('log')
  ax.set_xlim(kp[k0idx], kp[-1])
  ax.set_ylim(kp[k0idx], kp[-1])

  ax.tick_params(labelbottom=False)

  ax.set_ylabel(r'$K$')

  # add horizontal cutline
  for cut in cuts:
    idx = np.argmin(np.abs(kp - cut))  
    ax.plot([kp[k0idx], kp[-1]], [kp[idx], kp[idx]], '--')


  # 1D cut
  ax = plt.subplot2grid((6, 2), (3*row+2, col))
  for cut in cuts:
    idx = np.argmin(np.abs(kp - cut))  
    ax.plot(kp[k0idx:], trans[k0idx:, idx], 'o-', lw=2, label=r'$K=%.1f$' % kp[idx])

  ax.minorticks_on()
  ax.set_xscale('log')
  ax.set_xlabel(r'$Q$')
  ax.set_xlim(kp[k0idx], kp[-1])
  ax.set_ylim(-umax, umax)

  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  cb = plt.colorbar(im, cax=cax)
  cb.formatter.set_powerlimits((0, 0))
  cb.ax.yaxis.set_major_locator(MaxNLocator(5,prune='upper'))
  cb.ax.yaxis.set_offset_position('left')  
  cb.locator = ticker.MaxNLocator(nbins=9)
  cb.update_ticks()

plt.tight_layout(pad=0.0, w_pad=0.0, h_pad=0.0)
plt.savefig('trans.pdf')



#--------------------------------------------------------#
#                         Fluxes                         #
#--------------------------------------------------------#
def plot_semilogx1d_many(xs, ys, xlab, legends, ls, xmin='', xmax='', ymin='', ymax='', legendloc='upper left', title='', ylab='', term=True, save=False):
  fig = plt.figure(figsize=(9, 8))

  for i, y in enumerate(ys):
    plt.semilogx(xs[i], ys[i], ls[i], lw=2, label=legends[i])

  legenaxd = plt.legend(fontsize=20, loc=legendloc)
  plt.xlabel(xlab)
  if xmin != '':
    plt.xlim(xmin=xmin)
  if xmax != '':
    plt.xlim(xmax=xmax)
  if ymin != '':
    plt.ylim(ymin=ymin)
  if ymax != '':
    plt.ylim(ymax=ymax)
  if len(ylab) > 0:
    plt.ylabel(ylab)
  if len(title) > 0:
    plt.title(title)

  # force exponential yticks
  ax = plt.gca()
  ax.tick_params(which='both', direction='in')


  if save:
    print ('--   ' + title + '   --')
    plt.savefig(save)

  if term:
    plt.show()
  return fig

flx_uu = np.asarray([np.sum(trans_uu.T[i:, :i]) for i in np.arange(0, kp.size)])
flx_bb = np.asarray([np.sum(trans_bb.T[i:, :i]) for i in np.arange(0, kp.size)])
flx_ub = np.asarray([np.sum(trans_ub.T[i:, :i]) for i in np.arange(0, kp.size)])
flx_bu = np.asarray([np.sum(trans_bu.T[i:, :i]) for i in np.arange(0, kp.size)])

ys = [ 
      flx_uu[:],
      flx_bb[:],
      flx_ub[:],
      flx_bu[:],
      (flx_uu[:] + flx_bb[:] + flx_bu[:] + flx_ub[:]),
     ]
xs = [ 
      kp[:], 
      kp[:], 
      kp[:], 
      kp[:], 
      kp[:], 
     ]
ls = [ 
        '', 
        '', 
        '', 
        '', 
        '', 
     ]
legends = [ 
            r'$\Pi^{u^<}_{u^>}$', 
            r'$\Pi^{B^<}_{B^>}$', 
            r'$\Pi^{u^<}_{B^>}$', 
            r'$\Pi^{B^<}_{u^>}$', 
            r'$\Pi_\mr{tot}$', 
          ]
plot_semilogx1d_many(xs, ys, xlab='k', legends=legends, ls=ls, legendloc='lower left', ylab='', term=False, save='flux.pdf')

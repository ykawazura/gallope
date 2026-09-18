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
from cycler import cycler
plt.rcParams['axes.prop_cycle'] = cycler(
        color=["#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7", "#999999"]) # Okabe-Ito colorscheme
rcParams.update({'figure.autolayout': True})
from latex_preamble import preamble
rcParams['text.latex.preamble'] = preamble  
tab10 = plt.get_cmap('tab10').colors
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
#-------------------------------------------------------------#

grids = ['512']

import numpy as np
import os

fig, ax = plt.subplots(1, 1, figsize=(10, 9))

for igrid, grid in enumerate(grids):

    for lib in ['cuFFTMp-slab', 'cuFFTMp-pencil']:
        dir = lib + '/'+grid+'^3/'

        nodes = []
        time  = []

        sorted_files = sorted(os.listdir(dir), key=lambda x: int(x.replace("node", "")))
        for file in sorted_files:
            node = int(file.split('node')[1])
            nodes.append(node)

            lines = [line.rstrip('\n') for line in open(dir+file)]
            line = [x for x in lines if 'cpu time per loop =' in x]

            time.append(float(line[0].split('=')[2]))

        if lib == 'cuFFTMp-slab':
            marker = 'o'; mfc=colors[igrid]; ls = '-' ; lw = 4
        if lib == 'cuFFTMp-pencil':
            marker = 's'; mfc='none'; ls = '--'; lw = 2

        # if lib == 'cuFFTMp' and grid == '2048':
        #     plt.loglog([nodes[0], nodes[-1]], [time[0], time[0]*nodes[0]/nodes[-1]], 'k-')


        # ax.semilogx(nodes[1:], time[1]/np.asarray(time[1:]), lw=lw, marker=marker, ls=ls , ms=15, mew=4, mfc=mfc, color=colors[igrid], label=grid + ' (' + lib + ')')
        ax.semilogx(nodes[1:], time[1:], lw=lw, marker=marker, ls=ls , ms=15, mew=4, mfc=mfc, color=colors[igrid], label=grid + ' (' + lib + ')')

leg = ax.legend(frameon=False, loc='upper right', fontsize=18, handlelength=1.1, columnspacing=0.6)


ax.set_xlim([1, 400])
ax.set_ylim([0, 1])
ax.set_xlabel(r'Number of Nodes')
ax.set_ylabel(r'Time per FFT + iFFT [s]')

plt.savefig('slab_vs_pencil.pdf')

###############################################################################################################

grids = ['512', '1024']

fig, ax = plt.subplots(1, 1, figsize=(10, 9))

for igrid, grid in enumerate(grids):

    for lib in ['cuFFTMp-slab']:
        dir = lib + '/'+grid+'^3/'

        nodes = []
        time  = []

        sorted_files = sorted(os.listdir(dir), key=lambda x: int(x.replace("node", "")))
        for file in sorted_files:
            node = int(file.split('node')[1])
            nodes.append(node)

            lines = [line.rstrip('\n') for line in open(dir+file)]
            line = [x for x in lines if 'cpu time per loop =' in x]

            time.append(float(line[0].split('=')[2]))

        if lib == 'cuFFTMp-slab':
            marker = 'o'; mfc=colors[igrid]; ls = '-' ; lw = 4

        # if lib == 'cuFFTMp' and grid == '2048':
        #     plt.loglog([nodes[0], nodes[-1]], [time[0], time[0]*nodes[0]/nodes[-1]], 'k-')


        ax.semilogx(nodes[1:], time[1]/np.asarray(time[1:]), lw=lw, marker=marker, ls=ls , ms=15, mew=4, mfc=mfc, color=colors[igrid], label=grid + ' (' + lib + ')')

leg = ax.legend(frameon=False, loc='upper right', fontsize=18, handlelength=1.1, columnspacing=0.6)


ax.set_xlim([1, 400])
ax.set_ylim(ymin=0)
ax.set_xlabel(r'Number of Nodes')
ax.set_ylabel(r'Time per FFT + iFFT [s]')

plt.savefig('plot.pdf')

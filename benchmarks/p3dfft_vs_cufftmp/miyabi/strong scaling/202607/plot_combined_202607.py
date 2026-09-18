# -*- coding: utf-8 -*-
# Fig. 1 (2026-07 campaign), 2x2 figure: left column N=1024, right column N=2048.
#   top row    -- absolute time per FFT + iFFT (the originally submitted panels)
#   bottom row -- GPU-to-CPU speedup, i.e. the ratio of the P3DFFT time to the
#                 cuFFTMp time at the same node count
# The ratio row exists because the absolute panels answer "how fast" but not
# "how much faster than the CPU baseline", which is the quantity the text quotes
# and the reader has to reconstruct by eye from two log curves otherwise.
# Median over trials with min-max whiskers (referee R1-4 variability).
# Output: cuFFTMp-scaling.pdf (drop-in for the paper).
import warnings
warnings.filterwarnings('ignore')
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import *
from parula import parula_map

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
plt.rc('xtick', direction='in')
plt.rc('ytick', direction='in')
from cycler import cycler
plt.rcParams['axes.prop_cycle'] = cycler(
        color=["#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7", "#999999"])  # Okabe-Ito
mpl.rcParams['figure.autolayout'] = False   # we place the shared labels manually
from latex_preamble import preamble
rcParams['text.latex.preamble'] = preamble
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

import os
import re
from collections import defaultdict
from statistics import median

SUFFIX = '-202607'
grids  = ['1024', '2048']
LIBS   = ['cuFFTMp-slab', 'cuFFTMp-pencil', 'P3DFFT']
MARKER = {'cuFFTMp-slab': 'o', 'cuFFTMp-pencil': '^', 'P3DFFT': 's'}


def read_run(path):
    txt = open(path, errors='ignore').read()
    m = re.search(r'nx\s*=\s*(\d+)', txt)
    nx = int(m.group(1)) if m else None
    hit = [x for x in txt.splitlines() if 'cpu time per loop =' in x]
    if not hit:
        raise ValueError('no loop-time marker')
    return nx, float(hit[0].split('=')[2])


def collect(d, grid):
    out = defaultdict(list)
    for f in os.listdir(d):
        m = re.fullmatch(r'node(\d+)(?:-(\d+))?', f)
        if not m:
            continue
        try:
            nx, t = read_run(os.path.join(d, f))
        except (IndexError, ValueError):
            continue
        if nx is not None and nx != int(grid):
            continue
        out[int(m.group(1))].append(t)
    return out


fig = plt.figure(figsize=(16, 12.0))
# The ratio row is the shorter one: it carries a single dimensionless number per
# point, whereas the absolute row has to resolve three curves over three decades.
gs  = fig.add_gridspec(2, 2, height_ratios=[1.5, 0.8],
                       left=0.105, right=0.99, bottom=0.115, top=0.99,
                       hspace=0.05, wspace=0.04)
# absolute row first, so the ratio row below can inherit its x-axis
axes  = [fig.add_subplot(gs[0, 0])]
axes += [fig.add_subplot(gs[0, 1], sharey=axes[0])]
axr   = [fig.add_subplot(gs[1, 0], sharex=axes[0])]
axr  += [fig.add_subplot(gs[1, 1], sharex=axes[1], sharey=axr[0])]

for ax, axrat, grid in zip(axes, axr, grids):
    med_of = {}          # lib -> {node: median time}, for the ratio panel
    for ilib, lib in enumerate(LIBS):
        data = collect(f'{lib}{SUFFIX}/{grid}^3', grid)
        if lib == 'cuFFTMp-pencil':
            data = {n: v for n, v in data.items()
                    if int(n ** 0.5 + 0.5) ** 2 == n}
        nodes = sorted(data)
        # The protocol claimed in the paper is three independent submissions per
        # point.  read_run() already drops files whose header does not match the
        # grid, so a launch failure leaves a name on disk but no sample -- and a
        # point quietly plotted from one or two runs would still look like every
        # other point.  Say the count out loud for every point instead, so that
        # any shortfall has to be either fixed or disclosed, never overlooked.
        short = [(n, len(data[n])) for n in nodes if len(data[n]) != 3]
        if short:
            print(f'  !! {lib} N={grid}: trials != 3 at ' +
                  ', '.join(f'node{n}(n={c})' for n, c in short))
        med = [median(data[n]) for n in nodes]
        lo  = [median(data[n]) - min(data[n]) for n in nodes]
        hi  = [max(data[n]) - median(data[n]) for n in nodes]
        med_of[lib] = {n: (median(data[n]), min(data[n]), max(data[n]))
                       for n in nodes}

        if lib == 'cuFFTMp-slab':   # ideal strong-scaling guide
            ax.loglog([nodes[1], nodes[-1]],
                      [med[1], med[1] * nodes[1] / nodes[-1]], 'k-')

        # capsize is a half-length: matplotlib draws the cap as a '_' marker of
        # size 2*capsize, so capsize=12 gives a 24 pt bar against a 15 pt marker
        # and the cap stays visible when the spread is smaller than the symbol.
        ax.errorbar(nodes, med, yerr=[lo, hi], lw=3, ms=15, mew=3,
                    marker=MARKER[lib], ls='-', mfc=colors[ilib],
                    color=colors[ilib], capsize=12, capthick=linewidth,
                    label=lib)

    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlim([0.9, 400])
    ax.set_ylim([1e-3, 2.9])
    leg = ax.legend(frameon=False, loc='lower left', fontsize=22,
                    columnspacing=0.6, title=rf'$N = {grid}$')
    leg.get_title().set_fontsize(22)
    plt.setp(ax.get_xticklabels(), visible=False)

    # ---- ratio panel: how many times faster than the CPU baseline ----------
    # The whiskers are the widest ratio the samples permit -- slowest CPU over
    # fastest GPU and vice versa -- rather than a propagated error, because with
    # three samples the extremes are the honest statement of what was observed.
    cpu = med_of['P3DFFT']
    for ilib, lib in enumerate(['cuFFTMp-slab', 'cuFFTMp-pencil']):
        gpu = med_of[lib]
        ns  = [n for n in sorted(gpu) if n in cpu]
        if not ns:
            continue
        r   = [cpu[n][0] / gpu[n][0] for n in ns]
        rlo = [r[i] - cpu[n][1] / gpu[n][2] for i, n in enumerate(ns)]
        rhi = [cpu[n][2] / gpu[n][1] - r[i] for i, n in enumerate(ns)]
        print(f'  ratio P3DFFT/{lib} N={grid}: ' +
              ', '.join(f'node{n}={v:.3g}' for n, v in zip(ns, r)))
        # Pencil is computed and printed but not drawn: the solver uses slab,
        # and the pencil ratio hugs parity closely enough that plotting it
        # forces the ordinate below 1 and squashes the slab curve, which is the
        # one the text quotes.  The pencil-versus-slab comparison stays in the
        # top row, where it belongs.
        if lib != 'cuFFTMp-slab':
            continue
        axrat.errorbar(ns, r, yerr=[rlo, rhi], lw=3, ms=15, mew=3,
                       marker=MARKER[lib], ls='-', mfc=colors[ilib],
                       color=colors[ilib], capsize=12, capthick=linewidth,
                       label=rf'P3DFFT / {lib}')
    # No parity line: the ordinate now starts at 1, so the bottom spine is
    # parity and a dotted line there would be drawn underneath it.
    axrat.set_xscale('log'); axrat.set_yscale('log')
    axrat.set_ylim([1.0, 100.0])

axr[0].set_ylabel(r'Speedup over P3DFFT')
# The legend goes in the right-hand ratio panel: N=2048 has no data below eight
# nodes, so its left half is empty, whereas the N=1024 panel starts at its
# highest point and would collide with a legend anywhere along the top.
axr[1].legend(frameon=False, loc='upper left', fontsize=22, columnspacing=0.6)
plt.setp(axr[1].get_yticklabels(), visible=False)
plt.setp(axes[1].get_yticklabels(), visible=False)
axes[0].set_ylabel(r'Time per FFT + iFFT [s]')
# The two rows carry tick labels of different widths ('10^{-3}' against
# '10^{0}'), so matplotlib places their y-labels at different x by default and
# the pair reads as misaligned.  align_ylabels moves both to the leftmost of
# the two positions.
fig.align_ylabels([axes[0], axr[0]])
fig.text(0.55, 0.028, r'Number of nodes (= number of CPUs = number of GPUs)',
         ha='center', va='bottom', fontsize=fontsize)
plt.savefig('cuFFTMp-scaling.pdf')
print('wrote cuFFTMp-scaling.pdf')

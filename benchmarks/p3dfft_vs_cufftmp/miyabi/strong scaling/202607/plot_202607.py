# -*- coding: utf-8 -*-
# Fig. 1 (2026-07 campaign).  Reads nodeN-t files (t = trial index) and plots
# the median over trials with min-max whiskers, so the referee's request for
# run-to-run variability (R1-4) is answered by the figure itself.
import warnings
warnings.filterwarnings('ignore')
#-------------------------------------------------------------#
#                      Matplotlib setting                     #
#-------------------------------------------------------------#
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
rcParams.update({'figure.autolayout': True})
from latex_preamble import preamble
rcParams['text.latex.preamble'] = preamble
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
#-------------------------------------------------------------#

import os
import re
from collections import defaultdict
from statistics import median

SUFFIX = '-202607'
grids  = ['1024', '2048']
LIBS   = ['cuFFTMp-slab', 'cuFFTMp-pencil', 'P3DFFT']
MARKER = {'cuFFTMp-slab': 'o', 'cuFFTMp-pencil': '^', 'P3DFFT': 's'}


def read_run(path):
    """(nx, seconds/loop) from a micro-benchmark stdout.

    2048^3 does not fit on 1-2 GPUs, so those sub-runs silently fell back to
    1024^3; a few pencil sub-runs were mis-mapped the other way.  We therefore
    read the grid size actually run (the ``nx =`` header) and let collect()
    drop any file whose grid does not match the directory it sits in.
    """
    txt = open(path, errors='ignore').read()
    m = re.search(r'nx\s*=\s*(\d+)', txt)
    nx = int(m.group(1)) if m else None
    hit = [x for x in txt.splitlines() if 'cpu time per loop =' in x]
    if not hit:
        raise ValueError('no loop-time marker')
    return nx, float(hit[0].split('=')[2])


def collect(d, grid):
    """{node count: [time, ...]} over trials in d, keeping only grid-correct runs."""
    out = defaultdict(list)
    for f in os.listdir(d):
        m = re.fullmatch(r'node(\d+)(?:-(\d+))?', f)
        if not m:
            continue
        try:
            nx, t = read_run(os.path.join(d, f))
        except (IndexError, ValueError):
            print(f'  !! unparsable, skipped: {d}/{f}')
            continue
        if nx is not None and nx != int(grid):
            print(f'  !! wrong grid nx={nx} (want {grid}), skipped: {d}/{f}')
            continue
        out[int(m.group(1))].append(t)
    return out


for grid in grids:
    fig, ax = plt.subplots(1, 1, figsize=(9.5, 8.5))

    for ilib, lib in enumerate(LIBS):
        data = collect(f'{lib}{SUFFIX}/{grid}^3', grid)
        # cuFFTMp pencil uses a 2D process grid, which is only balanced at
        # perfect-square node counts; non-square counts (2,8,32,128) give a
        # lopsided grid and a misleading zig-zag.  Plot squares only, as in the
        # originally submitted figure.
        if lib == 'cuFFTMp-pencil':
            data = {n: v for n, v in data.items()
                    if int(n ** 0.5 + 0.5) ** 2 == n}
        # P3DFFT is reported up to 128 nodes, as in the originally submitted figure.
        if lib == 'P3DFFT':
            data = {n: v for n, v in data.items() if n <= 128}
        nodes = sorted(data)
        med = [median(data[n]) for n in nodes]
        lo  = [median(data[n]) - min(data[n]) for n in nodes]
        hi  = [max(data[n]) - median(data[n]) for n in nodes]

        for n in nodes:
            v = data[n]
            print(f'{grid:5s} {lib:16s} node{n:<4d} n={len(v)}  median={median(v):.5f}  '
                  f'min={min(v):.5f}  max={max(v):.5f}  spread={100*(max(v)-min(v))/median(v):.1f}%')

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

    leg = ax.legend(frameon=False, loc='lower left', fontsize=22,
                    columnspacing=0.6, title=rf'$N = {grid}$')
    leg.get_title().set_fontsize(22)

    ax.set_xlim([0.9, 400])
    ax.set_ylim([1e-3, 2.9])
    ax.set_xlabel(r'Number of nodes (= number of CPUs = number of GPUs)')
    ax.set_ylabel(r'Time per FFT + iFFT [s]')
    plt.savefig(f'{grid}{SUFFIX}.pdf')

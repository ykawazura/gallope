# -*- coding: utf-8 -*-
# Fig. 3 (2026-07 campaign).  Weak scaling of the cuFFTMp slab distributed FFT:
# the work per GPU is held fixed while the global grid grows with the node count.
# Reads nodeN-t files (t = trial index) and plots the median over trials with
# min-max whiskers, matching the strong-scaling figures (R1-4 variability).
#
# The ideal-weak-scaling guide is referenced to *2* nodes, not 1: on a single
# node the transform is entirely GPU-local and the all-to-all that a distributed
# FFT requires is simply absent, so node 1 is not a meaningful weak-scaling
# baseline (same effect discussed for the strong scaling in Comment 2.3).
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
BASES  = ['512', '1024']          # per-GPU grid held fixed along each series
MARKER = {'512': 'o', '1024': 's'}
REF    = 2                        # node count the ideal guide is referenced to


def read_run(path):
    """(points, seconds/loop) from a micro-benchmark stdout.

    The grid actually run is read from the header rather than inferred from the
    file name, so that a sub-run which silently fell back to another size cannot
    enter the figure (lesson from the strong-scaling campaign).
    """
    txt = open(path, errors='ignore').read()
    m = re.search(r'nx\s*=\s*(\d+),\s*ny\s*=\s*(\d+),\s*nz\s*=\s*(\d+),\s*nproc\s*=\s*(\d+)', txt)
    if not m:
        raise ValueError('no grid header')
    nx, ny, nz, nproc = (int(g) for g in m.groups())
    hit = [x for x in txt.splitlines() if 'cpu time per loop =' in x]
    if not hit:
        raise ValueError('no loop-time marker')
    return nx * ny * nz, nproc, float(hit[0].split('=')[2])


def collect(d, base):
    """{node count: [time, ...]}, keeping only runs whose work per GPU is correct."""
    want = int(base) ** 3
    out = defaultdict(list)
    for f in os.listdir(d):
        m = re.fullmatch(r'node(\d+)(?:-(\d+))?', f)
        if not m:
            continue
        node = int(m.group(1))
        try:
            pts, nproc, t = read_run(os.path.join(d, f))
        except (IndexError, ValueError) as e:
            print(f'  !! unparsable ({e}), skipped: {d}/{f}')
            continue
        if nproc != node:
            print(f'  !! nproc={nproc} != node{node}, skipped: {d}/{f}')
            continue
        if pts != want * node:
            print(f'  !! work/GPU {pts // node:,} != {want:,}, skipped: {d}/{f}')
            continue
        out[node].append(t)
    return out


fig, ax = plt.subplots(1, 1, figsize=(9.5, 8.5))

for ibase, base in enumerate(BASES):
    data = collect(f'cuFFTMp-slab{SUFFIX}/{base}^3', base)
    nodes = sorted(data)
    med = [median(data[n]) for n in nodes]
    lo  = [median(data[n]) - min(data[n]) for n in nodes]
    hi  = [max(data[n]) - median(data[n]) for n in nodes]

    for n in nodes:
        v = data[n]
        print(f'{base:>4s}^3/GPU  node{n:<4d} n={len(v)}  median={median(v):.5f}  '
              f'min={min(v):.5f}  max={max(v):.5f}  '
              f'spread={100*(max(v)-min(v))/median(v):.1f}%')

    ref = median(data[REF])
    growth = median(data[nodes[-1]]) / ref
    print(f'{base:>4s}^3/GPU  {REF} -> {nodes[-1]} nodes: x{growth:.2f}  '
          f'({REF} -> 128 nodes: x{median(data[128]) / ref:.2f})\n')

    # ideal weak scaling = constant time, referenced to REF nodes
    ax.loglog([REF, nodes[-1]], [ref, ref], '--', color=colors[ibase],
              lw=linewidth, alpha=0.9, dashes=(6, 4))
    # capsize is a half-length: matplotlib draws the cap as a '_' marker of
    # size 2*capsize, so capsize=12 gives a 24 pt bar against a 15 pt marker
    # and the cap stays visible when the spread is smaller than the symbol.
    ax.errorbar(nodes, med, yerr=[lo, hi], lw=3, ms=15, mew=3,
                marker=MARKER[base], ls='-', mfc=colors[ibase],
                color=colors[ibase], capsize=12, capthick=linewidth,
                label=rf'${base}^3$ per GPU')

leg = ax.legend(frameon=False, loc='lower right', fontsize=22,
                columnspacing=0.6)

ax.set_xlim([0.9, 400])
ax.set_ylim([4e-3, 3.0])
ax.set_xlabel(r'Number of nodes (= number of GPUs)')
ax.set_ylabel(r'Time per FFT + iFFT [s]')
plt.savefig(f'weak-scaling{SUFFIX}.pdf')
print('wrote weak-scaling%s.pdf' % SUFFIX)

# -*- coding: utf-8 -*-
# Fig. 3, TSUBAME4.0 version (2026-07 campaign).  Weak scaling of the cuFFTMp
# slab FFT: the work per GPU is held fixed while the global grid grows with the
# GPU count.
#
# On Miyabi-G the abscissa could be called "nodes = GPUs" and the ideal guide had
# to be referenced to 2 nodes, because at 1 node the transform is GPU-local and
# the all-to-all a distributed FFT needs is simply absent.  Here that same
# distinction is visible along the axis instead of hidden at its left end: 1-4
# GPUs are one node, so the all-to-all stays on NVLink, and only x > 4 forces it
# onto InfiniBand.  The guide is therefore referenced to 4 GPUs -- the largest
# allocation that is still a single node -- and the step between 4 and 8 GPUs
# measures the cost of leaving the node with the per-GPU work unchanged.
# The shaded band is that single-node region.
#
# Median over trials with min-max whiskers (referee R1-4 variability).
# Output: weak-scaling-tsubame.pdf
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
rcParams.update({'figure.autolayout': False})   # the top axis needs the margin
from latex_preamble import preamble
rcParams['text.latex.preamble'] = preamble
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
#-------------------------------------------------------------#

import os
import re
from collections import defaultdict
from statistics import median

ROOT   = 'cufftmp-slab/results-202607/tsubame/weak_scaling'
BASES  = ['512', '1024']          # per-GPU grid held fixed along each series
MARKER = {'512': 'o', '1024': 's'}
GPN    = 4                        # H100 per node
INTRA  = 4                        # x <= INTRA is one node: NVLink only
REF    = 4                        # GPUs the ideal guide is referenced to


def read_run(path):
    """(points, nproc, seconds/loop) from a micro-benchmark stdout.

    The grid actually run is read from the header rather than inferred from the
    file name, so a sub-run that silently fell back to another size cannot enter
    the figure (lesson L18).
    """
    txt = open(path, errors='ignore').read().replace('\0', '')
    m = re.search(r'nx\s*=\s*(\d+),\s*ny\s*=\s*(\d+),\s*nz\s*=\s*(\d+),'
                  r'\s*nproc\s*=\s*(\d+)', txt)
    if not m:
        raise ValueError('no grid header')
    nx, ny, nz, nproc = (int(g) for g in m.groups())
    hit = [x for x in txt.splitlines() if 'cpu time per loop =' in x]
    if not hit:
        raise ValueError('no loop-time marker')
    return nx * ny * nz, nproc, float(hit[0].split('=')[2])


def collect(base):
    """{GPUs: [time, ...]}, keeping only runs whose work per GPU is correct."""
    d = os.path.join(ROOT, f'{base}^3')
    want = int(base) ** 3
    out = defaultdict(list)
    if not os.path.isdir(d):
        return out
    for f in os.listdir(d):
        m = re.fullmatch(r'gpu(\d+)-node(\d+)-(\d+)', f)
        if not m:
            continue
        x, nd = int(m.group(1)), int(m.group(2))
        # Packed placement only, same rule as the strong-scaling figure: a run
        # spread thinner than 4 GPUs per node is a placement control, not a
        # point on this curve.
        if nd != -(-x // GPN):
            continue
        try:
            pts, nproc, t = read_run(os.path.join(d, f))
        except (IndexError, ValueError) as e:
            print(f'  !! unparsable ({e}), skipped: {d}/{f}')
            continue
        if nproc != x:
            print(f'  !! nproc={nproc} != gpu{x}, skipped: {d}/{f}')
            continue
        if pts != want * x:
            print(f'  !! work/GPU {pts // x:,} != {want:,}, skipped: {d}/{f}')
            continue
        out[x].append(t)
    return out


fig = plt.figure(figsize=(9.5, 8.8))
# The top spine carries both ticks and a label, so the axes has to stop well
# short of the canvas edge or the node label is clipped.
ax  = fig.add_axes([0.175, 0.125, 0.81, 0.755])

for ibase, base in enumerate(BASES):
    data = collect(base)
    xs = sorted(data)
    if not xs:
        print(f'  -- {base}^3/GPU: no data')
        continue
    short = [(n, len(data[n])) for n in xs if len(data[n]) != 3]
    if short:
        print(f'  !! {base}^3/GPU: trials != 3 at ' +
              ', '.join(f'gpu{n}(n={c})' for n, c in short))
    med = [median(data[n]) for n in xs]
    lo  = [median(data[n]) - min(data[n]) for n in xs]
    hi  = [max(data[n]) - median(data[n]) for n in xs]

    for n in xs:
        v = data[n]
        print(f'{base:>4s}^3/GPU  gpu{n:<4d} n={len(v)}  median={median(v):.5f}  '
              f'min={min(v):.5f}  max={max(v):.5f}  '
              f'spread={100*(max(v)-min(v))/median(v):.1f}%')
    if REF in data:
        ref = median(data[REF])
        print(f'{base:>4s}^3/GPU  {REF} -> {xs[-1]} GPUs: '
              f'x{median(data[xs[-1]]) / ref:.2f}\n')
        # ideal weak scaling = constant time, referenced to REF GPUs
        ax.loglog([REF, xs[-1]], [ref, ref], '--', color=colors[ibase],
                  lw=linewidth, alpha=0.9, dashes=(6, 4))

    # capsize is a half-length: matplotlib draws the cap as a '_' marker of size
    # 2*capsize, so capsize=12 gives a 24 pt bar against a 15 pt marker and the
    # cap stays visible when the spread is smaller than the symbol.
    ax.errorbar(xs, med, yerr=[lo, hi], lw=3, ms=15, mew=3,
                marker=MARKER[base], ls='-', mfc=colors[ibase],
                color=colors[ibase], capsize=12, capthick=linewidth,
                label=rf'${base}^3$ per GPU')

ax.axvspan(0.5, INTRA, color='0.85', zorder=0, lw=0)
ax.legend(frameon=False, loc='lower right', fontsize=22, columnspacing=0.6)
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlim([0.8, 180])
ax.set_ylim([4e-3, 3.0])
ax.set_xlabel(r'Number of GPUs')
ax.set_ylabel(r'Time per FFT + iFFT [s]')

sec = ax.secondary_xaxis('top', functions=(lambda g: g / GPN,
                                           lambda n: n * GPN))
sec.set_xscale('log')
sec.xaxis.set_major_locator(FixedLocator([1, 2, 4, 8, 16, 32]))
sec.xaxis.set_major_formatter(FixedFormatter(['1', '2', '4', '8', '16', '32']))
sec.xaxis.set_minor_formatter(NullFormatter())
sec.tick_params(direction='in', width=linewidth, length=17, labelsize=fontsize)
sec.set_xlabel(r'Number of nodes', labelpad=10)

# Last: set_xscale() resets the locators, so these have to come after it.
ax.xaxis.set_major_locator(FixedLocator([1, 4, 16, 64]))
ax.xaxis.set_major_formatter(FixedFormatter(['1', '4', '16', '64']))
ax.xaxis.set_minor_locator(FixedLocator([2, 8, 32, 128]))
ax.xaxis.set_minor_formatter(NullFormatter())

plt.savefig('weak-scaling-tsubame.pdf')
print('wrote weak-scaling-tsubame.pdf')

# -*- coding: utf-8 -*-
# Weak scaling of the cuFFTMp slab distributed FFT on both machines, as one
# two-panel figure: the work per GPU is held fixed while the global grid grows
# with the allocation.  This replaces the two separate figures produced by
#   miyabi/weak scaling/202607/plot_202607.py   -> weak-scaling-202607.pdf
#   tsubame/plot_weak_tsubame.py                -> weak-scaling-tsubame.pdf
# with no change to how the data are read or reduced; only the layout differs.
#
# (a) Miyabi-G.  One GPU per node, so "nodes" and "GPUs" are the same axis.  The
#     ideal guide is referenced to *2* nodes, not 1: on a single node the
#     transform is entirely GPU-local and the all-to-all that a distributed FFT
#     requires is simply absent, so node 1 is not a weak-scaling baseline.
# (b) TSUBAME 4.0.  Four H100 per node, so that same distinction is visible
#     along the axis instead of hidden at its left end: 1-4 GPUs are one node
#     and the all-to-all stays on NVLink; only x > 4 forces it onto InfiniBand.
#     The guide is therefore referenced to 4 GPUs -- the largest allocation that
#     is still a single node -- and the shaded band is that single-node region.
#
# The two ordinates are identical and shared, which is the point of putting the
# panels side by side: the cross-machine comparison of Section 4.4 can be read
# off directly.  Both panels use the same two series, so the legend is drawn
# once, in (a).
#
# Median over trials with min-max whiskers (referee R1-4 variability).
# Output: weak-scaling-combined.pdf
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

MIYABI  = '../miyabi/weak scaling/202607/cuFFTMp-slab-202607'
TSUBAME = '../tsubame/cufftmp-slab/results-202607/tsubame/weak_scaling'
BASES   = ['512', '1024']         # per-GPU grid held fixed along each series
MARKER  = {'512': 'o', '1024': 's'}
GPN     = 4                       # H100 per TSUBAME node
REF_M   = 2                       # nodes the Miyabi ideal guide is referenced to
REF_T   = 4                       # GPUs the TSUBAME ideal guide is referenced to
YLIM    = [4e-3, 3.0]

# The TSUBAME campaign reaches 32 nodes, but this panel stops at 16, so that
# every TSUBAME figure in the paper covers the same allocations -- the
# strong-scaling figures stop there because the solver and the transform both
# turn over at 32 nodes, and a weak-scaling panel that ran further would leave
# the section quoting three different reaches for one machine.  The 32-node
# runs are on disk and are skipped here rather than deleted, and the skip is
# printed, so that a hole in the data set is never mistaken for coverage.
XMAX_T  = 64                      # GPUs; the TSUBAME panel stops at 16 nodes


def read_run(path):
    """(points, nproc, seconds/loop) from a micro-benchmark stdout.

    The grid actually run is read from the header rather than inferred from the
    file name, so a sub-run that silently fell back to another size cannot enter
    the figure (lesson from the strong-scaling campaign).
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


def collect_miyabi(base):
    """{node count: [time, ...]}, keeping only runs whose work per GPU is right."""
    d = os.path.join(MIYABI, f'{base}^3')
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


def collect_tsubame(base):
    """{GPUs: [time, ...]}, same checks plus the packed-placement rule."""
    d = os.path.join(TSUBAME, f'{base}^3')
    want = int(base) ** 3
    out = defaultdict(list)
    if not os.path.isdir(d):
        return out
    for f in os.listdir(d):
        m = re.fullmatch(r'gpu(\d+)-node(\d+)-(\d+)', f)
        if not m:
            continue
        x, nd = int(m.group(1)), int(m.group(2))
        if x > XMAX_T:
            print(f'  -- beyond {XMAX_T} GPUs, see XMAX_T in the header: {d}/{f}')
            continue
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


def draw(ax, collect, ref, unit, tag):
    """Both panels are the same plot of the same two series; only the abscissa
    and the reference allocation of the ideal guide differ."""
    print(f'--- {tag} ---')
    for ibase, base in enumerate(BASES):
        data = collect(base)
        xs = sorted(data)
        if not xs:
            print(f'  -- {base}^3/GPU: no data')
            continue
        short = [(n, len(data[n])) for n in xs if len(data[n]) != 3]
        if short:
            print(f'  !! {base}^3/GPU: trials != 3 at ' +
                  ', '.join(f'{unit}{n}(n={c})' for n, c in short))
        med = [median(data[n]) for n in xs]
        lo  = [median(data[n]) - min(data[n]) for n in xs]
        hi  = [max(data[n]) - median(data[n]) for n in xs]

        for n in xs:
            v = data[n]
            print(f'{base:>4s}^3/GPU  {unit}{n:<4d} n={len(v)}  '
                  f'median={median(v):.5f}  min={min(v):.5f}  max={max(v):.5f}  '
                  f'spread={100*(max(v)-min(v))/median(v):.1f}%')
        if ref in data:
            r = median(data[ref])
            print(f'{base:>4s}^3/GPU  {ref} -> {xs[-1]} {unit}s: '
                  f'x{median(data[xs[-1]]) / r:.2f}\n')
            # ideal weak scaling = constant time, referenced to `ref`
            ax.loglog([ref, xs[-1]], [r, r], '--', color=colors[ibase],
                      lw=linewidth, alpha=0.9, dashes=(6, 4))

        # capsize is a half-length: matplotlib draws the cap as a '_' marker of
        # size 2*capsize, so capsize=12 gives a 24 pt bar against a 15 pt marker
        # and the cap stays visible when the spread is smaller than the symbol.
        ax.errorbar(xs, med, yerr=[lo, hi], lw=3, ms=15, mew=3,
                    marker=MARKER[base], ls='-', mfc=colors[ibase],
                    color=colors[ibase], capsize=12, capthick=linewidth,
                    label=rf'${base}^3$ per GPU')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_ylim(YLIM)


fig = plt.figure(figsize=(19, 8.8))
# add_axes rather than subplots: the right panel carries a labelled top axis, so
# the axes have to stop well short of the canvas edge or the label is clipped.
axL = fig.add_axes([0.070, 0.115, 0.415, 0.755])
axR = fig.add_axes([0.545, 0.115, 0.415, 0.755])

draw(axL, collect_miyabi, REF_M, 'node', 'Miyabi-G')
draw(axR, collect_tsubame, REF_T, 'gpu', 'TSUBAME 4.0')

# (a) Miyabi-G
axL.set_xlim([0.9, 400])
axL.set_xlabel(r'Number of nodes (= number of GPUs)')
axL.set_ylabel(r'Time per FFT + iFFT [s]')
# One legend for the whole figure: the two panels show the same two series.
axL.legend(frameon=False, loc='lower right', fontsize=22, columnspacing=0.6)
# Powers of four rather than the default powers of ten, so that the two panels
# are read the same way even though their abscissae are different quantities.
axL.xaxis.set_major_locator(FixedLocator([1, 4, 16, 64, 256]))
axL.xaxis.set_major_formatter(FixedFormatter(['1', '4', '16', '64', '256']))
axL.xaxis.set_minor_locator(FixedLocator([2, 8, 32, 128]))
axL.xaxis.set_minor_formatter(NullFormatter())

# (b) TSUBAME 4.0
axR.set_xlim([0.8, 90])
axR.set_xlabel(r'Number of GPUs')
axR.axvspan(0.5, REF_T, color='0.85', zorder=0, lw=0)
axR.tick_params(labelleft=False)        # the ordinate is shared with (a)

sec = axR.secondary_xaxis('top', functions=(lambda g: g / GPN,
                                            lambda n: n * GPN))
sec.set_xscale('log')
sec.xaxis.set_major_locator(FixedLocator([1, 2, 4, 8, 16]))
sec.xaxis.set_major_formatter(FixedFormatter(['1', '2', '4', '8', '16']))
sec.xaxis.set_minor_formatter(NullFormatter())
sec.tick_params(direction='in', width=linewidth, length=17, labelsize=fontsize)
sec.set_xlabel(r'Number of nodes', labelpad=10)

# Last: set_xscale() resets the locators, so these have to come after it.
axR.xaxis.set_major_locator(FixedLocator([1, 4, 16, 64]))
axR.xaxis.set_major_formatter(FixedFormatter(['1', '4', '16', '64']))
axR.xaxis.set_minor_locator(FixedLocator([2, 8, 32]))
axR.xaxis.set_minor_formatter(NullFormatter())

# Panel tags go inside the axes at the upper left, which is empty in both.
for ax, tag in ((axL, r'(a) Miyabi-G'), (axR, r'(b) TSUBAME 4.0')):
    ax.text(0.04, 0.955, tag, transform=ax.transAxes, va='top', ha='left',
            fontsize=fontsize)

plt.savefig('weak-scaling-combined.pdf')
print('wrote weak-scaling-combined.pdf')

# -*- coding: utf-8 -*-
# Fig. 1, TSUBAME4.0 version (2026-07 campaign).  Same 2x2 layout as the Miyabi
# figure -- top row absolute time per FFT + iFFT, bottom row the speedup over
# P3DFFT -- but the x-axis means something different here and that is the point.
#
# On Miyabi-G a node holds one GPU, so "nodes", "GPUs" and "CPUs" are one axis
# and the node1 -> node2 jump conflates twice-the-GPUs with the all-to-all
# leaving the node.  Here a node holds 4 H100 on NVLink (~900 GB/s) and nodes
# are joined by 4x NDR200 (~9:1 slower), so those two things separate: 1-4 GPUs
# stay inside one node, 5+ GPUs cannot.
#
# The axis is therefore "GPUs allocated", with P3DFFT placed at 4 x (its node
# count).  That is not a fudge factor -- it is the allocation: one TSUBAME node
# is 4 H100 to cuFFTMp and 160 EPYC cores to P3DFFT, and 1 and 2 GPUs are the
# node_q and node_h resource types, i.e. literally a quarter and a half node.
# So a given x is the same hardware rental either way, which is what makes the
# ratio panel a fair GPU-to-CPU comparison.  The shaded band is x <= 4, the
# single-node region.
#
# Median over trials with min-max whiskers (referee R1-4 variability).
# Output: cuFFTMp-scaling-tsubame.pdf
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

ROOT   = 'results-202607/tsubame/strong_scaling'
# The grids are the powers of two Miyabi used, so the two machines can be read
# against each other directly.  The first campaign ran 768^3 / 1536^3 to buy the
# 3- and 12-GPU points, and that turned out to cost more than the points were
# worth: measured on one H100 with the same binary, 768^3 sustains 10.1 Gpoint/s
# against 15.8 at 512^3 and 16.1 at 1024^3, i.e. the single factor of 3 makes it
# 1.66x more expensive per unit N log2 N.  Worse, the penalty is not a constant
# that divides out of a strong-scaling curve: it lands on the arithmetic only,
# so once the all-to-all dominates it falls to ~1.15x, and the 1 -> 64 GPU
# speedup read off a 768^3 curve is inflated by about 1.4x for that reason
# alone.  The 768^3 / 1536^3 data stays on disk as the measurement of that.
grids  = ['1024', '2048']
LIBS   = ['cuFFTMp-slab', 'cuFFTMp-pencil', 'P3DFFT']
DIR    = {'cuFFTMp-slab': 'cufftmp-slab',
          'cuFFTMp-pencil': 'cufftmp-pencil',
          'P3DFFT': 'p3dfft'}
MARKER = {'cuFFTMp-slab': 'o', 'cuFFTMp-pencil': '^', 'P3DFFT': 's'}
# Per panel, because the two grids leave different corners empty: at N=1024 every
# curve runs along the bottom, at N=2048 the CPU curve runs along the top.  Not
# 'best' -- that would move the legend whenever a point is added, and the two
# panels have to be read side by side.
LEGLOC = {'1024': 'upper right', '2048': 'lower left'}
GPN    = 4              # H100 per node
INTRA  = 4              # x <= INTRA is one node: NVLink only, no InfiniBand


def read_run(path):
    txt = open(path, errors='ignore').read().replace('\0', '')
    m = re.search(r'nx\s*=\s*(\d+)', txt)
    nx = int(m.group(1)) if m else None
    hit = [x for x in txt.splitlines() if 'cpu time per loop =' in x]
    if not hit:
        raise ValueError('no loop-time marker')
    return nx, float(hit[0].split('=')[2])


def collect(lib, grid):
    """key = GPUs allocated.  P3DFFT has no GPUs, so its node count is
    multiplied by GPN: both codes are then indexed by the same rental."""
    d = os.path.join(DIR[lib], ROOT, f'{grid}^3')
    out = defaultdict(list)
    if not os.path.isdir(d):
        return out
    for f in os.listdir(d):
        m = re.fullmatch(r'gpu(\d+)-node(\d+)-(\d+)', f)
        if m:
            x, nd = int(m.group(1)), int(m.group(2))
            # The placement control runs (step 14 of the plan) share this
            # directory and this naming: gpu4-node2 is four GPUs spread two per
            # node, which is a different machine from gpu4-node1 even though the
            # GPU count matches -- 0.0466 s against 0.0133 s when this was
            # measured at N=768.  The strong-scaling curve is the packed
            # placement only, so keep the
            # files whose node count is the minimum that holds x GPUs and let
            # the placement study read the rest.
            if nd != -(-x // GPN):
                continue
        else:
            m = re.fullmatch(r'rank(\d+)-node(\d+)-(\d+)', f)
            if not m:
                continue
            x = GPN * int(m.group(2))
        try:
            nx, t = read_run(os.path.join(d, f))
        except (IndexError, ValueError):
            continue
        # Name the grid, verify the grid: the header is written by the binary,
        # the file name by the driver, and only one of the two can be wrong
        # without anybody noticing (lessons L18).
        if nx is not None and nx != int(grid):
            continue
        out[x].append(t)
    return out


fig = plt.figure(figsize=(16, 12.0))
# The ratio row is the shorter one: it carries a single dimensionless number per
# point, whereas the absolute row has to resolve three curves over three decades.
gs  = fig.add_gridspec(2, 2, height_ratios=[1.5, 0.8],
                       left=0.105, right=0.99, bottom=0.125, top=0.905,
                       hspace=0.05, wspace=0.04)
# absolute row first, so the ratio row below can inherit its x-axis
axes  = [fig.add_subplot(gs[0, 0])]
axes += [fig.add_subplot(gs[0, 1], sharey=axes[0])]
axr   = [fig.add_subplot(gs[1, 0], sharex=axes[0])]
axr  += [fig.add_subplot(gs[1, 1], sharex=axes[1], sharey=axr[0])]

for ax, axrat, grid in zip(axes, axr, grids):
    med_of = {}          # lib -> {x: (median, min, max)}, for the ratio panel
    for ilib, lib in enumerate(LIBS):
        data = collect(lib, grid)
        if lib == 'cuFFTMp-pencil':
            data = {n: v for n, v in data.items()
                    if int(n ** 0.5 + 0.5) ** 2 == n}
        xs = sorted(data)
        if not xs:
            print(f'  -- {lib} N={grid}: no data')
            continue
        # The protocol claimed in the paper is three independent submissions per
        # point.  read_run() already drops files whose header does not match the
        # grid, so a launch failure leaves a name on disk but no sample -- and a
        # point quietly plotted from one or two runs would still look like every
        # other point.  Say the count out loud for every point instead, so that
        # any shortfall has to be either fixed or disclosed, never overlooked.
        short = [(n, len(data[n])) for n in xs if len(data[n]) != 3]
        if short:
            print(f'  !! {lib} N={grid}: trials != 3 at ' +
                  ', '.join(f'gpu{n}(n={c})' for n, c in short))
        med = [median(data[n]) for n in xs]
        lo  = [median(data[n]) - min(data[n]) for n in xs]
        hi  = [max(data[n]) - median(data[n]) for n in xs]
        med_of[lib] = {n: (median(data[n]), min(data[n]), max(data[n]))
                       for n in xs}
        print(f'  {lib} N={grid}: ' +
              ', '.join(f'{n}={v:.4g}s' for n, v in zip(xs, med)))

        if lib == 'cuFFTMp-slab':   # ideal strong-scaling guide
            ax.loglog([xs[0], xs[-1]],
                      [med[0], med[0] * xs[0] / xs[-1]], 'k-')

        # capsize is a half-length: matplotlib draws the cap as a '_' marker of
        # size 2*capsize, so capsize=12 gives a 24 pt bar against a 15 pt marker
        # and the cap stays visible when the spread is smaller than the symbol.
        ax.errorbar(xs, med, yerr=[lo, hi], lw=3, ms=15, mew=3,
                    marker=MARKER[lib], ls='-', mfc=colors[ilib],
                    color=colors[ilib], capsize=12, capthick=linewidth,
                    label=lib)

    for a in (ax, axrat):
        a.axvspan(0.5, INTRA, color='0.85', zorder=0, lw=0)
    ax.set_xscale('log'); ax.set_yscale('log')
    # The campaign stops at 64 GPUs / 16 nodes, so the x range no longer needs
    # to reach 128; 90 keeps the last point off the frame edge.  The y range is
    # taken from the measurements rather than inherited: 0.021 s (slab 1024^3
    # on 32 GPUs) to 2.97 s (P3DFFT 2048^3 on two nodes).  The old [2e-3, 5]
    # was set for 768^3/1536^3 and left the bottom decade empty.
    ax.set_xlim([0.8, 90])
    ax.set_ylim([8e-3, 5])
    leg = ax.legend(frameon=False, loc=LEGLOC[grid], fontsize=22,
                    columnspacing=0.6, title=rf'$N = {grid}$')
    leg.get_title().set_fontsize(22)
    plt.setp(ax.get_xticklabels(), visible=False)

    # ---- ratio panel: how many times faster than the CPU baseline ----------
    # The whiskers are the widest ratio the samples permit -- slowest CPU over
    # fastest GPU and vice versa -- rather than a propagated error, because with
    # three samples the extremes are the honest statement of what was observed.
    cpu = med_of.get('P3DFFT', {})
    for ilib, lib in enumerate(['cuFFTMp-slab', 'cuFFTMp-pencil']):
        gpu = med_of.get(lib, {})
        ns  = [n for n in sorted(gpu) if n in cpu]
        if not ns:
            print(f'  -- ratio P3DFFT/{lib} N={grid}: no shared allocation')
            continue
        r   = [cpu[n][0] / gpu[n][0] for n in ns]
        rlo = [r[i] - cpu[n][1] / gpu[n][2] for i, n in enumerate(ns)]
        rhi = [cpu[n][2] / gpu[n][1] - r[i] for i, n in enumerate(ns)]
        print(f'  ratio P3DFFT/{lib} N={grid}: ' +
              ', '.join(f'gpu{n}={v:.3g}' for n, v in zip(ns, r)))
        # Pencil is computed and printed but not drawn, as in the Miyabi figure:
        # the solver uses slab, and the pencil ratio hugs parity closely enough
        # that plotting it forces the ordinate below 1 and squashes the slab
        # curve, which is the one the text quotes.  The pencil-versus-slab
        # comparison stays in the top row, where it belongs.
        if lib != 'cuFFTMp-slab':
            continue
        axrat.errorbar(ns, r, yerr=[rlo, rhi], lw=3, ms=15, mew=3,
                       marker=MARKER[lib], ls='-', mfc=colors[ilib],
                       color=colors[ilib], capsize=12, capthick=linewidth,
                       label=rf'P3DFFT / {lib}')
    # No parity line: the ordinate starts at 1, so the bottom spine is parity and
    # a dotted line there would be drawn underneath it.  (Same convention as the
    # Miyabi figures, so the two machines can be read against each other.)
    axrat.set_xscale('log'); axrat.set_yscale('log')
    axrat.set_ylim([1.0, 40])

    # Node count on the top spine.  Both codes are billed by the node, so this
    # is the axis that says where the all-to-all has to cross InfiniBand; the
    # sub-node allocations (1-3 GPUs) fall below its first tick by construction.
    sec = ax.secondary_xaxis('top', functions=(lambda g: g / GPN,
                                               lambda n: n * GPN))
    sec.set_xscale('log')
    sec.xaxis.set_major_locator(FixedLocator([1, 2, 4, 8, 16, 32]))
    sec.xaxis.set_major_formatter(
            FixedFormatter(['1', '2', '4', '8', '16', '32']))
    sec.xaxis.set_minor_formatter(NullFormatter())
    sec.tick_params(direction='in', width=linewidth, length=17,
                    labelsize=fontsize)

    # Last, not with the other x-axis settings above: set_xscale() propagates
    # to every axes sharing the axis and resets its locators, so anything set
    # before axrat.set_xscale() is silently thrown away.  The GPU counts are now
    # powers of two throughout (the 3- and 12-GPU points went away with the
    # 768^3 grid), and the decade ticks a log scale defaults to would label none
    # of the points that matter here.
    ax.xaxis.set_major_locator(FixedLocator([1, 4, 16, 64]))
    ax.xaxis.set_major_formatter(FixedFormatter(['1', '4', '16', '64']))
    ax.xaxis.set_minor_locator(FixedLocator([2, 8, 32, 128]))
    ax.xaxis.set_minor_formatter(NullFormatter())

axr[0].set_ylabel(r'Speedup over P3DFFT')
axr[1].legend(frameon=False, loc='upper left', fontsize=22, columnspacing=0.6)
plt.setp(axr[1].get_yticklabels(), visible=False)
plt.setp(axes[1].get_yticklabels(), visible=False)
axes[0].set_ylabel(r'Time per FFT + iFFT [s]')
# The two rows carry tick labels of different widths, so matplotlib places
# their y-labels at different x by default and the pair reads as misaligned.
# align_ylabels moves both to the leftmost of the two positions.
fig.align_ylabels([axes[0], axr[0]])
fig.text(0.55, 0.028, r'Number of GPUs (P3DFFT: $4\times$ number of nodes)',
         ha='center', va='bottom', fontsize=fontsize)
fig.text(0.55, 0.962, r'Number of nodes', ha='center', va='bottom',
         fontsize=fontsize)
plt.savefig('cuFFTMp-scaling-tsubame.pdf')
print('wrote cuFFTMp-scaling-tsubame.pdf')

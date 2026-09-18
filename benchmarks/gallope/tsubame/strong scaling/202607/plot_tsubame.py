# -*- coding: utf-8 -*-
# Fig. 2, TSUBAME4.0 version (2026-07 campaign).  Strong scaling of the solver
# itself at 512^3: absolute time per step on top, speedup over Calliope below.
#
# 512^3, not the 384^3 of the first campaign.  384 = 2^7 * 3 was chosen because
# the factor of 3 admits 3- and 12-GPU allocations, and the factor of 3 turned
# out to cost more than those two points are worth: on one H100, the same FFT
# binary sustains 10.1 Gpoint/s at 768^3 against 15.8 at 512^3 and 16.1 at
# 1024^3.  The penalty lands on the arithmetic only, so it does not divide out
# of a strong-scaling curve -- it shrinks as the all-to-all takes over, which
# inflates the apparent speedup.  Powers of two throughout now, at the price of
# the 3- and 12-GPU points; 1, 2, 4 GPUs still resolve the intra-node region.
#
# The abscissa is "GPUs allocated", with Calliope placed at 4 x (its node count),
# for the same reason as the FFT figure: one TSUBAME node is 4 H100 to Gallope
# and 160 EPYC cores to Calliope, and the 1- and 2-GPU points are the node_q and
# node_h resource types, i.e. literally a quarter and a half node.  A given x is
# therefore the same hardware rental for both codes, which is what makes the
# ratio panel a fair GPU-to-CPU comparison.  The shaded band is x <= 4, the
# single-node region where the all-to-all never leaves NVLink.
#
# That band is the whole point of repeating the measurement here.  On Miyabi-G a
# node holds one GPU, so every Gallope point already paid the inter-node cost and
# the figure could not separate "more GPUs" from "communication left the node".
#
# Median over trials with min-max whiskers (referee R1-4 variability).
# Output: gallope-scaling-tsubame.pdf
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
rcParams.update({'figure.autolayout': False})   # margins are set on the gridspec
from latex_preamble import preamble
rcParams['text.latex.preamble'] = preamble
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
#-------------------------------------------------------------#

import os
import re
from collections import defaultdict
from statistics import median

ROOT  = '512^3-202607'
GRID  = 512
GPN   = 4               # H100 per node
INTRA = 4               # x <= INTRA is one node: NVLink only, no InfiniBand

STYLE = {
    'gallope':  dict(marker='o', label=r'\textsc{Gallope}'),
    'calliope': dict(marker='s', label=r'\textsc{Calliope}'),
}


def read_run(path):
    """(nlx, seconds per step) from a Calliope/Gallope stdout."""
    lines = [l.rstrip('\n') for l in open(path, errors='ignore')]
    m = re.search(r'nlx\s*=\s*(\d+)', '\n'.join(lines))
    nlx = int(m.group(1)) if m else None
    nloop = int([x for x in lines if '# of steps advanced' in x][0]
                .split('# of steps advanced')[1])
    adv = float([x for x in lines if 'Advance steps' in x][0]
                .split('Advance steps')[1].split('min')[0])
    return nlx, adv / nloop * 60.0


def collect(lib):
    """key = GPUs allocated.  Calliope has no GPUs, so its node count is
    multiplied by GPN: both codes are then indexed by the same rental."""
    d = os.path.join(ROOT, lib)
    out = defaultdict(list)
    if not os.path.isdir(d):
        return out
    for f in os.listdir(d):
        m = re.fullmatch(r'gpu(\d+)-node(\d+)-(\d+)', f)
        if m:
            x, nd = int(m.group(1)), int(m.group(2))
            # The placement control runs (step 14 of the plan) share this
            # directory and this naming: gpu4-node2 is four GPUs spread two per
            # node, a different machine from gpu4-node1 even though the GPU count
            # matches -- 0.193 s against 0.098 s per step.  The strong-scaling
            # curve is the packed placement only, so keep the files whose node
            # count is the minimum that holds x GPUs and let the placement study
            # read the rest.
            if nd != -(-x // GPN):
                continue
        else:
            m = re.fullmatch(r'rank(\d+)-node(\d+)-(\d+)', f)
            if not m:
                continue
            x = GPN * int(m.group(2))
        try:
            nlx, t = read_run(os.path.join(d, f))
        except (IndexError, ValueError):
            print(f'  !! unparsable, skipped: {d}/{f}')
            continue
        # Name the grid, verify the grid: the header is written by the solver,
        # the file name by the driver, and only one of the two can be wrong
        # without anybody noticing (lessons L18).
        if nlx is not None and nlx != GRID:
            print(f'  !! grid {nlx} != {GRID}, skipped: {d}/{f}')
            continue
        out[x].append(t)
    return out


fig = plt.figure(figsize=(9.5, 11.3))
# The ratio panel is the shorter one: it carries a single dimensionless number
# per point, whereas the absolute panel has to resolve two curves over decades.
gs  = fig.add_gridspec(2, 1, height_ratios=[1.55, 0.85],
                       left=0.165, right=0.985, bottom=0.145, top=0.925,
                       hspace=0.05)
ax    = fig.add_subplot(gs[0, 0])
axrat = fig.add_subplot(gs[1, 0], sharex=ax)

stat = {}          # lib -> {x: (median, min, max)}, for the ratio panel
for ilib, lib in enumerate(['gallope', 'calliope']):
    data = collect(lib)
    xs = sorted(data)
    if not xs:
        print(f'  -- {lib}: no data')
        continue
    # The protocol claimed in the paper is three independent submissions per
    # point.  read_run() already drops files whose header does not match the
    # grid, so a launch failure leaves a name on disk but no sample -- and a
    # point quietly plotted from one or two runs would look like every other
    # point.  Say the count out loud, so a shortfall has to be fixed or
    # disclosed, never overlooked.
    short = [(n, len(data[n])) for n in xs if len(data[n]) != 3]
    if short:
        print(f'  !! {lib}: trials != 3 at ' +
              ', '.join(f'gpu{n}(n={c})' for n, c in short))
    med = [median(data[n]) for n in xs]
    lo  = [median(data[n]) - min(data[n]) for n in xs]
    hi  = [max(data[n]) - median(data[n]) for n in xs]
    stat[lib] = {n: (median(data[n]), min(data[n]), max(data[n])) for n in xs}

    for n in xs:
        v = data[n]
        print(f'{lib:9s} gpu{n:<4d} n={len(v)}  median={median(v):.4f}  '
              f'min={min(v):.4f}  max={max(v):.4f}  '
              f'spread={100*(max(v)-min(v))/median(v):.1f}%')

    if lib == 'gallope':   # ideal strong-scaling guide
        ax.loglog([xs[0], xs[-1]], [med[0], med[0] * xs[0] / xs[-1]], 'k-')

    # capsize is a half-length: matplotlib draws the cap as a '_' marker of size
    # 2*capsize, so capsize=12 gives a 24 pt bar against a 15 pt marker and the
    # cap stays visible when the spread is smaller than the symbol.
    ax.errorbar(xs, med, yerr=[lo, hi], lw=3, ms=15, mew=3,
                marker=STYLE[lib]['marker'], ls='-', mfc=colors[ilib],
                color=colors[ilib], capsize=12, capthick=linewidth,
                label=STYLE[lib]['label'])

for a in (ax, axrat):
    a.axvspan(0.5, INTRA, color='0.85', zorder=0, lw=0)
leg = ax.legend(frameon=False, loc='lower left', fontsize=22,
                columnspacing=0.6, title=rf'$N = {GRID}$')
leg.get_title().set_fontsize(22)

# ---- ratio panel: how many times faster than the CPU code -------------------
# Whiskers are the widest ratio the samples permit (slowest CPU over fastest GPU
# and vice versa) rather than a propagated error: with three samples the observed
# extremes are the honest statement of the spread.
cpu = stat.get('calliope', {})
gpu = stat.get('gallope', {})
ns  = [n for n in sorted(gpu) if n in cpu]
if ns:
    r   = [cpu[n][0] / gpu[n][0] for n in ns]
    rlo = [r[i] - cpu[n][1] / gpu[n][2] for i, n in enumerate(ns)]
    rhi = [cpu[n][2] / gpu[n][1] - r[i] for i, n in enumerate(ns)]
    print('  ratio Calliope/Gallope: ' +
          ', '.join(f'gpu{n}={v:.3g}' for n, v in zip(ns, r)))
    axrat.errorbar(ns, r, yerr=[rlo, rhi], lw=3, ms=15, mew=3,
                   marker=STYLE['gallope']['marker'], ls='-', mfc=colors[0],
                   color=colors[0], capsize=12, capthick=linewidth)
else:
    print('  -- ratio: no shared allocation')
# No parity line: the ordinate starts at 1, so the bottom spine is parity and a
# dotted line there would be drawn underneath it.
axrat.set_xscale('log'); axrat.set_yscale('log')
axrat.set_ylim([1.0, 20])
axrat.set_ylabel(r'Speedup over \textsc{Calliope}')

ax.set_xscale('log'); ax.set_yscale('log')
# [2e-2, 2] was the 384^3 range and it clips this data: Calliope on one node is
# 2.60 s/step at 512^3, half a division above the old ceiling, so the leftmost
# CPU point was being drawn outside the axes.  The measured span is 0.094 s
# (Gallope on 32 GPUs) to 2.60 s.  x stops at 90 because the campaign stops at
# 64 GPUs -- 3 and 12 GPUs went away with 384^3, 128 with the 32-node tier.
ax.set_xlim([0.8, 90])
ax.set_ylim([4e-2, 4])
ax.set_ylabel(r'Time per step [s]')
plt.setp(ax.get_xticklabels(), visible=False)

# Node count on the top spine.  Both codes are billed by the node, so this is
# the axis that says where the all-to-all has to cross InfiniBand; the sub-node
# allocations (1-3 GPUs) fall below its first tick by construction.
sec = ax.secondary_xaxis('top', functions=(lambda g: g / GPN,
                                           lambda n: n * GPN))
sec.set_xscale('log')
sec.xaxis.set_major_locator(FixedLocator([1, 2, 4, 8, 16, 32]))
sec.xaxis.set_major_formatter(FixedFormatter(['1', '2', '4', '8', '16', '32']))
sec.xaxis.set_minor_formatter(NullFormatter())
sec.tick_params(direction='in', width=linewidth, length=17, labelsize=fontsize)

# Last, not with the other x-axis settings above: set_xscale() propagates to
# every axes sharing the axis and resets its locators, so anything set before
# axrat.set_xscale() is silently thrown away.  The GPU counts are powers of two
# throughout (3 and 12 went away with 384^3), and the decade ticks a log scale
# defaults to would label none of the points here.
ax.xaxis.set_major_locator(FixedLocator([1, 4, 16, 64]))
ax.xaxis.set_major_formatter(FixedFormatter(['1', '4', '16', '64']))
ax.xaxis.set_minor_locator(FixedLocator([2, 8, 32, 128]))
ax.xaxis.set_minor_formatter(NullFormatter())

axrat.set_xlabel(r'Number of GPUs (\textsc{Calliope}: $4\times$ nodes)')
fig.text(0.575, 0.962, r'Number of nodes', ha='center', va='bottom',
         fontsize=fontsize)
# The two panels carry tick labels of different widths, so matplotlib places
# their y-labels at different x by default and the pair reads as misaligned.
fig.align_ylabels([ax, axrat])
plt.savefig('gallope-scaling-tsubame.pdf')
print('wrote gallope-scaling-tsubame.pdf')

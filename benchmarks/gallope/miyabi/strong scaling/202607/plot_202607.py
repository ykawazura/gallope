# -*- coding: utf-8 -*-
# Fig. 2 (2026-07 campaign).  Reads nodeN-t files (t = trial index) and plots
# the median over trials with min-max whiskers, so the referee's request for
# run-to-run variability (R1-4) is answered by the figure itself.
# Two stacked panels: the upper one is the absolute time per step (as originally
# submitted), the lower one the GPU-to-CPU speedup, i.e. the ratio of the
# Calliope time to the Gallope time at the same node count.  The text quotes
# that ratio node by node (3.4x at 16 nodes to 6.6x at 64), so it should be
# readable off the figure instead of reconstructed by eye from two log curves.
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

ROOT  = '1024^3-202607'
grids = ['1024']

STYLE = {
    'gallope':         dict(marker='o', label=r'\textsc{Gallope} (separated)'),
    'gallope-unified': dict(marker='^', label=r'\textsc{Gallope} (unified)'),
    'calliope':        dict(marker='s', label=r'\textsc{Calliope}'),
}


def step_time(path):
    """Seconds per step from a Calliope/Gallope stdout."""
    lines = [l.rstrip('\n') for l in open(path)]
    nloop = int([x for x in lines if '# of steps advanced' in x][0]
                .split('# of steps advanced')[1])
    adv = float([x for x in lines if 'Advance steps' in x][0]
                .split('Advance steps')[1].split('min')[0])
    return adv / nloop * 60.0


def collect(d):
    """{node count: [time, ...]} over all trials found in directory d."""
    out = defaultdict(list)
    for f in os.listdir(d):
        m = re.fullmatch(r'node(\d+)(?:-(\d+))?', f)
        if not m:
            continue
        try:
            out[int(m.group(1))].append(step_time(os.path.join(d, f)))
        except (IndexError, ValueError):
            print(f'  !! unparsable, skipped: {d}/{f}')
    return out


for grid in grids:
    fig = plt.figure(figsize=(9.5, 11.3))
    # The ratio panel is the shorter one: it carries a single dimensionless
    # number per point, whereas the absolute panel has to resolve three curves.
    gs  = fig.add_gridspec(2, 1, height_ratios=[1.55, 0.85],
                           left=0.155, right=0.985, bottom=0.135, top=0.99,
                           hspace=0.05)
    ax    = fig.add_subplot(gs[0, 0])          # absolute time (as submitted)
    axrat = fig.add_subplot(gs[1, 0], sharex=ax)

    stat = {}          # lib -> {node: (median, min, max)}, for the ratio panel
    for ilib, lib in enumerate(['gallope', 'gallope-unified', 'calliope']):
        data = collect(os.path.join(ROOT, lib))
        nodes = sorted(data)
        med = [median(data[n]) for n in nodes]
        lo  = [median(data[n]) - min(data[n]) for n in nodes]
        hi  = [max(data[n]) - median(data[n]) for n in nodes]
        stat[lib] = {n: (median(data[n]), min(data[n]), max(data[n]))
                     for n in nodes}

        for n in nodes:
            v = data[n]
            print(f'{lib:16s} node{n:<4d} n={len(v)}  median={median(v):.4f}  '
                  f'min={min(v):.4f}  max={max(v):.4f}  spread={100*(max(v)-min(v))/median(v):.1f}%')

        if lib == 'gallope':   # ideal strong-scaling guide, anchored at nodes[1]
            ax.loglog([nodes[1], nodes[-1]],
                      [med[1], med[1] * nodes[1] / nodes[-1]], 'k-')

        # capsize is a half-length: matplotlib draws the cap as a '_' marker of
        # size 2*capsize, so capsize=12 gives a 24 pt bar against a 15 pt marker
        # and the cap stays visible when the spread is smaller than the symbol.
        ax.errorbar(nodes, med, yerr=[lo, hi], lw=3, ms=15, mew=3,
                    marker=STYLE[lib]['marker'], ls='-', mfc=colors[ilib],
                    color=colors[ilib], capsize=12, capthick=linewidth,
                    label=STYLE[lib]['label'])

    leg = ax.legend(frameon=False, loc='lower left', fontsize=22,
                    columnspacing=0.6, title=rf'$N = {grid}$')
    leg.get_title().set_fontsize(22)

    # ---- ratio panel: how many times faster than the CPU code ---------------
    # Whiskers are the widest ratio the samples permit (slowest CPU over fastest
    # GPU and vice versa) rather than a propagated error: with three samples the
    # observed extremes are the honest statement of the spread.
    cpu = stat['calliope']
    for ilib, lib in enumerate(['gallope', 'gallope-unified']):
        gpu = stat[lib]
        ns  = [n for n in sorted(gpu) if n in cpu]
        if not ns:
            continue
        r   = [cpu[n][0] / gpu[n][0] for n in ns]
        rlo = [r[i] - cpu[n][1] / gpu[n][2] for i, n in enumerate(ns)]
        rhi = [cpu[n][2] / gpu[n][1] - r[i] for i, n in enumerate(ns)]
        print('  ratio Calliope/' + lib + ': ' +
              ', '.join(f'node{n}={v:.3g}' for n, v in zip(ns, r)))
        axrat.errorbar(ns, r, yerr=[rlo, rhi], lw=3, ms=15, mew=3,
                       marker=STYLE[lib]['marker'], ls='-', mfc=colors[ilib],
                       color=colors[ilib], capsize=12, capthick=linewidth,
                       label=r'\textsc{Calliope} / ' + STYLE[lib]['label'])
    # No parity line: the ordinate now starts at 1, so the bottom spine is
    # parity and a dotted line there would be drawn underneath it.
    axrat.set_xscale('log'); axrat.set_yscale('log')
    axrat.set_ylim([1.0, 55])
    axrat.set_ylabel(r'Speedup over \textsc{Calliope}')
    axrat.legend(frameon=False, loc='upper left', fontsize=20,
                 columnspacing=0.6, labelspacing=0.3)

    ax.set_xlim([0.9, 400])
    ax.set_ylim([1e-1, 7])
    ax.set_ylabel(r'Time per step [s]')
    plt.setp(ax.get_xticklabels(), visible=False)
    axrat.set_xlabel('Number of nodes\n(= number of CPUs = number of GPUs)')
    # The two panels carry tick labels of different widths, so matplotlib
    # places their y-labels at different x by default and the pair reads as
    # misaligned.  align_ylabels moves both to the leftmost of the two.
    fig.align_ylabels([ax, axrat])
    plt.savefig(grid + '-202607.pdf')

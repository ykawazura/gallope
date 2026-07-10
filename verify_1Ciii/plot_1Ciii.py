#!/usr/bin/env python3
# ------------------------------------------------------------------
# 1-C-iii g-forcing validation figures (run LOCALLY after pulling the
# NetCDF from Miyabi). Produces PNGs + a numeric summary for the report:
#   fig_W_free.png : W_free(t) for lin64 & echo64 -> driven steady state
#                    (plateau, NOT the passive run's monotonic decay).
#   fig_W_m.png    : steady-avg Hermite spectrum W_m(m), log-log, with
#                    m^-1/2 (phase-mixing) and m^-1 (fluidized) guides.
#   fig_Gamma_m.png: cumulative Gamma_m(m) at mid-inertial k_perp and
#                    ECHO_RATIO(k_perp) -> forward cascade (lin, ~1) vs
#                    stochastic echo (echo, <<1).
#
# Usage:
#   plot_1Ciii.py OUTDIR T_STEADY lin64=PATH.nc echo64=PATH.nc [...]
# ------------------------------------------------------------------
import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from netCDF4 import Dataset

COLORS = {'lin64': '#c2410c', 'echo64': '#1d4ed8'}   # burnt orange / indigo


def load(path):
    ds = Dataset(path)
    v = ds.variables
    out = {
        'tt': np.asarray(v['tt'][:], float),
        'mm': np.asarray(v['mm'][:], float),
        'kpbin': np.asarray(v['kpbin'][:], float),
        'W_free': np.asarray(v['W_free'][:], float),
    }
    if 'W_m' in v:
        out['W_m'] = np.asarray(v['W_m'][:], float)          # (tt, mm)
    if 'Gamma_m' in v:
        out['Gamma_m'] = np.asarray(v['Gamma_m'][:], float)  # (tt, kz, mm, kpbin)
    if 'Gamma_m_kint' in v:
        out['Gamma_m_kint'] = np.asarray(v['Gamma_m_kint'][:], float)  # (tt, mm)
    ds.close()
    return out


def steady_idx(tt, tsteady):
    idx = np.where(tt >= tsteady)[0]
    if idx.size < 3:
        idx = np.arange(max(1, len(tt) // 2), len(tt))
    return idx


def inertial_bins(kpbin):
    kp_max = float(kpbin[-1])
    ib = np.where((kpbin > 2.0) & (kpbin < kp_max / 3.0))[0]
    if ib.size == 0:
        ib = np.arange(2, max(3, len(kpbin) // 3))
    return ib


def main():
    outdir = sys.argv[1]
    tsteady = float(sys.argv[2])
    runs = {}
    for a in sys.argv[3:]:
        label, path = a.split('=', 1)
        runs[label] = load(path)
    os.makedirs(outdir, exist_ok=True)

    summary = []

    # ---------- fig 1: W_free(t) ----------
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    for label, d in runs.items():
        tt, Wf = d['tt'], d['W_free']
        ax.plot(tt, Wf, color=COLORS.get(label, None), lw=1.6, label=label)
        idx = steady_idx(tt, tsteady)
        m, s = float(np.mean(Wf[idx])), float(np.std(Wf[idx]))
        relrms = s / max(m, 1e-300)
        ax.axhline(m, color=COLORS.get(label, None), ls='--', lw=0.8, alpha=0.6)
        summary.append((label, 'W_free', Wf[0], Wf[-1], m, s, relrms))
    ax.axvline(tsteady, color='0.5', ls=':', lw=1.0)
    ax.text(tsteady, ax.get_ylim()[1], ' steady window', va='top', ha='left',
            fontsize=8, color='0.4')
    ax.set_xlabel('time  $t$')
    ax.set_ylabel(r'compressive free energy  $W_{\rm free}$')
    ax.set_title(r'Driven steady state: $W_{\rm free}(t)$ plateaus (g-forcing on)')
    ax.legend(frameon=False)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, 'fig_W_free.png'), dpi=140)
    plt.close(fig)

    # ---------- fig 2: W_m(m) steady spectrum ----------
    if all('W_m' in d for d in runs.values()):
        fig, ax = plt.subplots(figsize=(6.4, 4.6))
        for label, d in runs.items():
            idx = steady_idx(d['tt'], tsteady)
            Wm = d['W_m'][idx].mean(axis=0)
            m = d['mm']
            ax.loglog(m[1:], Wm[1:], 'o-', ms=3.5, color=COLORS.get(label, None),
                      label=label)
        # reference slopes anchored at m=2
        mref = np.array([2.0, m[-1]])
        i2 = 1
        for lbl, p in [(r'$m^{-1/2}$', -0.5), (r'$m^{-1}$', -1.0)]:
            base = list(runs.values())[0]
            idx = steady_idx(base['tt'], tsteady)
            Wm0 = base['W_m'][idx].mean(axis=0)
            anchor = Wm0[i2]
            ax.loglog(mref, anchor * (mref / mref[0])**p, 'k', ls='--' if p == -1 else ':',
                      lw=1.0, alpha=0.7, label=lbl)
        ax.set_xlabel('Hermite index  $m$')
        ax.set_ylabel(r'$W_m$  (steady-avg)')
        ax.set_title('Hermite free-energy spectrum')
        ax.legend(frameon=False, fontsize=9)
        ax.grid(alpha=0.25, which='both')
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, 'fig_W_m.png'), dpi=140)
        plt.close(fig)

    # ---------- fig 3: Gamma_m cumulative + ECHO_RATIO ----------
    if all('Gamma_m' in d for d in runs.values()):
        fig, axs = plt.subplots(1, 2, figsize=(11.0, 4.4))
        for label, d in runs.items():
            idx = steady_idx(d['tt'], tsteady)
            G = d['Gamma_m'][idx].mean(axis=0).sum(axis=0)   # (mm, kpbin)
            kpbin, mm = d['kpbin'], d['mm']
            ib = inertial_bins(kpbin)
            # left: cumulative Gamma_m(m) at mid-inertial k_perp
            ikmid = ib[ib.size // 2]
            prof = G[:, ikmid]
            axs[0].plot(mm, prof, 'o-', ms=3, color=COLORS.get(label, None),
                        label='%s ($k_\\perp$=%.1f)' % (label, float(kpbin[ikmid])))
            # right: ECHO_RATIO(k_perp) over inertial band
            ratios = []
            for ik in ib:
                cum = G[:, ik]
                src = np.diff(cum, prepend=0.0)
                gross = float(np.sum(np.abs(src)))
                net = float(np.max(np.abs(cum)))
                ratios.append(net / gross if gross > 1e-300 else np.nan)
            ratios = np.array(ratios)
            axs[1].plot(kpbin[ib], ratios, 'o-', ms=3.5, color=COLORS.get(label, None),
                        label=label)
            summary.append((label, 'ECHO_RATIO_mean',
                            float(np.nanmean(ratios)), float(np.nanmin(ratios)),
                            0, 0, 0))
        axs[0].axhline(0, color='0.5', lw=0.8)
        axs[0].set_xlabel('Hermite index  $m$')
        axs[0].set_ylabel(r'cumulative flux  $\Gamma_m(k_\perp)$')
        axs[0].set_title('Hermite flux profile (mid-inertial $k_\\perp$)')
        axs[0].legend(frameon=False, fontsize=8)
        axs[0].grid(alpha=0.25)
        axs[1].axhline(1.0, color='0.6', ls=':', lw=1.0)
        axs[1].set_xlabel(r'$k_\perp$')
        axs[1].set_ylabel(r'ECHO_RATIO $=\max_m|\Gamma_m|/\sum_m|{\rm src}_m|$')
        axs[1].set_title('Echo metric (forward cascade $\\to$1, echo $\\to$0)')
        axs[1].set_ylim(0, 1.15)
        axs[1].legend(frameon=False, fontsize=9)
        axs[1].grid(alpha=0.25)
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, 'fig_Gamma_m.png'), dpi=140)
        plt.close(fig)

    # ---------- numeric summary ----------
    print("LABEL           METRIC              first/mean    last/min      mean       rms       relrms")
    for row in summary:
        label, metric, a, b, m, s, r = row
        print("%-14s  %-18s  %.4e  %.4e  %.4e  %.4e  %.4f" % (label, metric, a, b, m, s, r))


if __name__ == '__main__':
    main()

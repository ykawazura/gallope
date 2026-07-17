#!/usr/bin/env python3
# ------------------------------------------------------------------
# Alfvenic (electromagnetic) perpendicular energy spectra, from the
# SAME NetCDF as the g diagnostics -- no extra run needed.
#   E_u(k_perp) = 1/2 |u_perp|^2 = 1/2 k_perp^2 |Phi|^2  (kinetic / E-field side)
#   E_b(k_perp) = 1/2 |b_perp|^2 = 1/2 k_perp^2 |Psi|^2  (magnetic side)
#   E_zp/E_zm   = Elsasser z^+/z^- spectra
# kz is summed (full k_par) -> 1D perp spectrum, steady-time averaged.
#
# Usage:
#   plot_em_spectrum.py OUTDIR T_STEADY echo64=PATH.nc [lin64=PATH.nc ...]
# ------------------------------------------------------------------
import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from netCDF4 import Dataset

COLORS = {'lin64': '#c2410c', 'echo64': '#1d4ed8'}


def load(path):
    ds = Dataset(path)
    v = ds.variables
    out = {
        'tt': np.asarray(v['tt'][:], float),
        'kpbin': np.asarray(v['kpbin'][:], float),
        'upe2': np.asarray(v['upe2_bin'][:], float),   # (tt, kz, kpbin)
        'bpe2': np.asarray(v['bpe2_bin'][:], float),
        'zppe2': np.asarray(v['zppe2_bin'][:], float),
        'zmpe2': np.asarray(v['zmpe2_bin'][:], float),
    }
    ds.close()
    return out


def steady_mean_kperp(d, key, tsteady):
    tt = d['tt']
    idx = np.where(tt >= tsteady)[0]
    if idx.size < 3:
        idx = np.arange(max(1, len(tt) // 2), len(tt))
    E = d[key][idx].mean(axis=0).sum(axis=0)   # avg in t, sum over kz -> (kpbin,)
    # hide FFT round-off noise floor (empty modes ~1e-100) so it does not
    # blow up the log axis; a real spectrum never dips below ~1e-8 here.
    E = np.where(E < 1e-12, np.nan, E)
    return E


def add_slope(ax, kp, E, p, label, anchor_k=4.0, ls='--'):
    """Reference power law k^p anchored to E at k ~ anchor_k."""
    m = (kp >= anchor_k) & (E > 0)
    if not np.any(m):
        return
    k0 = kp[m][0]
    e0 = E[m][0]
    kk = kp[(kp >= k0)]
    ax.loglog(kk, e0 * (kk / k0) ** p, 'k', ls=ls, lw=1.0, alpha=0.6, label=label)


def main():
    outdir = sys.argv[1]
    tsteady = float(sys.argv[2])
    runs = {}
    for a in sys.argv[3:]:
        label, path = a.split('=', 1)
        runs[label] = load(path)
    os.makedirs(outdir, exist_ok=True)

    # ---- panel A: u vs b ;  panel B: Elsasser z+/z- ----
    fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.6))

    for label, d in runs.items():
        kp = d['kpbin']
        c = COLORS.get(label, None)
        Eu = steady_mean_kperp(d, 'upe2', tsteady)
        Eb = steady_mean_kperp(d, 'bpe2', tsteady)
        Ezp = steady_mean_kperp(d, 'zppe2', tsteady)
        Ezm = steady_mean_kperp(d, 'zmpe2', tsteady)
        mpos = kp > 0
        axs[0].loglog(kp[mpos], Eu[mpos], 'o-', ms=3.5, color=c, label=f'{label}  $E_u$')
        axs[0].loglog(kp[mpos], Eb[mpos], 's--', ms=3.0, color=c, alpha=0.65,
                      label=f'{label}  $E_b$')
        axs[1].loglog(kp[mpos], Ezp[mpos], 'o-', ms=3.5, color=c, label=f'{label}  $E^+$')
        axs[1].loglog(kp[mpos], Ezm[mpos], 's--', ms=3.0, color=c, alpha=0.65,
                      label=f'{label}  $E^-$')
        # numeric peek
        print(f'--- {label} (steady t>={tsteady}) ---')
        with np.printoptions(precision=3):
            print(' kperp :', kp[mpos][:12])
            print(' E_u   :', Eu[mpos][:12])
            print(' E_b   :', Eb[mpos][:12])

    # reference slopes anchored to the first run's E_u
    base = list(runs.values())[0]
    kp = base['kpbin']
    Eu = steady_mean_kperp(base, 'upe2', tsteady)
    add_slope(axs[0], kp, Eu, -3.0 / 2.0, r'$k_\perp^{-3/2}$', ls=':')
    add_slope(axs[0], kp, Eu, -5.0 / 3.0, r'$k_\perp^{-5/3}$', ls='--')

    axs[0].set_xlabel(r'$k_\perp$')
    axs[0].set_ylabel(r'$E(k_\perp)$  (steady-avg, $\sum_{k_z}$)')
    axs[0].set_title(r'Alfvenic energy: kinetic $E_u$ vs magnetic $E_b$')
    axs[0].set_ylim(1e-6, 3e2)
    axs[0].legend(frameon=False, fontsize=8)
    axs[0].grid(alpha=0.25, which='both')

    axs[1].set_xlabel(r'$k_\perp$')
    axs[1].set_ylabel(r'$E^\pm(k_\perp)$')
    axs[1].set_title(r'Elsasser spectra $E^+$, $E^-$')
    axs[1].set_ylim(1e-6, 3e2)
    axs[1].legend(frameon=False, fontsize=8)
    axs[1].grid(alpha=0.25, which='both')

    fig.tight_layout()
    out = os.path.join(outdir, 'fig_EM_kperp.png')
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print('wrote', out)


if __name__ == '__main__':
    main()

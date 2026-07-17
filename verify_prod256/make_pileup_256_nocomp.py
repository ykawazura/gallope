#!/usr/bin/env python3
# PILE-UP DIAGNOSTIC, raw spectra only (no compensated panels), latest data.
# Layout requested by user (3 rows x 2 cols, bottom-right blank):
#     Row 1:  (a) EM perp E(k_perp)      (b) EM parallel E(k_z)
#     Row 2:  (c) g  perp E_g(k_perp)    (d) g  parallel E_g(k_z)
#     Row 3:  (e) g  Hermite E_g(m)      [ blank ]
# Same raw-spectrum styling / guides / reductions as make_pileup_256.py.
# Time-mean = steady 2nd half of the records.
#
# Usage: python3 make_pileup_256_nocomp.py [NC] [outpng] [label]
import os
import sys
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from netCDF4 import Dataset

# Data + figure live next to this script (self-contained reproduction).
HERE  = os.path.dirname(os.path.abspath(__file__))
NC    = sys.argv[1] if len(sys.argv) > 1 else f"{HERE}/echo256_live_tt20.nc"
OUT   = sys.argv[2] if len(sys.argv) > 2 else f"{HERE}/fig_pileup_256_nocomp.png"
LABEL = sys.argv[3] if len(sys.argv) > 3 else r"echo256:  256$^2\times$64, nm=128, k$^4$ perp"

def load(f):
    v = Dataset(f).variables
    tt = np.asarray(v['tt'][:])
    sl = slice(len(tt)//2, len(tt))                    # steady 2nd half
    g2 = np.asarray(v['g2_bin'][:])[sl].mean(0)        # (kz,mm,kp)
    u  = np.asarray(v['upe2_bin'][:])[sl].mean(0)      # (kz,kp)  u_perp=|grad phi|^2
    b  = np.asarray(v['bpe2_bin'][:])[sl].mean(0)      # (kz,kp)  b_perp=|grad psi|^2
    kz = np.asarray(v['kz'][:]); kp = np.asarray(v['kpbin'][:])
    kxmax = float(np.max(np.abs(np.asarray(v['kx'][:]))))
    return dict(tt=tt, nrec=len(tt), twin=(float(tt[sl][0]), float(tt[sl][-1])),
                g_kp=g2.sum(axis=(0, 1)), g_kz=g2.sum(axis=(1, 2)),
                g_m=g2.sum(axis=(0, 2)),                       # E_g(m) Hermite spectrum
                u_kz=u.sum(axis=1), b_kz=b.sum(axis=1),
                u_kp=u.sum(axis=0), b_kp=b.sum(axis=0),
                kz=kz, kp=kp, kxmax=kxmax)

B = load(NC)
kp = B['kp']; kz = B['kz']
nm = len(B['g_m'])
kz_active = (B['u_kz'] + B['b_kz']) > 0
kp_cut = (2.0/3.0)*B['kxmax']
pp = (kp > 0) & (kp <= kp_cut); pz = (kz > 0) & kz_active

def hifrac(E, kk, kcut, thr=0.35):
    m = (kk > 0) & (kk <= kcut); kn = kk[m]/kk[m].max()
    return float(E[m][kn > thr].sum()/E[m].sum())
def mslope(E, m0=2, m1=16):
    mm = np.arange(len(E)); msk = (mm >= m0) & (mm <= m1) & (E > 0)
    p = np.polyfit(np.log(mm[msk]), np.log(E[msk]), 1)
    r = np.corrcoef(np.log(mm[msk]), np.log(E[msk]))[0, 1]**2
    return float(p[0]), float(r)

Cu = "#1f6feb"; Cb = "#e8890c"; Cg = "#2a9d5c"; GREY = '0.55'
def guide(ax, k0, k1, y0, s, txt):
    kk = np.array([k0, k1], float); y = y0*(kk/k0)**s
    ax.loglog(kk, y, color=GREY, lw=1, ls=':')
    ax.text(kk[1], y[1], txt, color=GREY, fontsize=9, va='center')

fig, ax = plt.subplots(3, 2, figsize=(14, 15))
fig.suptitle(LABEL + r"  —  pile-up diagnostic (raw spectra)"
             + f"\n(steady 2nd-half mean, tt={B['twin'][0]:.1f}-{B['twin'][1]:.1f}, {B['nrec']} recs)",
             fontsize=13.5, fontweight='bold')

# ================= ROW 1 LEFT : EM perp (k_perp) =================
a = ax[0, 0]
a.loglog(kp[pp], B['u_kp'][pp], color=Cu, lw=1.9, ls='-', label=r"$u_\perp$")
a.loglog(kp[pp], B['b_kp'][pp], color=Cb, lw=1.9, ls='-', label=r"$b_\perp$")
guide(a, 3, 20, B['u_kp'][pp][kp[pp] >= 3][0], -5/3, r"$k_\perp^{-5/3}$")
a.set_xlabel(r"$k_\perp$"); a.set_ylabel(r"$E(k_\perp)=\sum_{k_z}\frac{1}{2}|\cdot|^2$")
a.set_title(r"(a)  EM PERP $E(k_\perp)$", loc='left', fontsize=10.5, fontweight='bold')
a.legend(fontsize=9); a.grid(alpha=0.25, which='both')

# ================= ROW 1 RIGHT : EM parallel (k_z) =================
b = ax[0, 1]
b.loglog(kz[pz], B['u_kz'][pz], color=Cu, lw=1.9, ls='-', label=r"$u_\perp$")
b.loglog(kz[pz], B['b_kz'][pz], color=Cb, lw=1.9, ls='-', label=r"$b_\perp$")
guide(b, 2, kz[pz].max(), B['u_kz'][pz][kz[pz] >= 2][0], -2, r"$k_z^{-2}$")
b.set_xlabel(r"$k_z$"); b.set_ylabel(r"$E(k_z)=\sum_{k_\perp}\frac{1}{2}|\cdot|^2$")
b.set_title(r"(b)  EM PARALLEL $E(k_z)$", loc='left', fontsize=10.5, fontweight='bold')
b.legend(fontsize=9); b.grid(alpha=0.25, which='both')

# ================= ROW 2 LEFT : g perp (k_perp) — bottleneck check =================
thr = 0.35*kp[pp].max()
c = ax[1, 0]
c.axvspan(thr, kp[pp].max(), color='0.9', alpha=0.7, zorder=0)
c.loglog(kp[pp], B['g_kp'][pp], color=Cg, lw=1.9, ls='-', label="g")
guide(c, 3, 20, B['g_kp'][pp][kp[pp] >= 3][0], -5/3, r"$k_\perp^{-5/3}$")
fb = hifrac(B['g_kp'], kp, kp_cut)
c.set_xlabel(r"$k_\perp$"); c.set_ylabel(r"$E_g(k_\perp)=\sum_{k_z,m}\frac{1}{2}|g_m|^2$")
c.set_title(r"(c)  g PERP $E_g(k_\perp)$ (shaded: $k_\perp>0.35\,k_{\perp,max}$)", loc='left', fontsize=10.5, fontweight='bold')
c.text(0.02, 0.03, r"high-$k_\perp$ fraction ($>0.35\,k_{\perp,max}$): " + f"{100*fb:.1f}%",
       transform=c.transAxes, va='bottom', ha='left', fontsize=9,
       bbox=dict(boxstyle='round', fc='white', ec='0.7', alpha=0.9))
c.legend(fontsize=9, loc='upper right'); c.grid(alpha=0.25, which='both')

# ================= ROW 2 RIGHT : g parallel (k_z) =================
d = ax[1, 1]
d.loglog(kz[pz], B['g_kz'][pz], color=Cg, lw=1.9, ls='-', label="g")
guide(d, 2, kz[pz].max(), B['g_kz'][pz][kz[pz] >= 2][0], -2, r"$k_z^{-2}$")
d.set_xlabel(r"$k_z$"); d.set_ylabel(r"$E_g(k_z)=\sum_{k_\perp,m}\frac{1}{2}|g_m|^2$")
d.set_title(r"(d)  g PARALLEL $E_g(k_z)$  —  no $k_z$ sink by design (Meyrand eq.7)", loc='left', fontsize=10.5, fontweight='bold')
d.legend(fontsize=9); d.grid(alpha=0.25, which='both')

# ================= ROW 3 LEFT : g Hermite (m) =================
gm = B['g_m']; mm = np.arange(nm); mpos = mm[mm >= 1]
sl16 = mslope(gm, 2, 16); sl60 = mslope(gm, 2, min(60, nm-2))
e = ax[2, 0]
e.loglog(mpos, gm[mpos], color=Cg, lw=1.9, ls='-', marker='o', ms=3, label="g")
guide(e, 3, 40, gm[3], -0.5, r"$m^{-1/2}$")
e.set_xlabel(r"$m$  (Hermite moment)"); e.set_ylabel(r"$E_g(m)=\sum_{k_z,k_\perp}\frac{1}{2}|g_m|^2$")
e.set_ylim(bottom=max(gm[mpos].min()*0.5, 1e-3))
e.set_title(r"(e)  g HERMITE $E_g(m)$: all $m$, phase-mixing $m^{-1/2}$", loc='left', fontsize=10.5, fontweight='bold')
e.text(0.02, 0.03, f"inertial slope: $m$=2-16 {sl16[0]:+.2f} (R²{sl16[1]:.2f}), "
       + f"$m$=2-{min(60,nm-2)} {sl60[0]:+.2f}" + r"  ($\approx-1/2$)",
       transform=e.transAxes, va='bottom', ha='left', fontsize=8.5,
       bbox=dict(boxstyle='round', fc='white', ec='0.7', alpha=0.9))
e.legend(fontsize=9, loc='upper right'); e.grid(alpha=0.25, which='both')

# ================= ROW 3 RIGHT : blank =================
ax[2, 1].axis('off')

fig.tight_layout(rect=[0, 0, 1, 0.955])
fig.savefig(OUT, dpi=125)
print("saved", OUT)
print(f"nrec={B['nrec']} twin={B['twin'][0]:.2f}-{B['twin'][1]:.2f} nm={nm}")
print(f"g  kperp high-frac(>0.35kmax): {fb:.4f}")
print(f"g  m inertial slope m2-16={sl16[0]:+.3f}(R²{sl16[1]:.2f})  m2-{min(60,nm-2)}={sl60[0]:+.3f}")

#!/usr/bin/env python3
# ------------------------------------------------------------------
# 1-C-ii PRODUCTION echo analysis (128^2 x 64, nm=32).
#
# Deliverable: the time-averaged Hermite free-energy flux Gamma_m(k_perp)
# (eq 9, Meyrand 2019) should hover near 0 across the inertial range --
# the stochastic-echo / fluidization signature (forward phase-mixing
# largely cancelled by nonlinear anti-mixing).
#
# This is a richer, steady-window analysis than check_stage2.py:
#   - explicit steady window t > T_STEADY (default 10), not "second half"
#   - ECHO_RATIO reported per inertial k_perp (not just mean/min)
#   - steady-state confirmation from W_free(t) plateau (mean +/- rms)
#   - Hermite free-energy spectrum W_m(m) (steady-avg) if present
#   - telescoping top-slot Gamma_m_kint(nm-1) ~ 0 as a diagnostic check
#
# Usage:  analyze_prod128.py FILE.out.nc [T_STEADY]
# ------------------------------------------------------------------
import sys
import re
import subprocess
import numpy as np


def dump(nc, varlist):
    out = subprocess.run(['ncdump', '-p', '9,17', '-v', ','.join(varlist), nc],
                         capture_output=True, text=True, check=True).stdout
    data = out.split('data:', 1)[1]
    res = {}
    for v in varlist:
        m = re.search(r'(?m)^\s*' + re.escape(v) + r'\s*=\s*(.*?);', data, re.S)
        if not m:
            raise RuntimeError('variable %s not found in %s' % (v, nc))
        toks = m.group(1).replace('\n', ' ').replace(',', ' ').split()
        res[v] = np.array([float(x) for x in toks if x not in ('', '_')], float)
    return res


def have(nc, v):
    try:
        dump(nc, [v])
        return True
    except Exception:
        return False


nc = sys.argv[1]
T_STEADY = float(sys.argv[2]) if len(sys.argv) > 2 else 10.0

# ---- time base + steady state from W_free(t) ----
dW = dump(nc, ['W_free', 'tt'])
Wf, tt = dW['W_free'], dW['tt']
ntt = len(tt)
steady = np.where(tt >= T_STEADY)[0]
if steady.size < 3:                      # run too short: fall back to 2nd half
    steady = np.arange(max(1, ntt // 2), ntt)
    T_STEADY = float(tt[steady[0]])

print("NREC %d" % ntt)
print("T_REACHED %.4f" % tt[-1])
print("T_STEADY %.4f  (steady records: %d)" % (T_STEADY, steady.size))
Wsteady = Wf[steady]
print("W_FREE_FIRST %.6e" % Wf[0])
print("W_FREE_LAST %.6e" % Wf[-1])
print("W_FREE_STEADY_MEAN %.6e" % float(np.mean(Wsteady)))
print("W_FREE_STEADY_RMS %.6e" % float(np.std(Wsteady)))
print("W_FREE_STEADY_RELRMS %.4f" %
      (float(np.std(Wsteady)) / max(float(np.mean(Wsteady)), 1e-300)))

# ---- telescoping (binning-immune top slot) ----
if have(nc, 'Gamma_m_kint'):
    dk = dump(nc, ['Gamma_m_kint', 'mm'])
    nm = len(dk['mm'])
    Gk = dk['Gamma_m_kint'].reshape(ntt, nm)
    top = np.abs(Gk[steady, -1])
    peak = max(float(np.max(np.abs(Gk[steady]))), 1e-300)
    print("TELE_KINT %.6e" % (float(np.max(top)) / peak))

# ---- Hermite free-energy spectrum W_m(m), steady-avg ----
if have(nc, 'W_m'):
    dwm = dump(nc, ['W_m', 'mm'])
    nm = len(dwm['mm'])
    Wm = dwm['W_m'].reshape(ntt, nm)
    Wm_s = Wm[steady].mean(axis=0)
    print("WM_STEADY " + " ".join("%.3e" % v for v in Wm_s))

# ---- Gamma_m(k_perp): the echo deliverable ----
d = dump(nc, ['Gamma_m', 'kpbin', 'kz', 'mm'])
nkp, nkz, nm = len(d['kpbin']), len(d['kz']), len(d['mm'])
kpbin = d['kpbin']
G = d['Gamma_m'].reshape(ntt, nkz, nm, nkp)
# steady-avg over time -> (nkz, nm, nkp), then sum over kz -> (nm, nkp)
Gavg = G[steady].mean(axis=0).sum(axis=0)

kp_max = float(kpbin[-1])
inertial = np.where((kpbin > 2.0) & (kpbin < kp_max / 3.0))[0]
if inertial.size == 0:
    inertial = np.arange(2, max(3, nkp // 3))

print("INERTIAL_BINS %d  (k_perp %.2f..%.2f)" %
      (inertial.size, float(kpbin[inertial[0]]), float(kpbin[inertial[-1]])))

ratios = []
print("ECHO_BY_KPERP  (k_perp  echo_ratio=max|Gamma_m|/sum|src_m|)")
for ik in inertial:
    cum = Gavg[:, ik]
    src = np.diff(cum, prepend=0.0)
    gross = float(np.sum(np.abs(src)))
    net = float(np.max(np.abs(cum)))
    if gross > 1e-300:
        r = net / gross
        ratios.append(r)
        print("  %6.2f  %.4f" % (float(kpbin[ik]), r))
ratios = np.array(ratios)
if ratios.size:
    print("ECHO_RATIO_MEAN %.6e" % float(np.mean(ratios)))
    print("ECHO_RATIO_MEDIAN %.6e" % float(np.median(ratios)))
    print("ECHO_RATIO_MIN %.6e" % float(np.min(ratios)))
    print("ECHO_RATIO_MAX %.6e" % float(np.max(ratios)))

# representative mid-inertial k_perp: cumulative Gamma_m(m) profile
if inertial.size:
    ikmid = inertial[inertial.size // 2]
    prof = Gavg[:, ikmid]
    scale = max(float(np.max(np.abs(prof))), 1e-300)
    print("PROFILE_KPERP %.2f" % float(kpbin[ikmid]))
    print("PROFILE_GAMMA_M " + " ".join("%.3e" % v for v in prof))
    print("PROFILE_NORM " + " ".join("%.3f" % (v / scale) for v in prof))

#!/usr/bin/env python3
# ------------------------------------------------------------------
# 1-C-ii Stage 2 SMOKE analysis: does the driven-turbulence pipeline
# run end-to-end, does the Gamma_m diagnostic stay self-consistent in
# strong turbulence, and does the stochastic-echo signature appear?
#
#   PIPELINE   : run produced evolving turbulence (W_free grows then
#                saturates); report N records and reached time t.
#   TELE_KINT  : k-integrated top slot Gamma(nm-1)/peak over evolved
#                records -- must stay ~machine zero even in strong
#                turbulence (validates the diagnostic in the target
#                regime, binning-immune, same invariant as Stage 0).
#   ECHO_RATIO : the echo metric. At each inertial-range k_perp the
#                CUMULATIVE flux Gamma_m(k_perp) is compared to the
#                GROSS per-rung phase-mixing activity sum_m|src_m|:
#                  src_m(k) = Gamma_m(k) - Gamma_{m-1}(k)   (de-cumulate)
#                  ECHO_RATIO(k) = max_m|Gamma_m(k)| / sum_m|src_m(k)|
#                A pure forward Hermite cascade (Stage 1 linear) has all
#                src_m of one sign => cumulative reaches the full sum =>
#                ratio ~ 1. The stochastic echo cancels forward mixing
#                with anti-mixing => small net flux => ratio << 1.
#                Time-averaged over the steady (second-half) records.
# ------------------------------------------------------------------
import sys, re, subprocess
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


def load_gamma(nc):
    d = dump(nc, ['Gamma_m', 'kpbin', 'kz', 'mm', 'tt'])
    nkp, nkz, nm = len(d['kpbin']), len(d['kz']), len(d['mm'])
    ntt = len(d['tt'])
    G = d['Gamma_m'].reshape(ntt, nkz, nm, nkp)
    return G, d['kpbin'], d['tt']


def load_kint(nc):
    d = dump(nc, ['Gamma_m_kint', 'mm', 'tt'])
    nm, ntt = len(d['mm']), len(d['tt'])
    return d['Gamma_m_kint'].reshape(ntt, nm), d['tt']


nc = sys.argv[1]

# ---- PIPELINE: turbulence developed? ----
dW = dump(nc, ['W_free', 'tt'])
Wf, tt = dW['W_free'], dW['tt']
ntt = len(tt)
print("NREC %d" % ntt)
print("T_REACHED %.4f" % tt[-1])
print("W_FREE_FIRST %.6e" % Wf[0])
print("W_FREE_LAST %.6e" % Wf[-1])
print("W_FREE_MAX %.6e" % float(np.max(Wf)))

# ---- TELE_KINT: telescoping holds in strong turbulence? ----
Gk, _ = load_kint(nc)
evolved = np.arange(1, ntt)
top = np.abs(Gk[evolved, -1])
peak = max(float(np.max(np.abs(Gk[evolved]))), 1e-300)
print("TELE_KINT %.6e" % (float(np.max(top)) / peak))

# ---- ECHO_RATIO: net cumulative flux vs gross phase-mixing ----
G, kpbin, _ = load_gamma(nc)          # (ntt, nkz, nm, nkp)
nkp, nm = G.shape[3], G.shape[2]
# steady window = second half of the record (skip spin-up transient)
half = np.arange(max(1, ntt // 2), ntt)
Gavg = G[half].mean(axis=0).sum(axis=0)          # (nm, nkp): <Gamma_m(k_perp)>, sum over kz
# inertial band: above the forcing scale (k_perp>~2), below dissipation (< kmax/3)
kp_max = float(kpbin[-1])
inertial = np.where((kpbin > 2.0) & (kpbin < kp_max / 3.0))[0]
if inertial.size == 0:
    inertial = np.arange(2, max(3, nkp // 3))
ratios = []
for ik in inertial:
    cum = Gavg[:, ik]                             # cumulative flux Gamma_m at this k_perp
    src = np.diff(cum, prepend=0.0)              # per-rung source
    gross = float(np.sum(np.abs(src)))
    net = float(np.max(np.abs(cum)))
    if gross > 1e-300:
        ratios.append(net / gross)
ratios = np.array(ratios)
print("INERTIAL_BINS %d  (k_perp %.2f..%.2f)" %
      (inertial.size, float(kpbin[inertial[0]]), float(kpbin[inertial[-1]])))
print("ECHO_RATIO_MEAN %.6e" % (float(np.mean(ratios)) if ratios.size else -1))
print("ECHO_RATIO_MIN %.6e" % (float(np.min(ratios)) if ratios.size else -1))

# representative mid-inertial k_perp: print the cumulative Gamma_m(m) profile
# (echo => hovers near 0 across intermediate m; forward cascade => ramps up).
if inertial.size:
    ikmid = inertial[inertial.size // 2]
    prof = Gavg[:, ikmid]
    scale = max(float(np.max(np.abs(prof))), 1e-300)
    print("PROFILE_KPERP %.2f" % float(kpbin[ikmid]))
    print("PROFILE_GAMMA_M " + " ".join("%.3e" % v for v in prof))
    print("PROFILE_NORM " + " ".join("%.3f" % (v / scale) for v in prof))

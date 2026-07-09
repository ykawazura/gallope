#!/usr/bin/env python3
# ------------------------------------------------------------------
# Milestone 1-C-ii verification of the Hermite free-energy flux
# Gamma_m(k_perp) (Meyrand 2019 eq 9). netCDF4-python is absent on
# Miyabi, so we read the classic-format output via ncdump -p 9,17.
#
# TELESCOPING -- the key subtlety.
#   The per-rung source src_m(k) sums to zero over m ONLY after the
#   k-sum, and the two streaming pieces behave differently:
#     * LINEAR zi*kz*S_m is k-diagonal  -> sum_m src_m(k) = 0 at EVERY
#       fixed k (weighted antisymmetry w_m cp_m = w_{m+1} cm_{m+1} plus
#       ghost g_nm = 0). So the binned top-m slot is ~0 per k_perp.
#     * The {Psi,S_m} bracket transfers free energy ACROSS k_perp (a
#       perpendicular cascade), so sum_m src_m(k) != 0 at fixed k --
#       it cancels only after the full k-sum (int a{Psi,b}=-int b{Psi,a}).
#   Therefore the binned top-m slot is a valid telescoping check for the
#   LINEAR run but NOT for the nonlinear bracket. The k-integrated flux
#   Gamma(m) = sum_k src (Gamma_m_kint, all modes, no polar-corner drop)
#   carries the correct invariant: its top slot Gamma(nm-1) = sum_{m,k}
#   src_m ~ 0 (= -dW/dt at mu=nu=0) for BOTH linear and nonlinear.
#
#   TELE_LIN      : s0_lin binned top-m slot / peak, per k_perp. Isolates
#                   the Hermite coeff/weight/sign logic (pure spectral
#                   arithmetic, no FFT/aliasing). Expect ~1e-13.
#   TELE_LIN_KINT : s0_lin k-integrated top slot / peak. Sanity ~1e-13.
#   TELE_NL_KINT  : s0_nl_pm1 k-integrated top slot / peak over EVOLVED
#                   records (t>0; the t=0 raw broadband IC is reported
#                   separately). This is the binning-immune telescoping
#                   test of the {Psi,S_m} bracket. Expect ~1e-12.
#   BUDGET_NL     : max_evolved |dW_free/dt + Gamma_kint(nm-1)| (both ~0
#                   at mu=nu=0); informational closure of dW/dt=-Gamma_top.
#   PMREL         : max|Gamma_pm1 - Gamma_pm4| / max|Gamma| (split assembly).
#   ORDERING      : reshape self-check via stage1 single_mode (energy at
#                   k_perp=0 only).
#   SIGN_PEAK     : stage1 linear phase-mixing -> forward Hermite cascade,
#                   signed Gamma_m at the peak-|Gamma| box must be > 0.
# ------------------------------------------------------------------
import sys, re, subprocess
import numpy as np


def dump(nc, varlist):
    """Return {var: 1D float ndarray} read from ncdump's data section."""
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
    """Return (G, tt) with G shaped (ntt, nkz, nm, nkpolar).

    Gamma_m NetCDF dims are Fortran (kpbin, mm, kz, tt); ncdump prints
    C-order (tt, kz, mm, kpbin), so per record the fastest axis is kpbin.
    """
    d = dump(nc, ['Gamma_m', 'kpbin', 'kz', 'mm', 'tt'])
    nkp, nkz, nm = len(d['kpbin']), len(d['kz']), len(d['mm'])
    tt = d['tt']
    ntt = len(tt)
    rec = nkz * nm * nkp
    assert len(d['Gamma_m']) == ntt * rec, \
        'Gamma_m size %d != ntt*nkz*nm*nkpolar %d' % (len(d['Gamma_m']), ntt * rec)
    G = d['Gamma_m'].reshape(ntt, nkz, nm, nkp)
    return G, tt


def load_kint(nc):
    """Return (Gk, tt) with Gk shaped (ntt, nm).

    Gamma_m_kint NetCDF dims are Fortran (mm, tt); ncdump prints C-order
    (tt, mm), fastest axis mm.
    """
    d = dump(nc, ['Gamma_m_kint', 'mm', 'tt'])
    nm, ntt = len(d['mm']), len(d['tt'])
    assert len(d['Gamma_m_kint']) == ntt * nm
    return d['Gamma_m_kint'].reshape(ntt, nm), d['tt']


def tele(G, recs):
    """Binned top-m cumulative residual / peak |Gamma| over given records."""
    sub = G[recs]
    top = np.abs(sub[:, :, -1, :])              # cumulative over ALL m => ~0 (per k)
    peak = float(np.max(np.abs(sub)))
    return float(np.max(top)) / max(peak, 1e-300)


def tele_kint(Gk, recs):
    """k-integrated top-slot residual / peak |Gamma(m)| over given records."""
    sub = Gk[recs]
    top = np.abs(sub[:, -1])                     # Gamma(nm-1) = sum_{m,k} src ~ 0
    peak = float(np.max(np.abs(sub)))
    return float(np.max(top)) / max(peak, 1e-300)


s0_lin, s0_nl1, s0_nl4, stage1 = sys.argv[1:5]

# ---- TELE_LIN : linear telescoping, binned per-k (coeff/weight/sign) ----
Gl, tl = load_gamma(s0_lin)
Glk, _ = load_kint(s0_lin)
allrec = np.arange(len(tl))
print("NREC_LIN %d" % len(tl))
print("TELE_LIN %.6e" % tele(Gl, allrec))
print("TELE_LIN_KINT %.6e" % tele_kint(Glk, allrec))

# ---- TELE_NL_KINT : bracket telescoping, k-integrated, EVOLVED records ----
Gn, tn = load_gamma(s0_nl1)
Gnk, _ = load_kint(s0_nl1)
dn = dump(s0_nl1, ['W_free', 'tt'])
Wf = dn['W_free']
evolved = np.arange(1, len(tn))
print("NREC_NL %d" % len(tn))
print("TELE_NL_KINT_T0 %.6e" % tele_kint(Gnk, np.array([0])))
print("TELE_NL_KINT %.6e" % tele_kint(Gnk, evolved))
# reference-only: the binned nonlinear top slot is O(1) BY DESIGN (perp
# transfer), reported so the contrast with TELE_NL_KINT is on the record.
print("TELE_NL_BINNED_REF %.6e" % tele(Gn, evolved))

# ---- BUDGET_NL : informational closure dW/dt = -Gamma_kint(nm-1) --------
# Skip rec 0 and 1 (the raw-IC 2/3-filter transient contaminates the finite
# difference); over the remaining interior both sides are ~0 at mu=nu=0.
if len(tn) >= 5:
    dWdt = np.gradient(Wf, tn)
    gtop = Gnk[:, -1]
    ii = np.arange(2, len(tn) - 1)
    resid = np.abs(dWdt[ii] + gtop[ii])
    scale = max(float(np.max(np.abs(Gnk))), 1e-300)
    print("BUDGET_NL %.6e" % (float(np.max(resid)) / scale))
else:
    print("BUDGET_NL -1")

# ---- PMREL : Gamma_m P_m-invariance (pm1 vs pm4) --------------------
Gn4, tn4 = load_gamma(s0_nl4)
k = min(Gn.shape[0], Gn4.shape[0])
absd = np.abs(Gn[:k] - Gn4[:k])
peak = max(float(np.max(np.abs(Gn[:k]))), 1e-300)
pmrel = float(np.max(absd)) / peak
print("NREC_PM1 %d NREC_PM4 %d" % (Gn.shape[0], Gn4.shape[0]))
print("PMREL %.6e" % pmrel)

# ---- ORDERING + SIGN : stage1 (single_mode, linear forward cascade) --
d = dump(stage1, ['g2_bin', 'kpbin', 'kz', 'mm', 'tt'])
nkp, nkz, nm = len(d['kpbin']), len(d['kz']), len(d['mm'])
ntt = len(d['tt'])
g2 = d['g2_bin'].reshape(ntt, nkz, nm, nkp)
e_all = float(np.sum(np.abs(g2)))
e_off = float(np.sum(np.abs(g2[:, :, :, 1:])))
ordering = e_off / max(e_all, 1e-300)
print("ORDERING %.6e" % ordering)

Gs, ts = load_gamma(stage1)
# The whole run stays in the forward-mixing band (t < recurrence null
# ~2*sqrt(2*nm)). Gamma_m(k_perp=0) over intermediate m and all records:
# read the SIGNED value at the peak-|Gamma| box -- the dominant forward
# flux, which a sign error would flip negative.
gam0 = Gs[:, :, 1:nm - 1, 0].sum(axis=1)        # (ntt, nm-2): sum over kz
ridx, midx = np.unravel_index(int(np.argmax(np.abs(gam0))), gam0.shape)
print("SIGN_TREC %.4f" % ts[ridx])
print("SIGN_PEAK %.6e" % float(gam0[ridx, midx]))
print("SIGN_ABSPEAK %.6e" % float(np.max(np.abs(gam0))))
sig = np.abs(gam0) > 0.05 * np.max(np.abs(gam0))
print("SIGN_MIN_SIG %.6e" % float(np.min(gam0[sig])))

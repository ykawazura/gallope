#!/usr/bin/env python3
# ------------------------------------------------------------------
# Milestone 1-C-i verification of the new g NetCDF diagnostics.
# netCDF4 python is absent on Miyabi, so we read the classic-format
# output through ncdump -p 9,17 (full double precision round-trip).
#
#   DRIFT     : max_t |W_free(t) - W_free(0)| / W_free(0)
#               collisionless single-mode run -> eq(6) free energy is
#               conserved, so this must stay ~0 (1-B saw ~2e-9).
#   INTERNAL  : max_t |sum_bins g2_bin - sum_m W_m| / sum_m W_m
#               the single-mode IC keeps ALL energy at kperp=0 (bin 1),
#               so the binned spectrum must reproduce the independently
#               reduced per-m free energy exactly (no top-bin tail).
#   WFREE_REL : max rel diff of W_free between P_m=1 and P_m=4 runs
#   G2_MAXREL : max rel diff of g2_bin between P_m=1 and P_m=4 runs
#               g is bit-identical across the split (1-B); only the MPI
#               reduction order differs -> ULP-level, must be < 1e-9.
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


landau, pm1, pm4 = sys.argv[1], sys.argv[2], sys.argv[3]

# ---- (i) drift + (ii) internal consistency (collisionless single mode) ----
d = dump(landau, ['W_free', 'W_m', 'g2_bin'])
Wf, Wm, g2 = d['W_free'], d['W_m'], d['g2_bin']
ntt = len(Wf)
drift = float(np.max(np.abs((Wf - Wf[0]) / Wf[0])))
nm = len(Wm) // ntt
sumWm = Wm.reshape(ntt, nm).sum(axis=1)
sumbin = g2.reshape(ntt, -1).sum(axis=1)
good = sumWm > 0
internal = float(np.max(np.abs(sumbin[good] - sumWm[good]) / sumWm[good]))
print("NREC_LANDAU %d" % ntt)
print("DRIFT %.6e" % drift)
print("INTERNAL %.6e" % internal)

# ---- (iii) P_m invariance of W_free and g2_bin ----
a = dump(pm1, ['W_free', 'g2_bin'])
b = dump(pm4, ['W_free', 'g2_bin'])
na, nb = len(a['W_free']), len(b['W_free'])
print("NREC_PM1 %d NREC_PM4 %d" % (na, nb))
k = min(na, nb)
wa, wb = a['W_free'][:k], b['W_free'][:k]
wfree_rel = float(np.max(np.abs(wa - wb) / np.maximum(np.abs(wb), 1e-300)))
chunk = len(a['g2_bin']) // na
assert chunk == len(b['g2_bin']) // nb, "g2_bin record size mismatch"
ga, gb = a['g2_bin'][:k * chunk], b['g2_bin'][:k * chunk]
absd = np.abs(ga - gb)
sel = np.abs(gb) > 1e-14
maxrel = float(np.max(absd[sel] / np.abs(gb[sel]))) if np.any(sel) else 0.0
print("WFREE_REL %.6e" % wfree_rel)
print("G2_MAXREL %.6e" % maxrel)
print("G2_MAXABS %.6e" % float(np.max(absd) if len(absd) else 0.0))

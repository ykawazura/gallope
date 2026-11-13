#!/usr/bin/env python3
"""
Evaluate the RRMHD validation runs created by tests/RRMHD/setup_runs.sh.

  python3 tests/RRMHD/check_runs.py <runs dir>          # main set
  python3 tests/RRMHD/check_runs.py <runs dir> linear   # V7

Needs numpy, scipy and netCDF4. Reads <case>/gallope.out.nc, <case>/gallope.in
and <case>/restart/*.dat (raw complex128, Fortran order (nkz, nky, nkx)).
"""
import os
import re
import sys

import numpy as np
from netCDF4 import Dataset
from scipy.linalg import expm

RUNS = sys.argv[1]
SET = sys.argv[2] if len(sys.argv) > 2 else 'main'
NLX, NLY, NLZ = 64, 64, 32
NKX, NKY, NKZ = NLX, NLY, NLZ//2 + 1
FIELDS = ['phi', 'omg', 'psi', 'upa', 'bpa']
ENERGY = ['upe2_sum', 'bpe2_sum', 'upa2_sum', 'bpa2_sum']
DOT = [e.replace('2_sum', '2dot_sum') for e in ENERGY]
DISSIP = [e.replace('2_sum', '2dissip_sum') for e in ENERGY]
SOURCES = ['p_aw_rot_sum', 'p_compr_rot_sum', 'p_aw_grd_sum', 'p_compr_grd_sum']
ELS = ['zpep2_sum', 'zpem2_sum', 'zpap2_sum', 'zpam2_sum']
ok = True


def report(name, passed, msg):
    global ok
    ok &= passed
    print(f"[{'PASS' if passed else 'FAIL'}] {name}: {msg}")


def series(case):
    with Dataset(os.path.join(RUNS, case, 'gallope.out.nc')) as nc:
        out = {v: np.array(nc.variables[v][:]) for v in nc.variables
               if nc.variables[v].dimensions == ('tt',)}
    return out


def params(case):
    with Dataset(os.path.join(RUNS, case, 'gallope.out.nc')) as nc:
        return {v: float(nc.variables[v][:]) for v in
                ['beta', 'gamma', 'q', 'kappa_b', 'kappa_p']}


def inputs(case, names):
    text = open(os.path.join(RUNS, case, 'gallope.in')).read()
    out = {}
    for n in names:
        m = re.search(r'^\s*' + n + r'\s*=\s*([-+0-9.dDeE]+)', text, re.M)
        out[n] = float(m.group(1).replace('d', 'e').replace('D', 'e'))
    return out


PARAMS = ['beta', 'gamma', 'q', 'kappa_b', 'kappa_p'] + \
    [f'{a}{b}_{c}{e}' for a in ['nu', 'eta'] for b in ['pe', 'pa'] for c in 'xz' for e in ['', '_exp']]
# base case of tests/RRMHD/gallope.in
BASE = dict(beta=1.0, gamma=1.666666666666667, q=1.5, kappa_b=0.2, kappa_p=0.3,
            **{f'{a}{b}_x': 10.0 for a in ['nu', 'eta'] for b in ['pe', 'pa']},
            **{f'{a}{b}_z': 1.0 for a in ['nu', 'eta'] for b in ['pe', 'pa']},
            **{f'{a}{b}_{c}_exp': 4.0 for a in ['nu', 'eta'] for b in ['pe', 'pa'] for c in 'xz'})


def check_params(cases):
    """The parameters recorded in gallope.out.nc must equal those in the
    case's gallope.in (a ';' in a namelist line silently ends the group with
    nvfortran), and the base input must hold the intended values."""
    base = inputs_all(os.path.join('tests', 'RRMHD', 'gallope.in')) \
        if os.path.exists(os.path.join('tests', 'RRMHD', 'gallope.in')) else None
    if base is not None:
        bad = [k for k in PARAMS if base[k] != BASE[k]]
        report('params tests/RRMHD/gallope.in', not bad, 'intended base values' if not bad else f'differ: {bad}')
    for c in cases:
        with Dataset(os.path.join(RUNS, c, 'gallope.out.nc')) as nc:
            rec = {k: float(nc.variables[k][:]) for k in PARAMS}
        want = inputs_all(os.path.join(RUNS, c, 'gallope.in'))
        bad = [f'{k}={rec[k]} (in: {want[k]})' for k in PARAMS if rec[k] != want[k]]
        report(f'params {c}', not bad, 'gallope.out.nc = gallope.in' if not bad else f'differ: {bad}')


def inputs_all(path):
    text = open(path).read()
    out = {}
    for n in PARAMS:
        m = re.findall(r'^\s*' + n + r'\s*=\s*([-+0-9.dDeE]+)\s*(?:!.*)?$', text, re.M)
        assert len(m) == 1, f'{path}: {n} must appear once, alone on its line'
        out[n] = float(m[0].replace('d', 'e').replace('D', 'e'))
    return out


def field(case, name):
    a = np.fromfile(os.path.join(RUNS, case, 'restart', name + '.dat'), dtype=np.complex128)
    return a.reshape(NKX, NKY, NKZ)   # [i, j, k]


def restart_time(case):
    with open(os.path.join(RUNS, case, 'restart', 'time.dat')) as f:
        return float(f.readlines()[1].split()[0])


def same_bytes(a, b, names):
    diff = [n for n in names
            if open(os.path.join(RUNS, a, 'restart', n), 'rb').read()
            != open(os.path.join(RUNS, b, 'restart', n), 'rb').read()]
    return diff


def max_rel_series(a, b, keys):
    sa, sb = series(a), series(b)
    n = min(len(sa['tt']), len(sb['tt']))
    worst = 0.0
    for k in keys:
        scale = np.max(np.abs(sa[k][:n])) or 1.0
        worst = max(worst, np.max(np.abs(sa[k][:n] - sb[k][:n]))/scale)
    return worst, n


def max_rel_fields(a, b, names):
    worst = 0.0
    for n in names:
        fa, fb = field(a, n), field(b, n)
        scale = np.max(np.abs(fa)) or 1.0
        worst = max(worst, np.max(np.abs(fa - fb))/scale)
    return worst


def budget(case, tol):
    """d/dt E = p_rot + p_grd - dissipation, sample by sample (t > 0)."""
    s = series(case)
    dot = sum(s[k] for k in DOT)[1:]
    src = sum(s[k] for k in SOURCES)[1:]
    rhs = (src - sum(s[k] for k in DISSIP)[1:])
    err = np.max(np.abs(dot - rhs))/np.max(np.abs(dot))
    report(f'budget {case}', err < tol,
           f'max|dE/dt - (p_rot + p_grd - dissip)|/max|dE/dt| = {err:.2e}, '
           f'max|p_rot + p_grd|/max|dE/dt| = {np.max(np.abs(src))/np.max(np.abs(dot)):.2e}')


def main_set():
    check_params([c for c in sorted(os.listdir(RUNS)) if c.startswith('rr_')
                  and not c.startswith('rr_lin') and c != 'rr_driven'])

    # V2 smoke: series finite, final state finite
    for c in sorted(os.listdir(RUNS)):
        if not c.startswith('rr_') or c.startswith('rr_lin') or c == 'rr_driven':
            continue
        s = series(c)
        finite = all(np.all(np.isfinite(v)) for v in s.values())
        report(f'V2 {c}', finite, f"{len(s['tt'])} samples, t_end = {s['tt'][-1]:.4f}")

    # V3 decomposition: 1 GPU vs 4 GPUs, same OT3
    keys = ENERGY + DISSIP + SOURCES + ELS
    d, n = max_rel_series('rr_ot3_np1', 'rr_ot3_np4', keys)
    report('V3 series np1 vs np4', d < 1e-12, f'max rel diff {d:.2e} over {n} samples')
    d = max_rel_fields('rr_ot3_np1', 'rr_ot3_np4', FIELDS)
    report('V3 restart np1 vs np4', d < 1e-12, f'max rel diff {d:.2e}')

    # V4 redundancy: P_m = 2 must reproduce P_m = 1 bit for bit
    files = [f + '.dat' for f in FIELDS]
    diff = same_bytes('rr_ot3_np2', 'rr_ot3_pm2', files)
    report('V4 restart np2 vs pm2 (bitwise)', not diff, 'identical' if not diff else f'differ: {diff}')
    pm2 = open(os.path.join(RUNS, 'rr_ot3_pm2', 'out.std')).read()
    report('V4 check_redundant_consistency', 'FAIL' not in pm2 and 'PASS' in pm2,
           f"{pm2.count('PASS')} PASS lines")

    # V5 restart: 100 + 100 steps (final state in rr_rs1/restart) vs 200 steps
    diff = same_bytes('rr_rs1', 'rr_ot3_np1', files)
    t1, t2 = restart_time('rr_rs1'), restart_time('rr_ot3_np1')
    report('V5 restart 100+100 vs 200 (bitwise)', not diff and t1 == t2,
           f'tt {t1} vs {t2}; ' + ('identical' if not diff else f'differ: {diff}'))

    # V6 energy conservation (q = kappa = 0, ideal): drift at dt and dt/2
    drift = []
    for c in ['rr_ideal_dt1', 'rr_ideal_dt2']:
        s = series(c)
        e = sum(s[k] for k in ENERGY)
        drift.append(abs(e[-1] - e[0])/abs(e[0]))
    r = drift[0]/drift[1] if drift[1] > 0 else np.inf
    report('V6 E drift ratio', 6 < r < 10 or max(drift) < 1e-13,
           f'drift dt={drift[0]:.3e}, dt/2={drift[1]:.3e}, ratio {r:.2f} (RK3: 8)')
    # with q = kappa = 0 the Coriolis exchange p_aw_rot = -p_compr_rot does no
    # net work and the gradient sources vanish
    s = series('rr_ideal_dt2')
    net = np.max(np.abs(s['p_aw_rot_sum'] + s['p_compr_rot_sum']))/np.max(np.abs(s['p_aw_rot_sum']))
    grd = np.max(np.abs(s['p_aw_grd_sum'])) + np.max(np.abs(s['p_compr_grd_sum']))
    report('V6 sources (q = kappa = 0)', net < 1e-12 and grd == 0.0,
           f'|p_aw_rot + p_compr_rot|/|p_aw_rot| = {net:.1e}, |p_grd| = {grd:.1e}')

    # Energy budget with rotation and gradients, eSSPIFRK3 and gear3
    budget('rr_ot3_np1', 1e-2)
    budget('rr_gear3', 1e-2)

    # Elsasser energies (current fields, no 1/2, as in Calliope) are state
    # quantities: zpep2 + zpem2 = 4 (upe2 + bpe2) and zpap2 + zpam2 =
    # 4 (upa2 + bpa2) (bpa2 includes A) hold to round-off, summed and per bin
    for c in ('rr_ot3_np1', 'rr_ot3_np4'):
        s = series(c)
        with Dataset(os.path.join(RUNS, c, 'gallope.out.nc')) as nc:
            b = {v: np.array(nc.variables[v][:]) for v in nc.variables if v.endswith('_bin')}
        errs = []
        for tag, d in (('sum', s), ('bin', b)):
            for zp, zm, u, w in (('zpep2', 'zpem2', 'upe2', 'bpe2'),
                                 ('zpap2', 'zpam2', 'upa2', 'bpa2')):
                rhs = 4*(d[f'{u}_{tag}'] + d[f'{w}_{tag}'])
                errs.append(np.max(np.abs(d[f'{zp}_{tag}'] + d[f'{zm}_{tag}'] - rhs))
                            / np.max(np.abs(rhs)))
        report(f'Elsasser identities {c}', max(errs) < 1e-12,
               'max|z+^2 + z-^2 - 4 (u^2 + b^2)|/max(4 (u^2 + b^2)): '
               f'perp sum {errs[0]:.1e}, par sum {errs[1]:.1e}, '
               f'perp bin {errs[2]:.1e}, par bin {errs[3]:.1e}')

    # forcing is refused
    out = open(os.path.join(RUNS, 'rr_driven', 'out.std')).read()
    done = os.path.exists(os.path.join(RUNS, 'rr_driven', 'restart', 'time.dat'))
    report('driven = T stops', 'RRMHD has no forcing' in out and not done,
           f"message {'found' if 'RRMHD has no forcing' in out else 'missing'}, "
           f"final restart {'written' if done else 'absent'}")


def operator(kx, ky, kz, p):
    """Real R with d/dt (phi, psi, upa, bpa) = i R (phi, psi, upa, bpa), from
    advance::get_ext_terms (omg = -K phi)."""
    K = kx*kx + ky*ky
    A = 2/(p['beta']*p['gamma']) + 1
    c = p['beta']*p['gamma']/2
    q, kb, kp = p['q'], p['kappa_b'], p['kappa_p']
    return np.array([[0, kz, 2*ky/K, -2*(kb + c*kp)*ky/K],
                     [kz, 0, 0, 0],
                     [(2 - q)*ky, -c*kp*ky, 0, kz],
                     [-(2*kb + (c - 1)*kp)*ky/A, q*ky/A, kz/A, 0]])


def quartic_roots(kx, ky, kz, p):
    """Roots of the dispersion relation derived from the Calliope equations
    (exp(-i omega t)); a1 uses A = 1 + 1/c, so it is independent of operator()."""
    K = kx*kx + ky*ky
    A = 2/(p['beta']*p['gamma']) + 1
    c = p['beta']*p['gamma']/2
    q, kb, kp = p['q'], p['kappa_b'], p['kappa_p']
    a2 = -(A*K*kz**2 + K*kz**2 + (4 - 2*q)*A*ky**2 + 2*c**2*kp**2*ky**2 + 6*c*kb*kp*ky**2
           - 2*c*kp**2*ky**2 + 4*kb**2*ky**2 - 2*kb*kp*ky**2)/(A*K)
    a1 = -8*ky**2*kz*(kb + c*kp)/(A*K)
    a0 = kz**2*(K*kz**2 - 2*q*ky**2 - 2*c**2*kp**2*ky**2 - 2*c*kb*kp*ky**2)/(A*K)
    return np.roots([1, 0, a2, a1, a0])


def sort_c(z):
    return np.array(sorted(z, key=lambda x: (round(x.real, 9), x.imag)))


def linear_set():
    LIN = ['rr_lin_ap', 'rr_lin_am', 'rr_lin_rot_a', 'rr_lin_rot_b', 'rr_lin_unst']
    check_params(LIN)
    for c in LIN:
        p = params(c)
        inp = inputs(c, ['mode_kx', 'mode_ky', 'mode_kz', 'amp_phi', 'amp_psi', 'amp_upa', 'amp_bpa'])
        mx, my, mz = (int(inp[k]) for k in ['mode_kx', 'mode_ky', 'mode_kz'])
        x0 = np.array([inp[k] for k in ['amp_phi', 'amp_psi', 'amp_upa', 'amp_bpa']])
        R = operator(mx, my, mz, p)
        w = np.linalg.eigvals(R)

        # the operator transcribed here against the dispersion relation of the spec
        spec = np.max(np.abs(sort_c(-w.astype(complex)) - sort_c(quartic_roots(mx, my, mz, p))))

        t = restart_time(c)
        idx = (mx % NKX, my % NKY, mz)
        f = np.array([field(c, n)[idx] for n in ['phi', 'psi', 'upa', 'bpa']])
        pred = expm(1j*R*t) @ x0
        err = np.max(np.abs(f - pred))/np.max(np.abs(pred))
        K = mx*mx + my*my
        omg_err = abs(field(c, 'omg')[idx] + K*f[0])/abs(K*f[0])
        leak = 0.0
        for n in FIELDS:
            g = field(c, n).copy()
            g[idx] = 0.0
            leak = max(leak, np.max(np.abs(g)))

        # frequency of each component, f = x0 exp(i w t) for an eigenvector
        om = np.angle(f/x0)/t
        if c in ('rr_lin_ap', 'rr_lin_am'):
            sgn = 1 if c == 'rr_lin_ap' else -1
            A = 2/(p['beta']*p['gamma']) + 1
            om_th = sgn*mz*np.array([1, 1, 1/np.sqrt(A), 1/np.sqrt(A)])
            what = 'Alfven/slow'
        elif c.startswith('rr_lin_rot'):
            res = np.abs(R @ x0[:, None] - x0[:, None]*w[None, :]).max(axis=0)
            om_th = np.full(4, w[np.argmin(res)].real)
            what = f'eigenmode (residual {res.min():.1e})'
        else:
            om_th = None
            what = f'max growth rate {np.max(w.imag):.6f}'
        msg = (f'({mx},{my},{mz}) {what}: |f - expm(iRt) x0|/|x0(t)| = {err:.1e}, '
               f'omg = -K phi err {omg_err:.1e}, leakage {leak:.1e}, spec roots diff {spec:.1e}, t = {t}')
        passed = err < 1e-6 and omg_err < 1e-12 and leak < 1e-12 and spec < 1e-9
        if om_th is not None:
            ferr = np.max(np.abs(om - om_th))/np.max(np.abs(om_th))
            msg += f'; omega {np.array2string(om, precision=9)} vs {np.array2string(om_th, precision=9)}'
            passed &= ferr < 1e-6
        report(f'V7 {c}', passed, msg)

        s = series(c)
        e = sum(s[k] for k in ENERGY)
        dot = sum(s[k] for k in DOT)[1:]
        src = sum(s[k] for k in SOURCES)[1:]
        pmax = max(np.max(np.abs(s[k])) for k in SOURCES)
        if pmax < 1e-10*e[0]:
            # ky = 0: every source is proportional to ky. Real eigenvector: all
            # components share one phase, so every Im(x y^*) source vanishes
            # and E stays constant up to the RK3 amplitude error
            drift = np.max(np.abs(e - e[0]))/e[0]
            report(f'energy {c}', drift < 1e-8,
                   f'max|p_*|/E0 = {pmax/e[0]:.1e}, max|E - E0|/E0 = {drift:.1e}')
        else:
            # generic initial vector: E changes by the sources
            err = np.max(np.abs(dot - src))/max(pmax, np.max(np.abs(dot)))
            report(f'budget {c}', err < 1e-6,
                   f'max|dE/dt - sum p_*|/max(|p_*|, |dE/dt|) = {err:.1e}, '
                   f'max|p_*| = {pmax:.3e}, max|dE/dt| = {np.max(np.abs(dot)):.3e}')


if SET == 'linear':
    linear_set()
else:
    main_set()
print('ALL PASS' if ok else 'SOME CHECKS FAILED')
sys.exit(0 if ok else 1)

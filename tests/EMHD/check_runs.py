#!/usr/bin/env python3
"""
Evaluate the EMHD validation runs created by tests/EMHD/setup_runs.sh.

  python3 tests/EMHD/check_runs.py <runs dir>          # main set
  python3 tests/EMHD/check_runs.py <runs dir> linear   # V7 (helical whistler modes)

Needs numpy and netCDF4. Reads <case>/gallope.out.nc and <case>/restart/*.dat
(raw complex128, Fortran order (nkz, nky, nkx)).

Model, as implemented in src/model/EMHD (per Fourier mode, d = 1 + de^2 k^2):
  d db/dt = -i k x (b* x u) + f - eta (k^2/k2_max)^eta_exp b,
  b* = d b,  u = -i k x b,
  W  = 1/2 sum d |b|^2                      (wmag_sum, invariant for eta = f = 0),
  H* = sum d^2 Re(A . conj b), A = i k x b/k^2  (hel_sum, invariant if also b0 = 0).
"""
import os
import sys

import numpy as np
from netCDF4 import Dataset

RUNS = sys.argv[1]
SET = sys.argv[2] if len(sys.argv) > 2 else 'main'
NLX, NLY, NLZ = 64, 64, 32
NKX, NKY, NKZ = NLX, NLY, NLZ//2 + 1
FIELDS = ['bx', 'by', 'bz']
FILES = [f + '.dat' for f in FIELDS]
ETA, ETA_EXP, DE, B0 = 100.0, 2.0, 0.1, (0.0, 0.0, 1.0)
ok = True


def report(name, passed, msg):
    global ok
    ok &= bool(passed)
    print(f"[{'PASS' if passed else 'FAIL'}] {name}: {msg}")


def missing(case, name):
    """Report FAIL and return True when <case> has no NetCDF output yet."""
    if os.path.exists(os.path.join(RUNS, case, 'gallope.out.nc')):
        return False
    report(f'{name} {case}', False, 'no gallope.out.nc')
    return True


def params(case):
    """Scalar parameters recorded in <case>/gallope.out.nc."""
    with Dataset(os.path.join(RUNS, case, 'gallope.out.nc')) as nc:
        return {v: float(nc.variables[v][...]) for v in nc.variables if nc.variables[v].size == 1
                and nc.variables[v].dimensions != ('tt',)}


def var(case, name):
    with Dataset(os.path.join(RUNS, case, 'gallope.out.nc')) as nc:
        return np.array(nc.variables[name][:])


def check_params(case, expected, b0):
    """V0: the namelists were read in full (nvfortran stops a group at ';' silently)."""
    if missing(case, 'V0'):
        return
    got = params(case)
    bad = {k: (got.get(k), v) for k, v in expected.items() if got.get(k) != v}
    b0_got = tuple(var(case, 'b0')[0])
    if b0_got != tuple(b0):
        bad['b0(t=0)'] = (b0_got, tuple(b0))
    report(f'V0 params {case}', not bad, 'as intended' if not bad else f'(got, want) {bad}')


def series(case):
    with Dataset(os.path.join(RUNS, case, 'gallope.out.nc')) as nc:
        out = {v: np.array(nc.variables[v][:]) for v in nc.variables
               if nc.variables[v].dimensions == ('tt',)}
    return out


def rdir(case):
    """Restart directory of <case>; e_rs2 continues e_rs1 in ../e_rs1/restart/."""
    return os.path.join(RUNS, 'e_rs1' if case == 'e_rs2' else case, 'restart')


def field(case, name):
    a = np.fromfile(os.path.join(rdir(case), name + '.dat'), dtype=np.complex128)
    return a.reshape(NKX, NKY, NKZ)   # [i, j, k]


def restart_time(case):
    with open(os.path.join(rdir(case), 'time.dat')) as f:
        return float(f.readlines()[1].split()[0])


def same_bytes(a, b, names):
    return [n for n in names
            if open(os.path.join(rdir(a), n), 'rb').read()
            != open(os.path.join(rdir(b), n), 'rb').read()]


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


def drift(case, key):
    q = series(case)[key]
    return abs(q[-1] - q[0])/abs(q[0])


def check_divb(cases):
    worst, where = 0.0, None
    for c in cases:
        if missing(c, 'div b'):
            continue
        d = np.max(series(c)['divb_rel'])
        if d >= worst:
            worst, where = d, c
    report('div b', worst < 1e-12, f'max divb_rel = {worst:.1e} ({where}) over {len(cases)} cases')


def main_set():
    base = {'eta': ETA, 'eta_exp': ETA_EXP, 'de': DE}
    check_params('e_ot3_np1', base, B0)
    check_params('e_frc_np1', base, B0)
    check_params('e_ideal_dt2', {**base, 'eta': 0.0}, B0)
    check_params('e_rnd_ideal_dt2', {**base, 'eta': 0.0}, (0.0, 0.0, 0.0))

    # V2 smoke: series finite, final state finite
    cases = [c for c in sorted(os.listdir(RUNS))
             if c.startswith('e_') and not c.startswith('e_lin') and c != 'e_pm2']
    for c in cases:
        if missing(c, 'V2'):
            continue
        s = series(c)
        finite = all(np.all(np.isfinite(v)) for v in s.values()) and \
            all(np.all(np.isfinite(field(c, f))) for f in FIELDS)
        report(f'V2 {c}', finite, f"{len(s['tt'])} samples, t_end = {s['tt'][-1]:.4f}")

    # V3 decomposition: 1 GPU vs 4 GPUs, same OT3
    keys = ['wmag_sum', 'wmag_dot_sum', 'wmag_dissip_sum', 'hel_sum']
    d, n = max_rel_series('e_ot3_np1', 'e_ot3_np4', keys)
    report('V3 series np1 vs np4', d < 1e-12, f'max rel diff {d:.2e} over {n} samples')
    d = max_rel_fields('e_ot3_np1', 'e_ot3_np4', FIELDS)
    report('V3 restart np1 vs np4', d < 1e-12, f'max rel diff {d:.2e}')

    # V4: EMHD runs with P_m = P_s = 1 only, so the redundancy test is N/A; check that
    # P_m = 2 is rejected before any step. With nm = 1 the common grid check fires
    # first; EMHD_specific has its own guard behind it.
    out = open(os.path.join(RUNS, 'e_pm2', 'out.std')).read()
    msgs = [m for m in ('nm has to divide evenly by the comm_m rank count (P_m)',
                        'EMHD supports P_m = P_s = 1 only') if m in out]
    report('V4 P_m = 2 rejected', msgs and 'Program completed' not in out,
           f'aborted with "{msgs[0]}"' if msgs else 'no abort message')

    # V5 restart: 100 + 100 steps (final state in e_rs1/restart) vs 200 steps
    diff = same_bytes('e_rs1', 'e_ot3_np1', FILES)
    t1, t2 = restart_time('e_rs1'), restart_time('e_ot3_np1')
    report('V5 restart 100+100 vs 200 (bitwise)', not diff and t1 == t2,
           f'tt {t1} vs {t2}; ' + ('identical' if not diff else f'differ: {diff}'))

    # V6 ideal invariants at dt and dt/2: third-order drift, ratio ~ 8
    for (c1, c2), key, label in [(('e_ideal_dt1', 'e_ideal_dt2'), 'wmag_sum', 'W (OT3, b0 = z)'),
                                 (('e_rnd_ideal_dt1', 'e_rnd_ideal_dt2'), 'wmag_sum', 'W (random, b0 = 0)'),
                                 (('e_rnd_ideal_dt1', 'e_rnd_ideal_dt2'), 'hel_sum', 'H* (random, b0 = 0)')]:
        if missing(c1, 'V6') or missing(c2, 'V6'):
            continue
        d1, d2 = drift(c1, key), drift(c2, key)
        r = d1/d2 if d2 > 0 else np.inf
        report(f'V6 {label} drift ratio', 6 < r < 10 or max(d1, d2) < 1e-13,
               f'drift dt={d1:.3e}, dt/2={d2:.3e}, ratio {r:.2f} (RK3: 8); '
               f'{key}(0) = {series(c2)[key][0]:.6e}')

    # Energy budget dW/dt = p_ext - dissip. For the decaying run, the dissipation-free
    # twin's dW/dt (the eSSPIFRK3 drift) is subtracted and the residual is measured
    # against the dissipation; the forced run is measured against max|dW/dt|.
    for c, ideal in [('e_ot3_np1', 'e_ideal_dt2'), ('e_frc_np1', None)]:
        if missing(c, 'budget') or (ideal and missing(ideal, 'budget')):
            continue
        s = series(c)
        dot, dis, inj = s['wmag_dot_sum'][1:], s['wmag_dissip_sum'][1:], s['p_ext_sum'][1:]
        res = dot + dis - inj
        if ideal:
            res0 = np.max(np.abs(res))/np.max(np.abs(dis))
            res -= series(ideal)['wmag_dot_sum'][1:]
            err = np.max(np.abs(res))/np.max(np.abs(dis))
            report(f'budget {c}', err < 1e-3,
                   f'max |dW/dt + dissip - dW/dt({ideal})|/max dissip = {err:.2e} '
                   f'(without the twin {res0:.2e}; max dissip {np.max(dis):.2e})')
        else:
            err = np.max(np.abs(res))/np.max(np.abs(dot))
            report(f'budget {c}', err < 1e-2, f'max |dW/dt - (inj - dissip)|/max|dW/dt| = {err:.2e}')

    # V9 forcing: injected power vs ene_inj = 1, decomposition independence
    if not (missing('e_frc_np1', 'V9') or missing('e_frc_np4', 'V9')):
        inj = series('e_frc_np1')['p_ext_sum'][1:]
        report('V9 injection', np.all(np.abs(inj - 1.0) < 0.1),
               f'p_ext in [{inj.min():.4f}, {inj.max():.4f}] (ene_inj = 1)')
        d, n = max_rel_series('e_frc_np1', 'e_frc_np4', keys + ['p_ext_sum'])
        report('V9 series np1 vs np4', d < 1e-10, f'max rel diff {d:.2e} over {n} samples')

    check_divb(cases)


def helical(mx, my, mz, s, amp):
    """e_s of fields.F90 init_single_mode (box 2 pi, so k = mode numbers)."""
    sgn = -1 if mz < 0 else 1
    k = sgn*np.array([mx, my, mz], dtype=float)
    kh = k/np.linalg.norm(k)
    a = np.array([0.0, 0.0, 1.0]) if abs(kh[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    e1 = a - a.dot(kh)*kh
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(kh, e1)
    es = amp*(e1 - 1j*s*e2)/np.sqrt(2.0)
    # independent check that e_s is the helical eigenvector: k x e_s = i s |k| e_s
    assert np.allclose(np.cross(k, es), 1j*s*np.linalg.norm(k)*es, atol=1e-15)
    return k, es


def linear_set():
    amp = 1e-2
    runs = [('e_lin_001_de0', (0, 0, 1), +1, 0.0), ('e_lin_001_de01', (0, 0, 1), -1, 0.1),
            ('e_lin_3m42_de01', (3, -4, 2), +1, 0.1), ('e_lin_3m42_de1', (3, -4, 2), -1, 1.0)]
    for c, _, _, de in runs:
        check_params(c, {'eta': ETA, 'eta_exp': ETA_EXP, 'de': de}, B0)
    for c, (mx, my, mz), s, de in runs:
        if missing(c, 'V7'):
            continue
        k, es = helical(mx, my, mz, s, amp)
        kx, ky, kz = var(c, 'kx'), var(c, 'ky'), var(c, 'kz')
        k2_max = np.max(np.abs(kx))**2 + np.max(np.abs(ky))**2 + np.max(np.abs(kz))**2
        k2 = k.dot(k)
        d = 1.0 + de*de*k2
        om = s*k.dot(B0)*np.sqrt(k2)/d
        gam = ETA*(k2/k2_max)**ETA_EXP/d
        t = restart_time(c)
        idx = (mx % NKX, my % NKY, mz)
        f = np.array([field(c, n)[idx] for n in FIELDS])
        # f = es exp(i om t - gam t); unwrap the phase around the expected one
        j = np.argmax(np.abs(es))
        ph = np.angle(f[j]/es[j]*np.exp(-1j*om*t))
        om_num = om + ph/t
        err_om = abs(om_num - om)/abs(om)
        err_amp = abs(np.abs(f[j])/np.abs(es[j]) - np.exp(-gam*t))/np.exp(-gam*t)
        err_vec = np.max(np.abs(f - es*np.exp((1j*om - gam)*t)))/amp
        leak, b0_err = 0.0, 0.0
        for n, b0 in zip(FIELDS, B0):
            g = field(c, n).copy()
            b0_err = max(b0_err, abs(g[0, 0, 0] - b0))
            g[idx] = 0.0
            g[0, 0, 0] = 0.0
            leak = max(leak, np.max(np.abs(g))/amp)
        divb = np.max(series(c)['divb_rel'])
        report(f'V7 {c}', err_om < 1e-6 and err_amp < 1e-6 and err_vec < 1e-6 and leak < 1e-10
               and b0_err < 1e-14 and divb < 1e-12,
               f'omega {om_num:.12f} vs {om:.12f} (rel err {err_om:.1e}), '
               f'decay {gam:.4e}: |amp| rel err {err_amp:.1e}, vector err {err_vec:.1e}, '
               f'leakage {leak:.1e}, |b0 - b0(0)| {b0_err:.1e}, divb_rel {divb:.1e}, t = {t}')


if SET == 'linear':
    linear_set()
else:
    main_set()
print('ALL PASS' if ok else 'SOME CHECKS FAILED')
sys.exit(0 if ok else 1)

#!/usr/bin/env python3
"""
Evaluate the HRMHD validation runs created by tests/HRMHD/setup_runs.sh.

  python3 tests/HRMHD/check_runs.py <runs dir>          # main set
  python3 tests/HRMHD/check_runs.py <runs dir> linear   # V7 and rho -> 0

Needs numpy and netCDF4. Reads <case>/gallope.out.nc and <case>/restart/*.dat
(raw complex128, Fortran order (nkz, nky, nkx)).
"""
import os
import sys

import numpy as np
from netCDF4 import Dataset

RUNS = sys.argv[1]
SET = sys.argv[2] if len(sys.argv) > 2 else 'main'
NLX, NLY, NLZ = 64, 64, 32
NKX, NKY, NKZ = NLX, NLY, NLZ//2 + 1
FIELDS = ['phi', 'omg', 'psi', 'upa', 'bpa']
ENERGY = ['upe2_sum', 'bpe2_sum', 'upa2_sum', 'bpa2_sum']
DOT = [e.replace('2_sum', '2dot_sum') for e in ENERGY]
DISSIP = [e.replace('2_sum', '2dissip_sum') for e in ENERGY]
ok = True


def report(name, passed, msg):
    global ok
    ok &= passed
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


DISS = [f'{a}{b}_{c}' for a in ('nu', 'eta') for b in ('pe', 'pa') for c in ('x', 'z')]


def check_params(case, expected):
    """V0: the namelist was read in full (nvfortran stops a group at ';' silently)."""
    if missing(case, 'V0'):
        return
    got = params(case)
    bad = {k: (got.get(k), v) for k, v in expected.items() if got.get(k) != v}
    report(f'V0 params {case}', not bad, 'as intended' if not bad else f'(got, want) {bad}')


def series(case):
    with Dataset(os.path.join(RUNS, case, 'gallope.out.nc')) as nc:
        out = {v: np.array(nc.variables[v][:]) for v in nc.variables
               if nc.variables[v].dimensions == ('tt',)}
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


def main_set():
    # V0: the input reached the code (base and ideal cases)
    base = {'rho': 0.5, 'sgm': 1.0, **{d: (10.0 if d.endswith('x') else 1.0) for d in DISS},
            **{d + '_exp': 4.0 for d in DISS}}
    check_params('h_ot3_np1', base)
    check_params('h_frc_els_np1', base)
    check_params('h_ideal_dt2', {**base, **{d: 0.0 for d in DISS}})

    # V2 smoke: series finite, final state finite
    for c in sorted(os.listdir(RUNS)):
        if not c.startswith('h_') or c.startswith('h_lin'):
            continue
        if missing(c, 'V2'):
            continue
        s = series(c)
        finite = all(np.all(np.isfinite(v)) for v in s.values())
        report(f'V2 {c}', finite, f"{len(s['tt'])} samples, t_end = {s['tt'][-1]:.4f}")

    # V3 decomposition: 1 GPU vs 4 GPUs, same OT3
    keys = ENERGY + DISSIP + ['xhl_sum', 'KAW_dissip_prp_sum', 'ICW_dissip_prp_sum']
    d, n = max_rel_series('h_ot3_np1', 'h_ot3_np4', keys)
    report('V3 series np1 vs np4', d < 1e-12, f'max rel diff {d:.2e} over {n} samples')
    d = max_rel_fields('h_ot3_np1', 'h_ot3_np4', FIELDS)
    report('V3 restart np1 vs np4', d < 1e-12, f'max rel diff {d:.2e}')

    # V4 redundancy: P_m = 2 must reproduce P_m = 1 bit for bit
    files = [f + '.dat' for f in FIELDS]
    diff = same_bytes('h_ot3_np2', 'h_ot3_pm2', files)
    report('V4 restart np2 vs pm2 (bitwise)', not diff, 'identical' if not diff else f'differ: {diff}')
    pm2 = open(os.path.join(RUNS, 'h_ot3_pm2', 'out.std')).read()
    report('V4 check_redundant_consistency', 'FAIL' not in pm2 and 'PASS' in pm2,
           f"{pm2.count('PASS')} PASS lines")

    # V5 restart: 100 + 100 steps (final state in h_rs1/restart) vs 200 steps
    diff = same_bytes('h_rs1', 'h_ot3_np1', files)
    t1, t2 = restart_time('h_rs1'), restart_time('h_ot3_np1')
    report('V5 restart 100+100 vs 200 (bitwise)', not diff and t1 == t2,
           f'tt {t1} vs {t2}; ' + ('identical' if not diff else f'differ: {diff}'))

    # V6 ideal invariants: drift at dt and dt/2
    for inv in ['W', 'X']:
        drift = []
        for c in ['h_ideal_dt1', 'h_ideal_dt2']:
            s = series(c)
            q = sum(s[e] for e in ENERGY) if inv == 'W' else s['xhl_sum']
            drift.append(abs(q[-1] - q[0])/abs(q[0]))
        r = drift[0]/drift[1] if drift[1] > 0 else np.inf
        report(f'V6 {inv} drift ratio', 6 < r < 10 or max(drift) < 1e-13,
               f'drift dt={drift[0]:.3e}, dt/2={drift[1]:.3e}, ratio {r:.2f} (RK3: 8)')

    # Energy budget: d/dt sum E = -sum dissipation + injection.
    # The hyper-dissipation of the decaying runs (~1e-9) is as small as the eSSPIFRK3
    # truncation drift at dt = 1e-3 (~3.5e-9 per unit time), so for those the drift of
    # the same run without dissipation is subtracted and the residual is measured
    # against the dissipation itself. Forced runs are measured against dW/dt ~ ene_inj.
    for c, ideal in [('h_ot3_np1', 'h_ideal_dt2'), ('h_rho0', 'h_rho0_ideal'),
                     ('h_frc_np1', None), ('h_frc_els_np1', None)]:
        if missing(c, 'budget') or (ideal and missing(ideal, 'budget')):
            continue
        s = series(c)
        dot = sum(s[k] for k in DOT)[1:]
        dis = sum(s[k] for k in DISSIP)[1:]
        res = dot + dis - (s['p_aw_sum'] + s['p_compr_sum'])[1:]
        if ideal:
            res -= sum(series(ideal)[k] for k in DOT)[1:]
            err = np.max(np.abs(res))/np.max(np.abs(dis))
            report(f'budget {c}', err < 1e-3,
                   f'max |dW/dt + dissip - dW/dt({ideal})|/max dissip = {err:.2e}')
        else:
            err = np.max(np.abs(res))/np.max(np.abs(dot))
            report(f'budget {c}', err < 1e-2, f'max |dW/dt - (inj - dissip)|/max|dW/dt| = {err:.2e}')

    # V9 forcing: injected power vs ene_inj = 1, decomposition independence
    s = series('h_frc_np1')
    inj = (s['p_aw_sum'] + s['p_compr_sum'])[1:]
    report('V9 injection', np.all(np.abs(inj - 1.0) < 0.1),
           f'p_aw + p_compr in [{inj.min():.4f}, {inj.max():.4f}] (ene_inj = 1)')
    d, n = max_rel_series('h_frc_np1', 'h_frc_np4', ENERGY + ['p_aw_sum', 'p_compr_sum'])
    report('V9 series np1 vs np4', d < 1e-10, f'max rel diff {d:.2e} over {n} samples')

    # Elsasser forcing: P+ = ene_inj(1+xhl_inj)/2, P- = ene_inj(1-xhl_inj)/2, so the
    # energy injection is ene_inj = 1 and the X injection P+ - P- is ene_inj*xhl_inj = 0.5.
    # upa, bpa are zeroed at the stirred modes, so the -(rho/sgm) kprp^2 upa phi term
    # of X gets no injection there and p_xhl is the whole X injection.
    if missing('h_frc_els_np1', 'V9 Elsasser') or missing('h_frc_els_np4', 'V9 Elsasser'):
        return
    s = series('h_frc_els_np1')
    inj, xinj, pc = s['p_aw_sum'][1:], s['p_xhl_sum'][1:], s['p_compr_sum'][1:]
    report('V9 Elsasser injection', np.all(np.abs(inj - 1.0) < 0.1) and np.all(np.abs(xinj - 0.5) < 0.1)
           and np.all(pc == 0.0),
           f'p_aw in [{inj.min():.4f}, {inj.max():.4f}], p_xhl in [{xinj.min():.4f}, {xinj.max():.4f}] '
           f'(1, 0.5), max |p_compr| {np.abs(pc).max():.1e}')
    d, n = max_rel_series('h_frc_els_np1', 'h_frc_els_np4', ENERGY + ['p_aw_sum', 'p_xhl_sum'])
    report('V9 Elsasser series np1 vs np4', d < 1e-10, f'max rel diff {d:.2e} over {n} samples')
    # zeroing_at_forcing_wavenumber (device kernel) keeps the unforced upa, bpa at zero on the
    # stirred modes |kx|, |ky| <= 1, kz = 1, 2 (kmin/kmax of the input); with 4 GPUs ky = -1
    # lives on the last rank.
    for c in ['h_frc_els_np1', 'h_frc_els_np4']:
        z = max(abs(field(c, n)[kx % NKX, ky % NKY, kz]) for n in ['upa', 'bpa']
                for kx in (-1, 0, 1) for ky in (-1, 0, 1) for kz in (1, 2))
        report(f'V9 zeroing {c}', z == 0.0, f'max |upa|, |bpa| on the 18 stirred modes = {z:.1e}')


def linear_set():
    rho, sgm = 0.5, 1.0
    check_params('h_lin_kaw_401', {'rho': rho, 'sgm': sgm, **{d: 0.0 for d in DISS}})
    check_params('r_ot3', {**{d: (10.0 if d.endswith('x') else 1.0) for d in DISS if 'pe' in d},
                           **{d + '_exp': 4.0 for d in DISS if 'pe' in d}})
    for c, (mx, my, mz) in [('h_lin_kaw_401', (4, 0, 1)), ('h_lin_icw_401', (4, 0, 1)),
                            ('h_lin_kaw_3m42', (3, -4, 2)), ('h_lin_icw_3m42', (3, -4, 2))]:
        kp2 = mx*mx + my*my
        S = 1 + sgm**2 + rho**2*kp2
        D = np.sqrt(S*S - 4*sgm**2)
        lam = np.sqrt((S + D)/2 if 'kaw' in c else (S - D)/2)
        amp = np.array([1.0, lam, sgm*(lam**2 - 1)/(rho*lam), (lam**2 - 1)/rho])
        t = restart_time(c)
        idx = (mx % NKX, my % NKY, mz)
        f = np.array([field(c, n)[idx] for n in ['phi', 'psi', 'upa', 'bpa']])
        # omega from the phase of each component: f = amp*exp(i kz lam t)
        om = np.angle(f/amp)/t
        om_th = mz*lam
        err = np.max(np.abs(om - om_th))/om_th
        mag = np.max(np.abs(np.abs(f) - np.abs(amp))/np.abs(amp))
        leak = 0.0
        for n in ['phi', 'psi', 'upa', 'bpa']:
            g = field(c, n).copy()
            g[idx] = 0.0
            leak = max(leak, np.max(np.abs(g)))
        report(f'V7 {c}', err < 1e-6 and mag < 1e-6 and leak < 1e-12,
               f'omega {om.mean():.12f} vs {om_th:.12f} (rel err {err:.1e}), '
               f'|amp| rel err {mag:.1e}, leakage {leak:.1e}, t = {t}')

    # rho -> 0: HRMHD(rho = 0) phi, psi equal RMHD; upa, bpa stay zero
    d = max_rel_fields('h_rho0', 'r_ot3', ['phi', 'omg', 'psi'])
    zero = max(np.max(np.abs(field('h_rho0', n))) for n in ['upa', 'bpa'])
    report('rho->0 HRMHD vs RMHD', d < 1e-12 and zero == 0.0,
           f'phi/omg/psi max rel diff {d:.2e}, max |upa|,|bpa| = {zero:.1e}')


if SET == 'linear':
    linear_set()
else:
    main_set()
print('ALL PASS' if ok else 'SOME CHECKS FAILED')
sys.exit(0 if ok else 1)

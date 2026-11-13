#!/usr/bin/env python3
"""
Evaluate the v1.1 release checks created by tests/release/setup_runs.sh.

  python3 tests/release/check_runs.py <runs/rel dir> [r12] [r3] [r4] [r5]

With no set named, every set is evaluated. Needs numpy and netCDF4.
Reads <case>/gallope.in, <case>/out.std, <case>/gallope.out.nc,
<case>/restart/*.dat and <case>/out3d/*.dat (raw complex128, Fortran order
(nkz, nky, nkx); a 3D output file is a sequence of such records).
"""
import math
import os
import re
import sys

import numpy as np
from netCDF4 import Dataset

RUNS = sys.argv[1]
SETS = sys.argv[2:] or ['r12', 'r3', 'r4', 'r5']
TOL = 1e-12
ok = True


def report(name, passed, msg):
    """passed: True, False or None (information only)."""
    global ok
    if passed is not None:
        ok &= passed
    tag = 'INFO' if passed is None else ('PASS' if passed else 'FAIL')
    print(f'[{tag}] {name}: {msg}')


def path(case, *p):
    return os.path.join(RUNS, case, *p)


def missing(case, name):
    """Report FAIL and return True when <case> has no NetCDF output."""
    if os.path.exists(path(case, 'gallope.out.nc')):
        return False
    report(f'{name} {case}', False, 'no gallope.out.nc')
    return True


def out_std(case):
    p = path(case, 'out.std')
    return open(p, errors='replace').read() if os.path.exists(p) else ''


# ---------------------------------------------------------------- input --
def namelist_numbers(case):
    """Scalar numeric items of <case>/gallope.in (arrays, strings, logicals skipped)."""
    out = {}
    for line in open(path(case, 'gallope.in')):
        line = line.split('!')[0]
        if '=' not in line:
            continue
        k, v = (s.strip() for s in line.split('=', 1))
        if not re.fullmatch(r'[A-Za-z_]\w*', k):
            continue
        try:
            out[k.lower()] = float(v.rstrip(',').lower().replace('d', 'e'))
        except ValueError:
            pass
    return out


def kdims(case):
    """(nkx, nky, nkz) of <case>: the C-order shape of one complex field."""
    p = namelist_numbers(case)
    return int(p['nlx']), int(p['nly']), int(p['nlz'])//2 + 1


def params(case):
    """Scalar parameters recorded in <case>/gallope.out.nc."""
    with Dataset(path(case, 'gallope.out.nc')) as nc:
        return {v.lower(): (float(nc.variables[v][...]), nc.variables[v].dtype)
                for v in nc.variables
                if nc.variables[v].size == 1 and nc.variables[v].dimensions != ('tt',)}


def check_input(case, skip=()):
    """V0: every numeric input item recorded in the nc has the input value, and
    out.std has no namelist read failure other than the routine force one."""
    want, got = namelist_numbers(case), params(case)
    bad, n = {}, 0
    for k, v in want.items():
        if k in got and k not in skip:
            n += 1
            g, dtype = got[k]
            rtol = 1e-12 if dtype == np.float64 else 1e-6
            if not math.isclose(g, v, rel_tol=rtol, abs_tol=0.0):
                bad[k] = (g, v)
    std = out_std(case)
    failed = sorted({l.strip() for l in std.splitlines() if 'failed' in l})
    odd = [l for l in failed if l != 'Reading force_parameters failed']
    decreased = std.count('dt is decreased')
    passed = not bad and not odd
    msg = f'{n} items match' if not bad else f'(got, want) {bad}'
    if failed:
        msg += f'; out.std: {failed}'
    if decreased:
        msg += f'; dt decreased {decreased} times'
    report(f'V0 {case}', passed, msg)


# --------------------------------------------------------------- output --
def series(case):
    with Dataset(path(case, 'gallope.out.nc')) as nc:
        return {v: np.array(nc.variables[v][:], dtype=float) for v in nc.variables
                if nc.variables[v].dimensions == ('tt',)}


def restart_files(case):
    d = path(case, 'restart')
    return sorted(f for f in os.listdir(d) if f.endswith('.dat') and f != 'time.dat')


def restart_time(case):
    with open(path(case, 'restart', 'time.dat')) as f:
        return float(f.readlines()[1].split()[0])


def raw(case, *p):
    return open(path(case, *p), 'rb').read()


def differing(a, b, sub, names):
    return [n for n in names if raw(a, sub, n) != raw(b, sub, n)]


def max_rel_arrays(xa, xb):
    scale = np.max(np.abs(xa)) or 1.0
    return float(np.max(np.abs(xa - xb))/scale)


def max_rel_fields(a, b, names=None, sub='restart'):
    names = names or restart_files(a)
    worst = 0.0
    for n in names:
        fa = np.fromfile(path(a, sub, n), dtype=np.complex128)
        fb = np.fromfile(path(b, sub, n), dtype=np.complex128)
        if fa.size != fb.size:
            return float('inf')
        worst = max(worst, max_rel_arrays(fa, fb))
    return worst


def compare_series(a, b):
    """Max rel. diff over the common series and the first tt where they differ
    at all and by more than TOL."""
    sa, sb = series(a), series(b)
    keys = sorted(set(sa) & set(sb))
    n = min(len(sa['tt']), len(sb['tt']))
    per = np.zeros(n)
    worst_key = None
    for k in keys:
        scale = np.max(np.abs(sa[k][:n])) or 1.0
        d = np.abs(sa[k][:n] - sb[k][:n])/scale
        if np.max(d) > np.max(per):
            worst_key = k
        per = np.maximum(per, d)
    tt = sa['tt'][:n]
    first_any = tt[np.argmax(per > 0)] if np.any(per > 0) else None
    first_tol = tt[np.argmax(per > TOL)] if np.any(per > TOL) else None
    return float(np.max(per)), worst_key, n, len(keys), first_any, first_tol


def fmt_t(t):
    return 'never' if t is None else f'{t:.3f}'


def compare(name, a, b, bitwise, tol=TOL, info=False):
    """Restart (bitwise or rel. diff) and series of case a vs case b.
    info: report the numbers without a verdict."""
    if missing(a, name) or missing(b, name):
        return None
    names = restart_files(a)
    diff = differing(a, b, 'restart', names)
    ta, tb = restart_time(a), restart_time(b)
    d, key, n, nk, t_any, t_tol = compare_series(a, b)
    dr = max_rel_fields(a, b, names)
    msg = (f"restart {'identical' if not diff else f'differ ({len(diff)}/{len(names)}), max rel {dr:.2e}'}"
           f"; tt {ta!r} vs {tb!r}; series max rel {d:.2e} ({key}) over {n} samples x {nk} series"
           f"; first diff at t = {fmt_t(t_any)}, > {tol:.0e} at t = {fmt_t(t_tol)}")
    if info:
        passed = None
    elif bitwise:
        passed = not diff and ta == tb
    else:
        passed = d < tol and dr < tol and ta == tb
    report(f'{name} {a} vs {b}', passed, msg)
    return not diff


def smoke(case):
    if missing(case, 'smoke'):
        return
    s = series(case)
    finite = all(np.all(np.isfinite(v)) for v in s.values())
    fin = 'Finished the main loop' in out_std(case)
    report(f'smoke {case}', finite and fin,
           f"{len(s['tt'])} samples, t = {s['tt'][0]:.3f} .. {s['tt'][-1]:.4f}, "
           f"{'finished' if fin else 'NOT finished'}")


def redundancy(case):
    std = out_std(case)
    p, f = std.count('PASS'), std.count('FAIL')
    report(f'check_redundant_consistency {case}', p > 0 and f == 0, f'{p} PASS, {f} FAIL')


def cases(prefix):
    return sorted(c for c in os.listdir(RUNS) if c.startswith(prefix))


# ------------------------------------------------------------------ R1/R2 --
def set_r12():
    for c in cases('r1_') + cases('r2_'):
        if not missing(c, 'V0'):
            check_input(c)
            smoke(c)
    # R1: undriven, v1.0 vs v1.1 bitwise; np1 vs np4 in v1.1
    for k in ('mi', 'ms', 'rm'):
        for np_ in (1, 4):
            compare('R1 v1.0 vs v1.1', f'r1_{k}_v10_np{np_}', f'r1_{k}_v11_np{np_}', bitwise=True)
        compare('R1 np1 vs np4 (v1.1)', f'r1_{k}_v11_np1', f'r1_{k}_v11_np4', bitwise=False)
    # R2: driven. v1.1 np1 vs np4 < 1e-12; v1.0 differs (forcing bug)
    for k in ('mi', 'rm', 'rk0'):
        compare('R2 np1 vs np4 (v1.1)', f'r2_{k}_v11_np1', f'r2_{k}_v11_np4', bitwise=False)
    for k in ('mi', 'rm'):
        compare('R2 np1 vs np4 (v1.0)', f'r2_{k}_v10_np1', f'r2_{k}_v10_np4', bitwise=False, info=True)
    # R2: does a driven 1-GPU v1.1 run reproduce v1.0?
    for k in ('mi', 'rm', 'rk0'):
        compare('R2 v1.0 vs v1.1 on 1 GPU', f'r2_{k}_v10_np1', f'r2_{k}_v11_np1',
                bitwise=True, info=(k == 'rk0'))


# --------------------------------------------------------------------- R3 --
def set_r3():
    for c in cases('r3'):
        if not missing(c, 'V0'):
            check_input(c)
            smoke(c)
    compare('R3 np1 vs np4', 'r3_kf_np1', 'r3_kf_np4', bitwise=False)
    compare('R3 P_m = 1 vs 2 (P_fft = 2)', 'r3_kf_np2', 'r3_kf_pm2', bitwise=True)
    compare('R3 P_m = 1 vs 4 (P_fft = 1)', 'r3_kf_np1', 'r3_kf_pm4', bitwise=True)
    compare('R3 P_m = 1 vs 2 (P_fft = 2, random)', 'r3_kr_np2', 'r3_kr_pm2', bitwise=True)
    for c in ('r3_kf_pm2', 'r3_kf_pm4', 'r3_kr_pm2'):
        redundancy(c)
    # rs2 restarts from r3_kr_rs1/restart and writes its final state there
    if not (missing('r3_kr_rs1', 'R3 restart') or missing('r3_kr_full', 'R3 restart')):
        names = restart_files('r3_kr_full')
        diff = differing('r3_kr_rs1', 'r3_kr_full', 'restart', names)
        t1, t2 = restart_time('r3_kr_rs1'), restart_time('r3_kr_full')
        report('R3 restart 100+100 vs 200 (bitwise)', not diff and t1 == t2,
               f"{'identical' if not diff else f'differ: {diff}'}; tt {t1!r} vs {t2!r}")
    # R3': KRMHD with g undriven reproduces RMHD in phi and psi
    for k in ('ot2', 'els'):
        a, b = f'r3p_k_{k}', f'r3p_r_{k}'
        if missing(a, "R3'") or missing(b, "R3'"):
            continue
        names = ['phi.dat', 'omg.dat', 'psi.dat']
        diff = differing(a, b, 'restart', names)
        d = max_rel_fields(a, b, names)
        g = np.fromfile(path(a, 'restart', 'g.dat'), dtype=np.complex128)
        nm = g.size//math.prod(kdims(a))
        g0 = np.max(np.abs(g.reshape(nm, -1)[0]))
        g1 = np.max(np.abs(g.reshape(nm, -1)[1:]))
        dser, key, n, nk, _, _ = compare_series(a, b)
        report(f"R3' KRMHD vs RMHD {k}", d < TOL,
               f"phi, omg, psi {'identical' if not diff else f'max rel {d:.2e}'}; "
               f"common series max rel {dser:.2e} ({key}, {nk} series); "
               f"max|g| m = 0: {g0:.2e}, m >= 1: {g1:.2e}")


# --------------------------------------------------------------------- R4 --
def out3d_names(case):
    d = path(case, 'out3d')
    return sorted(f for f in os.listdir(d) if f.endswith('.dat') and f != 'time.dat')


def out3d_times(case):
    t = []
    for l in open(path(case, 'out3d', 'time.dat')):
        try:
            t.append(float(l.split()[0]))
        except (ValueError, IndexError):
            pass
    return t


def records(case, name):
    a = np.fromfile(path(case, 'out3d', name), dtype=np.complex128)
    shape = kdims(case)
    rec = math.prod(shape)
    return a[:a.size//rec*rec].reshape(-1, *shape), a.size % rec


def expected_times(case):
    p = namelist_numbers(case)
    dt, nstep, intvl = p['dt'], int(p['nstep']), p['write_intvl_3d']
    t, out = 0.0, [0.0]
    for _ in range(nstep):
        t = t + dt
        if math.floor(t/intvl) != math.floor((t + dt)/intvl):
            out.append(t)
    return out + [t]


# Series that equal, up to a constant factor, sum_k w_k |f_k|^2 over the
# named 3D fields (w_k: weight of the kz plane in the r2c layout).
ENERGY = (('u2_sum', ('ux', 'uy', 'uz')), ('b2_sum', ('bx', 'by', 'bz')),
          ('upa2_sum', ('upa',)), ('bpa2_sum', ('bpa',)))


def energy_groups(case, names, times):
    """For each ENERGY group whose fields are all in out3d/: compare the energy
    of the 3D records with the series at the record times. The group passes when
    E_3D/series is one constant over the records where the series is nonzero
    (rel. spread < 1e-10) and E_3D is exactly 0 where the series is 0.
    Returns {field: passed} and a message."""
    s = series(case)
    tt = s['tt']
    idx = [int(np.argmin(np.abs(tt - t))) for t in times]
    if any(abs(tt[i] - t) > 1e-12 for i, t in zip(idx, times)):
        return {}, 'energy check skipped: record times not in the series'
    nkz = kdims(case)[2]
    w = np.full(nkz, 2.0)
    w[0] = w[-1] = 1.0
    status, msg = {}, []
    for key, fields in ENERGY:
        if key not in s or not all(f'{f}.dat' in names for f in fields):
            continue
        e = sum(np.sum(w*np.abs(records(case, f'{f}.dat')[0])**2, axis=(1, 2, 3))
                for f in fields)
        ser = s[key][idx]
        if e.shape != ser.shape:
            passed, m = False, f'{e.size} records vs {ser.size} times'
        else:
            nz = ser != 0
            ratio = e[nz]/ser[nz]
            top = float(np.max(ratio)) if ratio.size else 0.0
            # (max - min)/max; a group whose records are all zero while the
            # series is not gets an infinite spread
            spread = ((top - float(np.min(ratio)))/top if top > 0
                      else (np.inf if ratio.size else 0.0))
            passed = bool(np.all(e[~nz] == 0)) and spread < 1e-10
            m = (f'E_3D/series = {top:.6g}, spread {spread:.1e}' if ratio.size else 'all zero')
            if np.any(~nz):
                m += f', series 0 at t = {[f"{t:.4f}" for t in np.array(times)[~nz]]}'
        msg.append(f'{key} [{" ".join(fields)}]: {m}')
        status.update({f'{f}.dat': passed for f in fields})
    return status, '; '.join(msg)


def r4_self(case):
    """(c) last record = restart file.
    (d) each out3d file is exactly len(time.dat) records of
        16*(nlz/2+1)*nly*nlx bytes; time.dat holds the expected times; no
        record is all zero unless the field is zero there, shown by the energy
        of its group in the series (ENERGY)."""
    if missing(case, 'R4'):
        return
    names = out3d_names(case)
    times = out3d_times(case)
    want = expected_times(case)
    rec_bytes = 16*math.prod(kdims(case))
    want_size = len(times)*rec_bytes
    sizes = {n: os.path.getsize(path(case, 'out3d', n)) for n in names}
    size_ok = bool(names) and all(v == want_size for v in sizes.values())
    t_ok = len(times) == len(want) and all(abs(a - b) <= 1e-14 for a, b in zip(times, want))
    last_ok, zero = [], {}
    for n in names:
        r, _ = records(case, n)
        if r.shape[0] and os.path.exists(path(case, 'restart', n)):
            last_ok.append(r[-1].tobytes() == raw(case, 'restart', n))
        z = [i for i, x in enumerate(r) if not np.any(x)]
        if z:
            zero[n] = z
    status, emsg = energy_groups(case, names, times) if size_ok else ({}, '')
    unexplained = [n for n in zero if not status.get(n, False)]
    egroups_ok = all(status.values())
    msg = (f'{len(names)} files, each {want_size} B = {len(times)} records x {rec_bytes} B'
           if size_ok else f'sizes {sizes}, want {want_size} B = {len(times)} x {rec_bytes} B')
    msg += f'; time.dat {[f"{t:.4f}" for t in times]}'
    if not t_ok:
        msg += f', expected {[f"{t:.4f}" for t in want]}'
    if zero:
        at = ((lambda z: f'@ t = {[f"{times[i]:.4f}" for i in z]}') if size_ok
              else (lambda z: f'@ record {z}'))
        msg += '; all-zero records ' + ', '.join(f'{n[:-4]} {at(z)}' for n, z in zero.items())
        msg += (' (field zero there)' if not unexplained else f' (NOT explained: {unexplained})')
    if emsg:
        msg += '; ' + emsg
    # Diagnose the layout written by a disp that advances two records per write.
    nrec = {v//rec_bytes for v in sizes.values()}
    if nrec == {2*len(times) - 1}:
        odd = all(not np.any(records(case, n)[0][1::2]) for n in names)
        if odd:
            msg += '; records 1, 3, ... are all-zero gaps: the data sit at a stride of 2 records'
    report(f'R4 (d) {case}', size_ok and t_ok and not unexplained and egroups_ok, msg)
    report(f'R4 (c) {case}', bool(last_ok) and all(last_ok) and len(last_ok) == len(names),
           f'last record == restart/<field>.dat for {sum(last_ok)}/{len(names)} fields')


def r4_pair(name, a, b):
    """(a), (b): out3d/*.dat and out3d/time.dat byte-identical. When the two
    solutions themselves differ (restart files not bitwise equal, as in a driven
    run on different decompositions), the 3D files can only agree to the same
    round-off; then the pass criterion is max rel < TOL."""
    if missing(a, name) or missing(b, name):
        return
    names = out3d_names(a) + ['time.dat']
    diff = differing(a, b, 'out3d', names)
    rdiff = differing(a, b, 'restart', restart_files(a))
    msg = 'identical' if not diff else f'differ: {diff}'
    worst = 0.0
    if diff:
        worst = 0.0
        for n in diff:
            if n == 'time.dat':
                continue
            ra, _ = records(a, n)
            rb, _ = records(b, n)
            if ra.shape != rb.shape:
                worst = float('inf')
                continue
            worst = max(worst, max_rel_arrays(ra, rb))
        msg += f', max rel {worst:.2e}'
        if rdiff:
            msg += (f'; the solutions differ too (restart max rel '
                    f'{max_rel_fields(a, b, restart_files(a)):.2e})')
    passed = not diff or (bool(rdiff) and 'time.dat' not in diff and worst < TOL)
    report(f'{name} {a} vs {b}', passed, f'{len(names)} files {msg}')


def r4_v10_pair(name, a, b):
    """v1.0 (a) vs v1.1 (b) on one GPU, same input. v1.0 advances the file
    offset by two records per write (fixed in a93e71c), so its record 2i must
    equal v1.1 record i bytewise and its odd records must be all zero;
    out3d/time.dat and the restart files must be identical."""
    if missing(a, name) or missing(b, name):
        return
    names = out3d_names(b)
    bad = [] if out3d_names(a) == names else [f'fields {out3d_names(a)} vs {names}']
    for n in names:
        ra, resta = records(a, n)
        rb, restb = records(b, n)
        if resta or restb or ra.shape[0] != 2*rb.shape[0] - 1:
            bad.append(f'{n}: {ra.shape[0]} (+{resta}) vs {rb.shape[0]} (+{restb}) records')
            continue
        if np.any(ra[1::2]):
            bad.append(f'{n}: an odd v1.0 record is nonzero')
        if ra[0::2].tobytes() != rb.tobytes():
            bad.append(f'{n}: data records differ, max rel {max_rel_arrays(ra[0::2], rb):.2e}')
    tdiff = raw(a, 'out3d', 'time.dat') != raw(b, 'out3d', 'time.dat')
    rdiff = differing(a, b, 'restart', restart_files(a))
    if tdiff:
        bad.append('time.dat differs')
    if rdiff:
        bad.append(f'restart differs: {rdiff}')
    nb = records(b, names[0])[0].shape[0] if names else 0
    report(f'{name} {a} vs {b}', not bad,
           f'{len(names)} fields: v1.0 records 0, 2, ..., {2*nb - 2} == v1.1 records 0..{nb - 1} '
           f'bytewise, v1.0 odd records all zero; time.dat and restart identical'
           if not bad else '; '.join(bad))


def bands(case, name, nband=4):
    """Which ky bands of the last 3D record equal the restart file."""
    r, _ = records(case, name)
    f = np.fromfile(path(case, 'restart', name), dtype=np.complex128).reshape(kdims(case))
    w = f.shape[1]//nband
    return [bool(np.array_equal(r[-1][:, j*w:(j + 1)*w, :], f[:, j*w:(j + 1)*w, :]))
            for j in range(nband)], r.shape[0]


def timer_3d(case):
    m = re.search(r'^\s*3D\s+([\d.]+) min\s+([\d.]+) %', out_std(case), re.M)
    return (float(m.group(1)), float(m.group(2))) if m else None


def set_r4():
    for c in cases('r4_'):
        if not missing(c, 'V0'):
            check_input(c)
            smoke(c)
    rc = [c for c in cases('r4_') if '_v10_' not in c]
    for c in rc:
        r4_self(c)
    for k in ('mi', 'rm', 'kf', 'hr', 'rr', 'em'):
        r4_pair('R4 (a) np1 vs np4', f'r4_{k}_np1', f'r4_{k}_np4')
    for k in ('rm', 'kf', 'hr', 'rr'):
        r4_pair('R4 (b) P_m = 1 vs 2', f'r4_{k}_np2', f'r4_{k}_pm2')
    # v1.0 with the 3D calls enabled
    for k in ('mi', 'rm'):
        for np_ in (1, 4):
            c = f'r4_{k}_v10_np{np_}'
            if missing(c, 'R4 v1.0'):
                continue
            info = []
            for n in out3d_names(c):
                b, nrec = bands(c, n)
                info.append(f'{n[:-4]}: {nrec} rec, bands {"".join("=" if x else "x" for x in b)}')
            report(f'R4 v1.0 {c}', None, '; '.join(info) + ' (last record vs restart by ky band)')
        r4_v10_pair('R4 v1.0 np1 vs v1.1 np1', f'r4_{k}_v10_np1', f'r4_{k}_np1')
    for c in ('r4_kf_np4', 'r4_mi_np4', 'r4_em_np1'):
        if missing(c, 'R4 cost'):
            continue
        size = sum(os.path.getsize(path(c, 'out3d', n)) for n in os.listdir(path(c, 'out3d')))
        t = timer_3d(c)
        report(f'R4 cost {c}', None,
               f'out3d {size/2**20:.1f} MiB; timer_io_3D '
               + (f'{t[0]*60:.2f} s ({t[1]:.1f} % of total)' if t else 'not found'))


# --------------------------------------------------------------------- R5 --
def set_r5():
    for c in cases('r5_'):
        if not missing(c, 'V0'):
            # dt of an A case is replaced by the dt_adjust file
            check_input(c, skip={'dt'} if '_A' in c else ())
            smoke(c)
    for k in ('mi', 'rm', 'kr', 'hr', 'em', 'rr'):
        for g in ('', '_g3'):
            a, b = f'r5_{k}_A{g}', f'r5_{k}_B{g}'
            if missing(a, 'R5'):
                continue
            std = out_std(a)
            gone = not os.path.exists(path(a, 'dt_adjust'))
            n = std.count('Manually dt is changed')
            report(f'R5 dt_adjust read {a}', n == 1 and gone,
                   f"message {n} times; file {'deleted' if gone else 'still there'}")
            compare('R5 A (dt 2e-3 -> dt_adjust 1e-3) vs B (dt 1e-3)', a, b, bitwise=True)
    for old, b in (('r5_mi_A_v10', 'r5_mi_B'), ('r5_rm_A_v10', 'r5_rm_B'),
                   ('r5_kr_A_f05', 'r5_kr_B')):
        compare('R5 old binary A vs B', old, b, bitwise=True, info=True)


if __name__ == '__main__':
    for s in SETS:
        print(f'===== {s}')
        {'r12': set_r12, 'r3': set_r3, 'r4': set_r4, 'r5': set_r5}[s]()
    print('ALL PASS' if ok else 'SOME CHECKS FAILED')
    sys.exit(0 if ok else 1)

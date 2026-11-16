#!/usr/bin/env python3
"""
Gate G3 of the MHD_COMP_ISOTH tests: compare a Gallope run with the same case
run by the final Calliope build s13d (outputs fetched from ps1).

  python3 compare_s13d.py RUN_G RUN_C [--tol 1e-11] [--tol-lines 2e-11] [--tol-final 1e-11]
                          [--tol-ts 1e-10]
                          [--glayout gallope|calliope] [--json OUT.json]

RUN_C is the s13d case directory (calliope.in, restart/, out1d/, cfl.dat, ...).
The grid is read from RUN_C/calliope.in. The Gallope restart is taken from
RUN_G/restart_g (Gallope layout, written by job-tsubame.sh) if it exists, else
from RUN_G/restart in the layout given by --glayout (default gallope; the two
layouts have the same size for a 1D grid, so it cannot be guessed).

Checks that decide the verdict:
  time     restart time, |t_g - t_c|/|t_c| < 1e-12 (dt sequences agree)
  restart  the seven spectral fields, expanded to the full Hermitian cube in
           both codes; combined relative L2 difference <= --tol
  out1d    every line_NNNNN.dat present in RUN_C: same file set, same rows, and the
           measure of Calliope's Stage 14 (analyze_stage3.maxdiff: max over rho, u_n,
           u_t1, u_t2, B_t1, B_t2, n the line direction, of max|a - b| / max(1, max|b|))
           <= --tol-lines on every line (D3 and CV2, against the older Fortran build:
           TOL_OLD1D = 2e-11) and <= --tol-final on the last line (D2, final state against
           the Python prototype: 1e-11; Calliope applies D2 to the suite and eta_suite
           cases only, compare_s13d_all.py passes inf for the others). Listed, not
           judged: the column-wise relative
           max_col max|a - b| / max|b| over all columns (zeta, D = div u and the
           components that vanish analytically included) except those at the rounding
           level (max|a| and max|b| both below ROUNDOFF = 1e-13)
Reported, not judged: cfl.dat, conservation.dat, budget.dat, out2d/plane_{y0,z0}.dat,
shock_selftest.dat (the rows must pair, see below; per-column normalized differences
are listed, and those above --tol-ts are flagged; rounding-level columns as for
out1d; out2d/plane_{y0,z0}.dat are binary, records of 10 real*8 arrays of (nlx, nlz) and (nlx, nly) values
whose t is in out2d/plane_time.dat, compared array by array over the paired records). The power-method
columns of cfl.dat (13-16) amplify round-off by the finite-difference Jacobian and are
listed separately. Comment lines must agree token by token (numbers within
--tol-ts).
A run continued from restarts keeps the outputs of its earlier pieces in
chunk_1/, chunk_2/, ... (job-tsubame.sh and the Calliope job script), and the
two codes split a run at different times. A piece stopped by its wall-time limit
writes a row of conservation.dat, a plane record and an out1d line at the stop
time, off the output schedule, and the next piece writes them again when it
starts. The time series and the plane records are therefore the concatenation of
chunk_1, chunk_2, ... and the case directory, a row whose t (column 1; column 2
of plane_time.dat) does not exceed the last one kept being dropped, and the rows
of the two runs are paired by t (|dt| <= 1e-12 max(1, |t|)); a row without a
partner must be the end of a piece. The out1d lines of all pieces are paired in
the same way by the t of their header. The record and step numbers (out1d
headers, plane_time.dat), counted from the start of each piece, are not compared.
Exit status 0 on PASS, 1 on FAIL.
"""
import argparse
import glob
import json
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from compare_calliope import full_cube  # noqa: E402

FIELDS = ['lnrho', 'mx', 'my', 'mz', 'bx', 'by', 'bz']
TS_FILES = ['cfl.dat', 'conservation.dat', 'budget.dat', 'shock_selftest.dat']
PLANE_BIN = ['out2d/plane_y0.dat', 'out2d/plane_z0.dat']   # binary: records of the arrays rho ... eta
PLANE_ARRAYS = 10
POWER_COLS = {'cfl.dat': [13, 14, 15, 16]}    # 1-based
ROUNDOFF = 1e-13                               # columns below it in both codes: rounding level
STEP = re.compile(r',?\s*step\s*=\s*\d+')
NUM = re.compile(r'[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?')


def grid_of(run_c):
    txt = open(os.path.join(run_c, 'calliope.in')).read()
    n = []
    for key in ('nlx', 'nly', 'nlz'):
        m = re.search(r'^\s*%s\s*=\s*(\d+)' % key, txt, re.M)
        if not m:
            sys.exit(f'{run_c}/calliope.in: no {key}')
        n.append(int(m.group(1)))
    return tuple(n)


def read_half(path, code, n):
    nlx, nly, nlz = n
    a = np.fromfile(path, dtype=np.complex128)
    if code == 'gallope':
        shp = (nlz//2 + 1, nly, nlx)
        return a.reshape(shp, order='F').transpose(2, 1, 0), 2
    shp = (nlx, nlz, nly//2 + 1)
    if a.size != np.prod(shp):
        sys.exit(f'{path}: {a.size} values, expected {np.prod(shp)} {shp}')
    return a.reshape(shp, order='F').transpose(0, 2, 1), 1


def restart_time(d):
    with open(os.path.join(d, 'time.dat')) as f:
        lines = [ln for ln in f.read().split('\n') if ln.strip()]
    return float(lines[-1].split()[0])


def segments(run):
    """The output directories of a run in time order: chunk_1, chunk_2, ... (earlier pieces of a run continued
    from restarts) and the run directory itself."""
    ks = sorted(int(m.group(1)) for m in (re.fullmatch(r'chunk_(\d+)', x) for x in os.listdir(run)) if m)
    return [os.path.join(run, f'chunk_{k}') for k in ks] + [run]


def split_file(path):
    """Comment lines and the numeric rows (as a 2D array) of an ASCII output."""
    com, rows = [], []
    with open(path) as f:
        for ln in f:
            if ln.lstrip().startswith('#'):
                com.append(ln.rstrip('\n'))
            elif ln.strip():
                rows.append([float(x.replace('D', 'E').replace('d', 'e')) for x in ln.split()])
    if rows and len({len(r) for r in rows}) != 1:
        return com, None
    return com, np.array(rows) if rows else np.zeros((0, 0))


def split_series(paths, tcol=0):
    """Comment lines of the first piece, the rows of all pieces of a time series in time order, and a flag per
    row for the end of a piece (the last row of a piece before the last, the first row of a piece after the
    first). A row whose t (column tcol + 1) does not exceed the last one kept is dropped (written again at the
    start of a piece) and passes its flag to the kept row."""
    com, parts, ends = None, [], []
    for k, path in enumerate(paths):
        c, a = split_file(path)
        if a is None:
            return c, None, None
        com = c if com is None else com
        if not a.size:
            continue
        if parts and parts[-1].shape[1] != a.shape[1]:
            return com, None, None
        e = np.zeros(len(a), bool)
        e[0] = k > 0
        e[-1] |= k < len(paths) - 1
        if parts:
            keep = a[:, tcol] > parts[-1][-1, tcol]
            ends[-1][-1] |= e[~keep].any()
            a, e = a[keep], e[keep]
        if a.size:
            parts.append(a)
            ends.append(e)
    if not parts:
        return com or [], np.zeros((0, 0)), np.zeros(0, bool)
    return com or [], np.concatenate(parts), np.concatenate(ends)


def pair_by_t(tg, tc):
    """Indices of the paired rows (|t_g - t_c| <= 1e-12 max(1, |t|)) of two increasing t sequences, and those of
    the rows without a partner in each."""
    i = j = 0
    ig, ic, ug, uc = [], [], [], []
    while i < len(tg) and j < len(tc):
        if abs(tg[i] - tc[j]) <= 1e-12*max(1.0, abs(tc[j])):
            ig.append(i)
            ic.append(j)
            i, j = i + 1, j + 1
        elif tg[i] < tc[j]:
            ug.append(i)
            i += 1
        else:
            uc.append(j)
            j += 1
    return ig, ic, ug + list(range(i, len(tg))), uc + list(range(j, len(tc)))


def pair_pieces(tg, eg, tc, ec, r):
    """Pair the rows of two runs by t (pair_by_t) and record the number of rows without a partner at the end of a
    piece in r; None and r['error'] set if a row without a partner is not at the end of a piece."""
    ig, ic, ug, uc = pair_by_t(tg, tc)
    lone = [('gallope', tg[i]) for i in ug if not eg[i]] + [('calliope', tc[j]) for j in uc if not ec[j]]
    if lone:
        r['error'] = f'{len(lone)} rows without a partner, first: {lone[0][0]} t = {lone[0][1]!r}'
        return None
    r['unpaired_ends'] = len(ug) + len(uc)
    return ig, ic


def line_time(path):
    """t in the header '# t = ..., dir = ..., step = ...' of an out1d line file."""
    with open(path) as f:
        m = re.match(r'#\s*t\s*=\s*(' + NUM.pattern + ')', f.readline())
    return float(m.group(1)) if m else None


def lines_of(run):
    """[(t, path, end_of_piece)] of the out1d lines of all pieces of a run in time order, one per t."""
    out, segs = [], segments(run)
    for i, seg in enumerate(segs):
        files = sorted(glob.glob(os.path.join(seg, 'out1d', 'line_*.dat')))
        for j, path in enumerate(files):
            t = line_time(path)
            end = (i > 0 and j == 0) or (i < len(segs) - 1 and j == len(files) - 1)
            if out and t is not None and out[-1][0] is not None and abs(t - out[-1][0]) <= 1e-12*max(1.0, abs(t)):
                out[-1] = (out[-1][0], out[-1][1], out[-1][2] or end)
                continue
            out.append((t, path, end))
    return out


def compare_comments(ca, cb, tol):
    """None if the comment lines agree token by token, else a short reason."""
    if len(ca) != len(cb):
        return f'{len(ca)} vs {len(cb)} comment lines'
    for a, b in zip(ca, cb):
        a, b = STEP.sub('', a), STEP.sub('', b)
        if NUM.sub('#', a) != NUM.sub('#', b):
            return f'comment differs: {a[:60]!r} vs {b[:60]!r}'
        na = [float(x) for x in NUM.findall(a)]
        nb = [float(x) for x in NUM.findall(b)]
        for x, y in zip(na, nb):
            if abs(x - y) > tol*max(abs(x), abs(y), 1.0):
                return f'comment number differs: {a[:60]!r} vs {b[:60]!r}'
    return None


def col_diffs(a, b):
    """Per-column max|a - b| / max|b| (absolute when the column of b is zero), the rounding-level columns
    (max|a| and max|b| both below ROUNDOFF: zero up to rounding in both codes, e.g. a residual or a momentum that
    vanishes by symmetry, whose relative difference means nothing) and their max|a - b|."""
    return col_ratio(np.abs(a - b).max(axis=0), np.abs(a).max(axis=0), np.abs(b).max(axis=0))


def col_ratio(d, sa, s):
    """col_diffs from the per-column max|a - b|, max|a| and max|b|."""
    rl = np.maximum(s, sa) < ROUNDOFF
    return np.where(s > 1e-300, d/np.where(s > 1e-300, s, 1.0), d), rl, d


def maxdiff_w(pg, pc):
    """Calliope's D2/D3 measure of one out1d line (analyze_stage3.maxdiff over the W of analyze_stage3.snaps_of):
    max over rho, u_n, u_t1, u_t2, B_t1, B_t2 (n: the direction in the header) of max|a - b| / max(1, max|b|)."""
    (_, ag), (cc, ac) = split_file(pg), split_file(pc)
    m = re.search(r'dir = ([xyz])', ' '.join(cc))
    i = 'xyz'.index(m.group(1)) if m else 0
    cols = [1, 2 + i, 2 + (i + 1) % 3, 2 + (i + 2) % 3, 5 + (i + 1) % 3, 5 + (i + 2) % 3]
    return max(float(np.abs(ag[:, j] - ac[:, j]).max())/max(1.0, float(np.abs(ac[:, j]).max())) for j in cols)


def compare_ascii(pg, pc, tol):
    """dict with rows, worst normalized column difference and the column, or an error. pg and pc are a path
    or a list of the paths of the pieces of a time series (rows paired by t, see the header)."""
    r = {'file': os.path.basename(pc[-1] if isinstance(pc, list) else pc)}
    if isinstance(pc, list):
        (cg, ag, eg), (cc, ac, ec) = split_series(pg), split_series(pc)
    else:
        (cg, ag), (cc, ac) = split_file(pg), split_file(pc)
    if ag is None or ac is None:
        r['error'] = 'ragged rows'
        return r
    if isinstance(pc, list) and ag.size and ac.size and ag.shape[1] == ac.shape[1]:
        p = pair_pieces(ag[:, 0], eg, ac[:, 0], ec, r)
        if p is None:
            return r
        ag, ac = ag[p[0]], ac[p[1]]
    if ag.shape != ac.shape:
        r['error'] = f'shape {ag.shape} vs {ac.shape}'
        return r
    r['rows'] = int(ac.shape[0])
    why = compare_comments(cg, cc, tol)
    if why:
        r['comment'] = why
    if ac.size == 0:
        r['worst'], r['col'] = 0.0, 0
        return r
    d, rl, dabs = col_diffs(ag, ac)
    r['cols'] = [float(x) for x in d]
    r['roundoff'] = [int(i) + 1 for i in np.flatnonzero(rl)]
    r['roundoff_abs'] = float(dabs[rl].max()) if rl.any() else 0.0
    dm = np.where(rl, 0.0, d)
    r['worst'], r['col'] = float(dm.max()), int(dm.argmax()) + 1
    if ac.shape[0] > 1:
        # rows 2..: the first row (t = 0) holds ratios of rounding-level quantities when the initial shock sensor
        # zeta vanishes analytically (div u = 0), e.g. cfl.dat 4 and 9, budget.dat 14
        d1, rl1, _ = col_diffs(ag[1:], ac[1:])
        r['cols_after1'] = [float(x) for x in d1]
        r['roundoff_after1'] = [int(i) + 1 for i in np.flatnonzero(rl1)]
    return r


def plane_series(run, name, npts):
    """t, end-of-piece flags and records (memory-mapped, PLANE_ARRAYS x npts each) of the plane file name
    (out2d/plane_y0.dat or plane_z0.dat) over the pieces of a run, as split_series with the t of
    out2d/plane_time.dat of each piece; a string on a size mismatch."""
    segs = [d for d in segments(run) if os.path.exists(os.path.join(d, name))]
    t, ends, recs = [], [], []
    for k, d in enumerate(segs):
        _, pt = split_file(os.path.join(d, 'out2d', 'plane_time.dat'))
        path = os.path.join(d, name)
        nrec = os.path.getsize(path)//(8*PLANE_ARRAYS*npts)
        if pt is None or nrec*8*PLANE_ARRAYS*npts != os.path.getsize(path) or nrec != len(pt):
            return f'{os.path.relpath(path, run)}: {os.path.getsize(path)} bytes, ' \
                   f'{0 if pt is None else len(pt)} records in plane_time.dat'
        if nrec == 0:
            continue
        v = np.memmap(path, dtype='<f8', mode='r', shape=(nrec, PLANE_ARRAYS, npts))
        for j in range(nrec):
            end = (k > 0 and j == 0) or (k < len(segs) - 1 and j == nrec - 1)
            if t and pt[j, 1] <= t[-1]:
                ends[-1] |= end
                continue
            t.append(pt[j, 1])
            ends.append(end)
            recs.append(v[j])
    return np.array(t), np.array(ends, bool), recs


def compare_plane(rg, rc, name, npts):
    """As compare_ascii for the binary plane file name of the runs rg and rc over their pieces (records of
    PLANE_ARRAYS real*8 arrays of npts values, the same layout in both codes, paired by t); the arrays take the
    place of the columns and the records of the rows."""
    r = {'file': os.path.basename(name)}
    sg, sc = plane_series(rg, name, npts), plane_series(rc, name, npts)
    for code, s in (('gallope', sg), ('calliope', sc)):
        if isinstance(s, str):
            r['error'] = f'{code} {s}'
            return r
    p = pair_pieces(sg[0], sg[1], sc[0], sc[1], r)
    if p is None:
        return r
    r['rows'] = len(p[1])
    if not p[1]:
        r['worst'], r['col'] = 0.0, 0
        return r
    d, sa, s = (np.zeros(PLANE_ARRAYS) for _ in range(3))
    for i, j in zip(*p):
        a, b = np.asarray(sg[2][i]), np.asarray(sc[2][j])
        d = np.maximum(d, np.abs(a - b).max(axis=1))
        sa = np.maximum(sa, np.abs(a).max(axis=1))
        s = np.maximum(s, np.abs(b).max(axis=1))
    d, rl, dabs = col_ratio(d, sa, s)
    r['cols'] = [float(x) for x in d]
    r['roundoff'] = [int(i) + 1 for i in np.flatnonzero(rl)]
    r['roundoff_abs'] = float(dabs[rl].max()) if rl.any() else 0.0
    dm = np.where(rl, 0.0, d)
    r['worst'], r['col'] = float(dm.max()), int(dm.argmax()) + 1
    return r


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('gallope')
    p.add_argument('calliope')
    p.add_argument('--tol', type=float, default=1e-11)
    p.add_argument('--tol-lines', type=float, default=2e-11)
    p.add_argument('--tol-final', type=float, default=1e-11)
    p.add_argument('--tol-ts', type=float, default=1e-10)
    p.add_argument('--glayout', choices=['gallope', 'calliope'], default='gallope')
    p.add_argument('--json')
    p.add_argument('--quiet', action='store_true')
    a = p.parse_args()
    n = grid_of(a.calliope)
    res = {'gallope': a.gallope, 'calliope': a.calliope, 'n': n, 'fail': []}
    say = (lambda *x: None) if a.quiet else print

    # time and restart fields
    gdir = os.path.join(a.gallope, 'restart_g')
    glay = 'gallope'
    if not os.path.isdir(gdir):
        gdir, glay = os.path.join(a.gallope, 'restart'), a.glayout
    cdir = os.path.join(a.calliope, 'restart')
    tg, tc = restart_time(gdir), restart_time(cdir)
    rel = abs(tg - tc)/max(abs(tc), 1e-300)
    res['time'] = {'gallope': tg, 'calliope': tc, 'rel': rel}
    if rel >= 1e-12:
        res['fail'].append('time')
    say(f'time      gallope {tg!r}  calliope {tc!r}  rel {rel:.1e}')
    num = den = 0.0
    res['fields'] = {}
    for f in FIELDS:
        hg, axg = read_half(os.path.join(gdir, f + '.dat'), glay, n)
        hc, axc = read_half(os.path.join(cdir, f + '.dat'), 'calliope', n)
        g, c = full_cube(hg, axg, n), full_cube(hc, axc, n)
        d, s = np.linalg.norm(g - c), np.linalg.norm(c)
        num, den = num + d*d, den + s*s
        res['fields'][f] = {'l2rel': d/s if s else float(d), 'norm': float(s)}
        say(f'field     {f:5s} L2 rel {d/s if s else d:.3e}   ||f|| {s:.6e}')
    comb = float(np.sqrt(num/den)) if den else 0.0
    res['restart_l2'] = comb
    if not comb <= a.tol:
        res['fail'].append('restart')
    say(f'restart   combined L2 rel {comb:.3e} (tol {a.tol:.0e})')

    # out1d lines, matched by t over all pieces
    lc, lg = lines_of(a.calliope), lines_of(a.gallope)
    worst, where, npair, ends = 0.0, '', 0, 0
    wrel, where_rel, final = 0.0, '', (-1.0, 0.0, '')
    rel = lambda p: os.path.relpath(p, a.calliope if p.startswith(a.calliope) else a.gallope)  # noqa: E731
    for code, mine, other in (('calliope', lc, lg), ('gallope', lg, lc)):
        for t, path, end in mine:
            if t is None:
                res['fail'].append('out1d header ' + rel(path))
                say(f'out1d     {rel(path)}: no t in the header')
            elif not any(u is not None and abs(t - u) <= 1e-12*max(1.0, abs(t)) for u, _, _ in other):
                if end:
                    ends += 1
                else:
                    res['fail'].append(f'out1d {code} t={t:.6g}')
                    say(f'out1d     {code} {rel(path)} (t = {t!r}) has no partner')
    for t, pc, _ in lc:
        pg = next((q for u, q, _ in lg if t is not None and u is not None and abs(t - u) <= 1e-12*max(1.0, abs(t))),
                  None)
        if pg is None:
            continue
        npair += 1
        r = compare_ascii(pg, pc, a.tol_ts)
        if 'error' in r:
            res['fail'].append('out1d ' + rel(pc))
            say(f'out1d     {rel(pc)}: {r["error"]}')
            continue
        if 'comment' in r:
            res['fail'].append('out1d header ' + rel(pc))
            say(f'out1d     {rel(pc)}: {r["comment"]}')
        if r['worst'] > wrel:
            wrel, where_rel = r['worst'], f'{rel(pc)} col {r["col"]}'
        w = maxdiff_w(pg, pc)
        if w > worst:
            worst, where = w, rel(pc)
        if t > final[0]:
            final = (t, w, rel(pc))
    res['out1d'] = {'n': npair, 'unpaired_ends': ends, 'worst': worst, 'where': where, 'final': final[1],
                    'where_final': final[2], 'worst_col_rel': wrel, 'where_col_rel': where_rel}
    if npair == 0 and lc:
        res['fail'].append('out1d none paired')
    if not worst <= a.tol_lines:
        res['fail'].append('out1d')
    if not final[1] <= a.tol_final:
        res['fail'].append('out1d final')
    note = f', {ends} unpaired ends of pieces' if ends else ''
    say(f'out1d     {npair} lines paired by t{note}, worst {worst:.3e} ({where}) (tol {a.tol_lines:.0e}), '
        f'last line {final[1]:.3e} (tol {a.tol_final:.0e}); '
        f'column-wise relative {wrel:.3e} ({where_rel}), not judged')

    # reported time series and planes
    res['ts'] = {}
    for name in TS_FILES + PLANE_BIN:
        if not os.path.exists(os.path.join(a.calliope, name)):
            continue
        pc = [q for q in (os.path.join(d, name) for d in segments(a.calliope)) if os.path.exists(q)]
        pg = [q for q in (os.path.join(d, name) for d in segments(a.gallope)) if os.path.exists(q)]
        if not pg:
            res['ts'][name] = {'error': 'missing in gallope'}
            res['fail'].append(name + ' missing')
            say(f'{name:16s} missing in the Gallope run')
            continue
        if name in PLANE_BIN:
            r = compare_plane(a.gallope, a.calliope, name, n[0]*n[2] if name.endswith('y0.dat') else n[0]*n[1])
        else:
            r = compare_ascii(pg, pc, a.tol_ts)
        res['ts'][name] = r
        if 'error' in r:
            res['fail'].append(name)
            say(f'{name:16s} {r["error"]}')
            continue
        pw = [c - 1 for c in POWER_COLS.get(name, []) if c <= len(r.get('cols', []))]
        ro = [c - 1 for c in r.get('roundoff', [])]
        rest = [(x, i + 1) for i, x in enumerate(r.get('cols', [])) if i not in pw and i not in ro]
        r['worst_main'], r['col_main'] = max(rest) if rest else (0.0, 0)
        r['worst_power'] = max((r['cols'][i] for i in pw if i not in ro), default=0.0)
        flag = ' *' if r['worst_main'] > a.tol_ts else ''
        extra = ''
        if 'cols_after1' in r:
            ro1 = [c - 1 for c in r['roundoff_after1']]
            rest1 = [(x, i + 1) for i, x in enumerate(r['cols_after1']) if i not in pw and i not in ro1]
            r['worst_main_after1'], r['col_main_after1'] = max(rest1) if rest1 else (0.0, 0)
            extra += f'  rows 2..: {r["worst_main_after1"]:.2e} (col {r["col_main_after1"]})'
        extra += f'  power-method cols {r["worst_power"]:.2e}' if pw else ''
        extra += f'  {r["unpaired_ends"]} unpaired ends of pieces' if r.get('unpaired_ends') else ''
        if ro:
            extra += f'  rounding-level cols {r["roundoff"]} max|diff| {r["roundoff_abs"]:.1e}'
        com = f'  [{r["comment"]}]' if 'comment' in r else ''
        say(f'{name:16s} {r["rows"]:6d} rows, worst {r["worst_main"]:.2e} (col {r["col_main"]}){flag}{extra}{com}')

    res['verdict'] = 'PASS' if not res['fail'] else 'FAIL'
    say(res['verdict'] + ('' if not res['fail'] else ': ' + ', '.join(res['fail'][:6])))
    if a.json:
        with open(a.json, 'w') as f:
            json.dump(res, f, indent=1)
    return 0 if not res['fail'] else 1


if __name__ == '__main__':
    sys.exit(main())

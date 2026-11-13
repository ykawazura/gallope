#!/usr/bin/env python3
"""
Deterministic Gallope <-> Calliope comparison of restart files (V8).

  python3 tests/compare_calliope.py --model HRMHD --gallope RUN_G --calliope RUN_C \
      --nlx 64 --nly 64 --nlz 32 [--tol 1e-10]

RUN_G/restart and RUN_C/restart hold raw complex(8) dumps of the spectral
fields. The two codes keep different half spaces of the same Hermitian cube:

  Gallope   Fortran (nkz, nky, nkx), nkz = nlz/2+1  (r2c along z)
  Calliope  Fortran (nkx, nkz, nky), nky = nly/2+1  (r2c along y; needs pruned = .false.)

Both are expanded to the full nlx x nly x nlz cube with F(-k) = conj(F(k))
and compared mode by mode; no FFT or normalization enters. Mode index m maps
to the signed wavenumber m for m <= n/2, else m - n, in both codes.

Gates, in order:
  1. restart times agree (rel < 1e-12), otherwise dt diverged;
  2. axis-mapping self-check on each reconstructed cube: omg + kprp^2 phi = 0
     for the potential models, k.b = 0 for EMHD (O(1) means a crossed axis);
  3. combined relative L2 difference of the independent fields < --tol.
Same method as gallope-calliope_comparison/gallope/determinism/compare.py,
extended to the HRMHD/RRMHD/EMHD field sets.
"""
import argparse
import os
import sys

import numpy as np

FIELDS = {'HRMHD': ['phi', 'omg', 'psi', 'upa', 'bpa'],
          'RRMHD': ['phi', 'omg', 'psi', 'upa', 'bpa'],
          'EMHD': ['bx', 'by', 'bz']}
# omg = -kprp^2 phi is redundant and is left out of the combined difference
COMBINED = {'HRMHD': ['phi', 'psi', 'upa', 'bpa'],
            'RRMHD': ['phi', 'psi', 'upa', 'bpa'],
            'EMHD': ['bx', 'by', 'bz']}


def signed(n):
    m = np.arange(n)
    return np.where(m <= n//2, m, m - n)


def read(path, shape_f):
    a = np.fromfile(path, dtype=np.complex128)
    if a.size != np.prod(shape_f):
        sys.exit(f'{path}: {a.size} values, expected {np.prod(shape_f)} {shape_f} '
                 '(Calliope must run with pruned = .false.)')
    return a.reshape(shape_f, order='F')


def full_cube(half, axis, n):
    """Fill the (nlx, nly, nlz) cube from a half space holding 0..n/2 on `axis`."""
    full = np.zeros(n, dtype=np.complex128)
    nh = n[axis]//2 + 1
    sl = [slice(None)]*3
    sl[axis] = slice(0, nh)
    full[tuple(sl)] = half
    src = [(-np.arange(m)) % m for m in n]
    src[axis] = (-np.arange(nh, n[axis])) % n[axis]
    sl[axis] = slice(nh, n[axis])
    full[tuple(sl)] = np.conj(half[np.ix_(*src)])
    return full


def load(run, code, name, n):
    nlx, nly, nlz = n
    p = os.path.join(run, 'restart', name + '.dat')
    if code == 'gallope':
        half = read(p, (nlz//2 + 1, nly, nlx)).transpose(2, 1, 0)   # -> [x, y, z]
        return full_cube(half, 2, n)
    half = read(p, (nlx, nlz, nly//2 + 1)).transpose(0, 2, 1)       # -> [x, y, z]
    return full_cube(half, 1, n)


def restart_time(run):
    with open(os.path.join(run, 'restart', 'time.dat')) as f:
        lines = [ln for ln in f.read().split('\n') if ln.strip()]
    return float(lines[-1].split()[0])


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--model', required=True, choices=list(FIELDS))
    p.add_argument('--gallope', required=True)
    p.add_argument('--calliope', required=True)
    p.add_argument('--nlx', type=int, required=True)
    p.add_argument('--nly', type=int, required=True)
    p.add_argument('--nlz', type=int, required=True)
    p.add_argument('--lx', type=float, default=1.0)
    p.add_argument('--ly', type=float, default=1.0)
    p.add_argument('--lz', type=float, default=1.0)
    p.add_argument('--tol', type=float, default=1e-10)
    a = p.parse_args()
    n = (a.nlx, a.nly, a.nlz)
    kx = (signed(a.nlx)/a.lx)[:, None, None]
    ky = (signed(a.nly)/a.ly)[None, :, None]
    kz = (signed(a.nlz)/a.lz)[None, None, :]
    ok = True

    tg, tc = restart_time(a.gallope), restart_time(a.calliope)
    rel = abs(tg - tc)/abs(tc)
    ok &= rel < 1e-12
    print(f'time      gallope {tg!r}  calliope {tc!r}  rel diff {rel:.1e}')

    g = {f: load(a.gallope, 'gallope', f, n) for f in FIELDS[a.model]}
    c = {f: load(a.calliope, 'calliope', f, n) for f in FIELDS[a.model]}

    for tag, cube in (('gallope', g), ('calliope', c)):
        if a.model == 'EMHD':
            r = np.abs(kx*cube['bx'] + ky*cube['by'] + kz*cube['bz']).sum() / \
                (np.sqrt(sum(np.abs(cube[f])**2 for f in FIELDS['EMHD'])).sum()
                 * np.sqrt(kx**2 + ky**2 + kz**2).sum())
            label = '<|k.b|>/(<|k|><|b|>)'
        else:
            kp2 = kx**2 + ky**2
            r = np.linalg.norm(cube['omg'] + kp2*cube['phi'])/np.linalg.norm(kp2*cube['phi'])
            label = '||omg + kprp^2 phi||/||kprp^2 phi||'
        ok &= r < 1e-12
        print(f'mapping   {tag:8s} {label} = {r:.1e}')

    num = den = 0.0
    for f in FIELDS[a.model]:
        d = np.linalg.norm(g[f] - c[f])
        s = np.linalg.norm(c[f])
        if f in COMBINED[a.model]:
            num += d*d
            den += s*s
        print(f'field     {f:4s} L2 rel {d/s if s else np.nan:.3e}   '
              f'max|df|/max|f| {np.abs(g[f] - c[f]).max()/(np.abs(c[f]).max() or 1):.3e}   '
              f'||f|| {s:.6e}')
    comb = np.sqrt(num/den)
    ok &= comb < a.tol
    print(f'combined  {",".join(COMBINED[a.model])}: {comb:.3e} (tol {a.tol:.0e})')
    print('PASS' if ok else 'FAIL')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())

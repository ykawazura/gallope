#!/usr/bin/env python3
"""
Convert MHD_COMP_ISOTH restart files between the Gallope and Calliope layouts.

  python3 gallope_layout.py to-calliope SRC DST --n NLX NLY NLZ [--check]
  python3 gallope_layout.py to-gallope  SRC DST --n NLX NLY NLZ [--check]
  python3 gallope_layout.py out3d-to-calliope DIR --n NLX NLY NLZ

SRC and DST are restart directories. Each field file (lnrho, mx, my, mz, bx, by,
bz) is a raw complex(8) dump of the spectral field, normalized by nlx*nly*nlz.
The two codes keep different halves of the same Hermitian cube:

  Gallope   Fortran (nkz, nky, nkx), nkz = nlz/2+1  (r2c along z)
  Calliope  Fortran (nkx, nkz, nky), nky = nly/2+1  (r2c along y, pruned = .false.)

The source half is expanded to the full cube with F(-k) = conj(F(k)) and the
other half is cut out of it, so the conversion only moves and conjugates
modes and adds no rounding. The planes that each layout stores twice (k = 0
and the Nyquist plane of its half axis) are taken as they are; a field that
is the r2c transform of real data is Hermitian on them, and c2r ignores the
anti-Hermitian part in both codes.

Other files in SRC (time.dat) are copied unchanged. With --check, both
layouts are transformed to real space with numpy (irfftn) and must agree to
1e-13 relative to max|f|, which verifies the axis mapping on the real data.

out3d-to-calliope converts the out3d/ files of a Gallope run (OUT3D_FIELDS,
records of the Gallope layout written back to back) in place to the Calliope
layout, record by record, so that the out3d readers of the Calliope analysis
scripts read them. out3d/time.dat is left as it is. Each file f.dat is first
renamed to f.dat.gallope, converted into f.dat.tmp, renamed to f.dat, listed
in DIR/.calliope_layout and its backup removed: a second call converts only
what is not listed yet (and redoes a file whose backup is still there), so it
can be repeated after an interruption. The run script removes .calliope_layout
before every run, whose new out3d/ files are in the Gallope layout again.
"""
import argparse
import os
import shutil
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from compare_calliope import full_cube  # noqa: E402

FIELDS = ['lnrho', 'mx', 'my', 'mz', 'bx', 'by', 'bz']
OUT3D_FIELDS = ['rho', 'mx', 'my', 'mz', 'bx', 'by', 'bz']
SENTINEL = '.calliope_layout'


def shape_of(code, n):
    nlx, nly, nlz = n
    return (nlz//2 + 1, nly, nlx) if code == 'gallope' else (nlx, nlz, nly//2 + 1)


def from_layout(a, code, n):
    """The half indexed [x, y, z] of a flat record in the Fortran layout of `code`."""
    a = a.reshape(shape_of(code, n), order='F')
    return a.transpose(2, 1, 0) if code == 'gallope' else a.transpose(0, 2, 1)


def to_layout(half, code):
    """Flat record in the Fortran layout of `code` of a half indexed [x, y, z]."""
    a = half.transpose(2, 1, 0) if code == 'gallope' else half.transpose(0, 2, 1)
    return np.asfortranarray(a).ravel(order='F').astype(np.complex128)


def read_half(path, code, n):
    """Return the stored half as an array indexed [x, y, z]."""
    shp = shape_of(code, n)
    a = np.fromfile(path, dtype=np.complex128)
    if a.size != np.prod(shp):
        sys.exit(f'{path}: {a.size} values, expected {np.prod(shp)} = {shp} ({code} layout)')
    return from_layout(a, code, n)


def write_half(path, half, code):
    """Write a half indexed [x, y, z] in the Fortran layout of `code`."""
    to_layout(half, code).tofile(path)


def half_axis(code):
    return 2 if code == 'gallope' else 1


def convert(half, src, dst, n):
    full = full_cube(half, half_axis(src), n)
    sl = [slice(None)]*3
    sl[half_axis(dst)] = slice(0, n[half_axis(dst)]//2 + 1)
    return full[tuple(sl)]


def to_real(half, code, n):
    """Real field [x, y, z] from a half indexed [x, y, z] (numpy irfftn on the half axis)."""
    nlx, nly, nlz = n
    if code == 'gallope':
        return np.fft.irfftn(half, s=(nlx, nly, nlz), axes=(0, 1, 2))*(nlx*nly*nlz)
    r = np.fft.irfftn(half.transpose(0, 2, 1), s=(nlx, nlz, nly), axes=(0, 1, 2))*(nlx*nly*nlz)
    return r.transpose(0, 2, 1)


def out3d_to_calliope(d, n):
    """Convert the out3d/ files of DIR d in place from the Gallope to the Calliope layout (see the docstring)."""
    rec = int(np.prod(shape_of('gallope', n)))
    sentinel = os.path.join(d, SENTINEL)
    done = set(open(sentinel).read().split()) if os.path.exists(sentinel) else set()
    tfile = os.path.join(d, 'time.dat')
    ntime = (sum(1 for ln in open(tfile) if ln.strip()) if os.path.exists(tfile) else None)
    status = 0
    for f in OUT3D_FIELDS:
        path = os.path.join(d, f + '.dat')
        bak, tmp = path + '.gallope', path + '.tmp'
        if os.path.exists(bak):
            src, what = bak, 'redone from ' + os.path.basename(bak)
        elif f in done:
            print(f'{path}: already in the Calliope layout')
            continue
        elif not os.path.exists(path):
            print(f'{path}: absent')
            continue
        else:
            src, what = path, 'converted'
        size = os.path.getsize(src)
        if size % (16*rec):
            print(f'{src}: {size} bytes is not a multiple of a record of {16*rec} bytes '
                  f'({shape_of("gallope", n)} complex(8)); left as it is', file=sys.stderr)
            status = 1
            continue
        if src == path:
            os.replace(path, bak)
        nrec = size//(16*rec)
        with open(bak, 'rb') as fi, open(tmp, 'wb') as fo:
            for _ in range(nrec):
                half = from_layout(np.fromfile(fi, dtype=np.complex128, count=rec), 'gallope', n)
                to_layout(convert(half, 'gallope', 'calliope', n), 'calliope').tofile(fo)
            fo.flush()
            os.fsync(fo.fileno())
        os.replace(tmp, path)
        with open(sentinel, 'a') as fs:
            fs.write(f + '\n')
            fs.flush()
            os.fsync(fs.fileno())
        os.remove(bak)
        note = '' if ntime is None or ntime == nrec else f' (WARNING: time.dat has {ntime} rows)'
        print(f'{path}: {what}, {nrec} records{note}')
    return status


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = p.add_subparsers(dest='direction')
    for name in ('to-calliope', 'to-gallope'):
        s = sp.add_parser(name)
        s.add_argument('src')
        s.add_argument('dst')
        s.add_argument('--n', type=int, nargs=3, required=True, metavar=('NLX', 'NLY', 'NLZ'))
        s.add_argument('--fields', default=','.join(FIELDS))
        s.add_argument('--check', action='store_true')
    s = sp.add_parser('out3d-to-calliope')
    s.add_argument('dir')
    s.add_argument('--n', type=int, nargs=3, required=True, metavar=('NLX', 'NLY', 'NLZ'))
    a = p.parse_args()
    if a.direction is None:
        p.error('give to-calliope, to-gallope or out3d-to-calliope')
    n = tuple(a.n)
    if a.direction == 'out3d-to-calliope':
        return out3d_to_calliope(a.dir, n)
    src, dst = ('gallope', 'calliope') if a.direction == 'to-calliope' else ('calliope', 'gallope')
    fields = a.fields.split(',')
    os.makedirs(a.dst, exist_ok=True)
    worst = 0.0
    for f in fields:
        h_src = read_half(os.path.join(a.src, f + '.dat'), src, n)
        h_dst = convert(h_src, src, dst, n)
        write_half(os.path.join(a.dst, f + '.dat'), h_dst, dst)
        if a.check:
            r_src, r_dst = to_real(h_src, src, n), to_real(h_dst, dst, n)
            scale = np.abs(r_src).max() or 1.0
            err = np.abs(r_src - r_dst).max()/scale
            worst = max(worst, err)
            print(f'{f:6s} max|f_{src} - f_{dst}|/max|f| = {err:.2e}')
    for name in sorted(os.listdir(a.src)):
        if name.endswith('.dat') and name[:-4] in fields:
            continue
        s = os.path.join(a.src, name)
        if os.path.isfile(s):
            shutil.copy2(s, os.path.join(a.dst, name))
    if a.check:
        ok = worst < 1e-13
        print('CHECK', 'PASS' if ok else 'FAIL', f'({worst:.2e})')
        return 0 if ok else 1
    return 0


if __name__ == '__main__':
    sys.exit(main())

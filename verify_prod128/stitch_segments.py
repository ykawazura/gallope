#!/usr/bin/env python3
"""
Stitch restart-chain segment NetCDFs into one continuous file.

The KRMHD code opens gallope.out.nc with NF90_CLOBBER, so each restart
segment writes a fresh file (nout resets to 1) while the physical time tt
continues (read from time.dat). job.pbs-echo128 preserves each segment as
echo128_seg_<jobid>.nc. This script concatenates all such files along the
unlimited 'tt' record dimension, sorts by tt, drops duplicate/overlapping
records (keeping the last occurrence), and writes a combined NetCDF that
analyze.py / plot_1Ciii.py / plot_em_spectrum.py read unchanged.

Usage:  python3 stitch_segments.py OUT.nc SEG1.nc SEG2.nc [SEG3.nc ...]
        python3 stitch_segments.py OUT.nc 'echo128_seg_*.nc'   (glob quoted)
"""
import sys
import glob
import numpy as np
from netCDF4 import Dataset


def expand(args):
    files = []
    for a in args:
        g = sorted(glob.glob(a))
        files.extend(g if g else [a])
    # de-dup while preserving order
    seen, out = set(), []
    for f in files:
        if f not in seen:
            seen.add(f)
            out.append(f)
    return out


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    out_path = sys.argv[1]
    seg_paths = expand(sys.argv[2:])
    if not seg_paths:
        sys.exit("no segment files matched")

    print("stitching %d segments -> %s" % (len(seg_paths), out_path))

    # Read tt from every segment; build the global sorted, de-duplicated axis.
    per_file_tt = []
    for p in seg_paths:
        with Dataset(p) as nc:
            tt = np.array(nc.variables['tt'][:], dtype=float)
        per_file_tt.append(tt)
        print("  %-40s  nrec=%4d  tt=[%.3f, %.3f]" %
              (p.split('/')[-1], tt.size, tt[0] if tt.size else np.nan,
               tt[-1] if tt.size else np.nan))

    # (tt_value, file_index, record_index) tuples; later files win on ties.
    triples = []
    for fi, tt in enumerate(per_file_tt):
        for ri, tv in enumerate(tt):
            triples.append((tv, fi, ri))
    # stable sort by tt; for equal tt keep the one from the later file
    triples.sort(key=lambda x: (x[0], x[1]))
    kept = []
    eps = 1e-9
    for tv, fi, ri in triples:
        if kept and abs(kept[-1][0] - tv) <= eps:
            kept[-1] = (tv, fi, ri)      # overwrite duplicate tt (later file)
        else:
            kept.append((tv, fi, ri))
    ntt_new = len(kept)
    print("  combined: %d unique records, tt=[%.3f, %.3f]" %
          (ntt_new, kept[0][0], kept[-1][0]))

    # Template = first segment (dims, variables, attributes).
    src0 = Dataset(seg_paths[0])
    dst = Dataset(out_path, 'w', format=src0.file_format)

    # dimensions (tt becomes fixed-size ntt_new; others copied verbatim)
    tt_dimname = src0.variables['tt'].dimensions[0]
    for name, dim in src0.dimensions.items():
        if name == tt_dimname:
            dst.createDimension(name, ntt_new)
        else:
            dst.createDimension(name, (None if dim.isunlimited() else len(dim)))

    # open all sources once for record gathering
    srcs = [Dataset(p) for p in seg_paths]

    for name, var in src0.variables.items():
        v = dst.createVariable(name, var.datatype, var.dimensions)
        v.setncatts({k: var.getncattr(k) for k in var.ncattrs()})
        if var.dimensions and var.dimensions[0] == tt_dimname:
            # time-dependent: gather record-by-record from the chosen files
            data = np.empty((ntt_new,) + var.shape[1:], dtype=var.dtype)
            for out_i, (tv, fi, ri) in enumerate(kept):
                data[out_i] = srcs[fi].variables[name][ri]
            v[:] = data
        else:
            v[:] = var[:]                # static: copy from template

    for a in src0.ncattrs():
        dst.setncattr(a, src0.getncattr(a))

    dst.close()
    for s in srcs:
        s.close()
    src0.close()
    print("wrote %s" % out_path)


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""
Calliope's checks D1, D4 (REPORT_stage14) and CV1 (REPORT_stage16) of the MHD1D folder on the Gallope runs.

  python3 check_mhd1d.py --calliope CALLIOPE_ROOT --runroot RUNROOT

Calliope's own commands (analyze_stage14.py mhd1d, analyze_stage16.py conv) also run D3 and CV2, which read the
runs of an earlier Calliope version, and stop where those runs are missing. This script calls the functions of
analyze_stage14.py that the two commands use for D1, D4 and CV1, with the same tolerances:

  D1   the 207 cases of stage14/mhd1d/{suite,eta_suite,decomp,eta_decomp} complete (lines_1d, complete_1d)
  CV1  the 8 cases of stage16/conv complete (the same functions)
  D4   every case of decomp and eta_decomp against its _x_np1 run: equal line steps and max difference
       <= TOL_DEC1D (decomp_1d)

The exit status is 0 if the three checks pass, 1 otherwise.
"""
import argparse
import glob
import os
import sys

N_1D, N_CONV = 207, 8
BATCHES = ('suite', 'eta_suite', 'decomp', 'eta_decomp')


def complete(a14, case):
    m, s = a14.lines_1d(case)
    return bool(s) and a14.complete_1d(case, m, s), len(s)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--calliope', required=True, help='Calliope repository (its tests/MHD_COMP_ISOTH/stages)')
    p.add_argument('--runroot', required=True)
    a = p.parse_args()
    stages = os.path.join(os.path.abspath(a.calliope), 'tests', 'MHD_COMP_ISOTH', 'stages')
    R = os.path.abspath(a.runroot)
    sys.path.insert(0, stages)
    os.chdir(stages)
    import analyze_stage14 as a14
    import analyze_stage3 as a3

    bad, n1, ok1 = [], 0, 0
    for b in BATCHES:
        cases = a3.cases_in(os.path.join(R, 'stage14', 'mhd1d', b))
        ok = 0
        for c in cases:
            comp, nl = complete(a14, c)
            ok += comp
            if not comp:
                bad.append('%s/%s (%d lines)' % (b, os.path.basename(c), nl))
        n1, ok1 = n1 + len(cases), ok1 + ok
        print('D1  %-11s %3d cases, %3d complete' % (b, len(cases), ok), flush=True)
    d1 = n1 == N_1D and ok1 == n1

    conv = [c for c in sorted(glob.glob(os.path.join(R, 'stage16', 'conv', '*'))) if os.path.isdir(c)]
    okc = 0
    for c in conv:
        comp, nl = complete(a14, c)
        okc += comp
        print('CV1 %-16s %6d lines, complete %s' % (os.path.basename(c), nl, comp), flush=True)
        if not comp:
            bad.append('conv/%s (%d lines)' % (os.path.basename(c), nl))
    cv1 = len(conv) == N_CONV and okc == N_CONV

    d4 = True
    for b in ('decomp', 'eta_decomp'):
        d = a14.decomp_1d(os.path.join(R, 'stage14', 'mhd1d', b))
        same = all(v['lines'][0] == v['lines'][1] for v in d.values())
        ok = bool(d) and same and all(v['max'] is not None and v['max'] <= a14.TOL_DEC1D for v in d.values())
        worst = max(d, key=lambda k: -1.0 if d[k]['max'] is None else d[k]['max']) if d else None
        print('D4  %-11s %d cases, max %s (%s), equal lines %s (tol %.0e): %s' % (
            b, len(d), 'n/a' if worst is None or d[worst]['max'] is None else '%.2e' % d[worst]['max'], worst,
            same, a14.TOL_DEC1D, 'PASS' if ok else 'FAIL'), flush=True)
        d4 = d4 and ok

    for x in bad:
        print('incomplete:', x)
    print(' '.join('%s: %s' % (k, 'PASS' if v else 'FAIL') for k, v in (('D1', d1), ('CV1', cv1), ('D4', d4))))
    return 0 if d1 and cv1 and d4 else 1


if __name__ == '__main__':
    sys.exit(main())

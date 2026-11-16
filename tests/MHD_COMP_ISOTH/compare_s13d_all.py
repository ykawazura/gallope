#!/usr/bin/env python3
"""
Run compare_s13d.py on a list of cases and tabulate the results.

  python3 compare_s13d_all.py run   CASEFILE --runroot RUNROOT --ref REF [--first I] [--count N]
  python3 compare_s13d_all.py table CASEFILE --runroot RUNROOT [--md]

CASEFILE lists one case per line (the first word, e.g. stage14/mhd1d/suite/C3_N256, as written by
make_cases.py list). `run` compares RUNROOT/<case> with REF/<case> (REF = the s13d tree of the folder, e.g.
.../s13d/MHD1D) for the cases first .. first+count-1 of the list (a login node limits a process to 10
minutes) and writes RUNROOT/<case>/compare_s13d.json. A case whose exit_status is not 'exit 0' is not
compared and is listed as such. The last-line check of out1d (Calliope's D2) applies to the cases of
stage14/mhd1d/suite and eta_suite only, as in Calliope (REPORT_stage14, section 3). `table` reads the JSON
files and prints one row per case and the worst value of every column over the cases; --md prints a markdown
table. The columns cfl, conservation and budget are the
worst normalized differences over the rows after the first (t > 0); 'all rows' is the worst of the three
including the first row, and 'planes' that of out2d/plane_{y0,z0}.dat.
"""
import argparse
import json
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
JSON = 'compare_s13d.json'


def cases_of(path):
    with open(path) as f:
        return [ln.split()[0] for ln in f if ln.strip() and not ln.lstrip().startswith('#')]


def status(run):
    p = os.path.join(run, 'exit_status')
    return open(p).read().strip() if os.path.exists(p) else '(none)'


def run(a):
    cases = cases_of(a.casefile)
    last = len(cases) if a.count is None else min(len(cases), a.first + a.count)
    for c in cases[a.first:last]:
        g, r = os.path.join(a.runroot, c), os.path.join(a.ref, c)
        st = status(g)
        if st != 'exit 0':
            print(f'{c}: not compared, exit_status {st}', flush=True)
            continue
        out = os.path.join(g, JSON)
        d2 = [] if re.match(r'stage14/mhd1d/(eta_)?suite/', c) else ['--tol-final', 'inf']
        p = subprocess.run([sys.executable, os.path.join(HERE, 'compare_s13d.py'), g, r, '--quiet', '--json', out]
                           + d2, capture_output=True, text=True)
        if p.returncode not in (0, 1) or not os.path.exists(out):
            print(f'{c}: compare_s13d.py failed ({p.returncode}): {(p.stderr or p.stdout).strip()[-300:]}', flush=True)
            continue
        print(f'{c}: {json.load(open(out))["verdict"]}', flush=True)
    return 0


def row(c, d):
    ts = d.get('ts', {})
    nan = float('nan')
    get = lambda name, key: ts.get(name, {}).get(key, nan)  # noqa: E731
    after1 = lambda name: ts.get(name, {}).get('worst_main_after1', get(name, 'worst_main'))  # noqa: E731
    main = ('cfl.dat', 'conservation.dat', 'budget.dat')
    allrows = max((get(x, 'worst_main') for x in main if x in ts), default=nan)
    planes = max((v.get('worst', nan) for k, v in ts.items() if k.startswith('out2d/plane_') and 'time' not in k),
                 default=nan)
    return [c, d['verdict'], d['time']['rel'], d['restart_l2'], d['out1d']['worst'], d['out1d'].get('final', nan),
            d['out1d'].get('worst_col_rel', nan), d['out1d']['n'],
            after1('cfl.dat'), get('cfl.dat', 'worst_power'), after1('conservation.dat'), after1('budget.dat'),
            allrows, planes, ';'.join(d['fail'][:3])]


# out1d, out1d last: Calliope's Stage 14 measure over all lines (D3, CV2) and on the last line (D2), judged;
# out1d col rel: column-wise relative difference, listed only
# cfl, conservation, budget: rows 2.. (t > 0); 'all rows' includes the first row (see compare_s13d.compare_ascii)
HEAD = ['case', 'verdict', 't rel', 'restart L2', 'out1d', 'out1d last', 'out1d col rel', 'lines', 'cfl',
        'cfl power', 'conservation', 'budget', 'all rows', 'planes', 'fail']


def table(a):
    rows, missing = [], []
    for c in cases_of(a.casefile):
        p = os.path.join(a.runroot, c, JSON)
        if os.path.exists(p):
            rows.append(row(c, json.load(open(p))))
        else:
            missing.append(f'{c} ({status(os.path.join(a.runroot, c))})')
    fmt = lambda x: f'{x:.1e}' if isinstance(x, float) else str(x)  # noqa: E731
    sep = ' | ' if a.md else '  '
    if a.md:
        print('| ' + ' | '.join(HEAD) + ' |')
        print('|' + '---|'*len(HEAD))
    for r in rows:
        line = sep.join(fmt(x) for x in r)
        print(f'| {line} |' if a.md else line)
    worst = ['worst', f'{sum(r[1] == "PASS" for r in rows)}/{len(rows)} PASS']
    for k in range(2, len(HEAD) - 1):
        if HEAD[k] == 'lines':
            worst.append(sum(r[k] for r in rows))
            continue
        vals = [r[k] for r in rows if r[k] == r[k]]          # drop nan (file absent)
        worst.append(max(vals) if vals else float('nan'))
    worst.append('')
    line = sep.join(fmt(x) for x in worst)
    print(f'| {line} |' if a.md else line)
    if missing:
        print(f'not compared ({len(missing)}): ' + ', '.join(missing))
    return 0 if rows and all(r[1] == 'PASS' for r in rows) and not missing else 1


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = p.add_subparsers(dest='cmd')
    s = sp.add_parser('run')
    s.add_argument('casefile')
    s.add_argument('--runroot', required=True)
    s.add_argument('--ref', required=True)
    s.add_argument('--first', type=int, default=0)
    s.add_argument('--count', type=int)
    s = sp.add_parser('table')
    s.add_argument('casefile')
    s.add_argument('--runroot', required=True)
    s.add_argument('--md', action='store_true')
    a = p.parse_args()
    if a.cmd is None:
        p.error('give run or table')
    return run(a) if a.cmd == 'run' else table(a)


if __name__ == '__main__':
    sys.exit(main())

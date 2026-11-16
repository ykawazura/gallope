#!/usr/bin/env python3
"""
Set up, check and list the MHD_COMP_ISOTH test cases of Calliope for runs with Gallope.

  make_cases.py setup    --calliope C --runroot R [--only GLOB ...]
  make_cases.py ic       --calliope C --runroot R [--python PY] [--only GLOB ...]
  make_cases.py check    --gallope G --runroot R [--only GLOB ...]
  make_cases.py list     --runroot R [--folder F ...] [--status S ...] [--max-ngpu N] [--long]
  make_cases.py stop     --runroot R
  make_cases.py folders  --calliope C --runroot R [--apply]
  make_cases.py ngpu-map --calliope C [--runroot R] [-o FILE]

C is the Calliope repository (its cases are C/tests/MHD_COMP_ISOTH/<Folder>/<stage>/<batch>/.../<case>, Folder
one of FOLDERS), R the run root and G the Gallope repository. A case runs in R/<stage>/<batch>/.../<case>, the
path of Calliope's runs/ (tests/MHD_COMP_ISOTH/runs, which can be a symbolic link to R, so that the analysis
scripts of Calliope find the runs at their default place). After `folders`, R/<stage>/.../<case> is a relative
symbolic link to the case in T/<Folder>/<stage>/.../<case> (T the parent of R), as in Calliope on ps1; every
command takes such links for cases.

setup     copies the input files of each case (INPUTS) to its run directory and writes
            ngpu           the number of GPUs (ranks) of the case, from ngpu_rule; an existing file is kept
            calliope_case  the path of the case in C/tests/MHD_COMP_ISOTH (with the Folder)
            calliope_extra the files of EXTRA that the case has (not copied: the job ignores them)
            terminate      for CHAIN_CASE, when the case has never run (the Calliope chain test stops at the
                           first step and continues from restart/, Euler1D/README.md)
          A file that exists is never overwritten; a different one is reported as a conflict (exit status 1).
          calliope.in and nproc are copied unchanged.
ic        writes profile1d.dat or ic/ of the cases with generated.md5 with stages/regen_inputs.py of C and
          reports the cases whose files differ from generated.md5 (exit status 1).
check     translates every calliope.in as the job does (job-tsubame.sh translate) and compares the groups and
          keys with the namelists declared in G/src/*.F90 and G/src/model/MHD_COMP_ISOTH/*.F90. Gallope reads
          every group with iostat, so an undeclared key silently leaves the whole group at its defaults and an
          undeclared group is silently ignored. Exit status 1 on any undeclared group or key.
list      prints the run directories (relative to R) by decreasing cost nlx*nly*nlz*nstep. --status takes
          'all' (default), 'todo' (what the job does: no exit_status or "incomplete", or an out3d/ conversion
          that the end of a job interrupted), 'done', or a prefix of exit_status (e.g. 'exit 0', 'dt_changed',
          'exit layout').
stop      sets t_end (calliope.in and the file t_end) of the cases of CALLIOPE_STOP, which Calliope stopped
          before t_end, to the t at which Calliope stopped; the nominal files go to calliope.in.nominal and
          t_end.nominal (setup then compares the Calliope files with these). Only cases that have never run
          are changed; exit status 1 on a case that has run or was stopped at another t.
folders   moves each case R/<run> that is a directory to T/<calliope_case> (R must be T/runs), leaves the
          relative link R/<run> -> T/<calliope_case>, and makes in T the symbolic links of the folders of C (the
          summaries in runs/ that the figure scripts read). Dry run unless --apply; links already made and cases
          already moved are skipped, a path that exists otherwise is a conflict (exit status 1, nothing moved).
          --apply writes T/attic/folders_<time>/{moves.txt,links.txt,rollback.sh}; `DRY=0 sh rollback.sh` undoes it.
          The moves are renames, which fail across file systems (or Lustre projects) instead of copying.
ngpu-map  writes the table of ngpu_rule (Markdown).
"""
import argparse
import collections
import filecmp
import fnmatch
import glob
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
import time

FOLDERS = ['Euler1D', 'MHD1D', 'Smooth_wave', 'OT2D', 'Turb3D']
INPUTS = ['calliope.in', 'nproc', 'case.json', 'generated.md5', 't_end']
EXTRA = ['exe', 'fftw_estimate']
CHAIN_CASE = 'stage1_es/smoke/collide_U2_N512_est_chain'
MAX_GPU = 4                     # GPUs of a TSUBAME4.0 node
HERE = os.path.dirname(os.path.abspath(__file__))
OUT3D_FIELDS = ['rho', 'mx', 'my', 'mz', 'bx', 'by', 'bz']     # gallope_layout.OUT3D_FIELDS (no numpy needed here)

# Decomposition tests: Calliope compares the members of a family with a reference run. The reference runs on
# 1 GPU; the other members run on 2, 4, ... GPUs in the order of their Calliope np (the j-th distinct np gets
# 2^j), capped by MAX_GPU and by the slab decomposition of Gallope (the number of ranks divides nlx and nly).
# Each entry: regular expression on the run path (groups fam, optional ax, np, suf), np and axis of the reference
# (np None: no _np suffix), and where the comparison is made.
DECOMP = [
    (r'(?P<fam>stage1_es/decomp/U2_N512)_(?P<ax>[xyz])_np(?P<np>\d+)', '1', 'x',
     'judge_stage1.judge_batch_decomp'),
    (r'(?P<fam>stage1_es/smoke/selftest)_np(?P<np>\d+)', '1', None, 'judge_stage1 (shock_selftest.dat)'),
    (r'(?P<fam>stage14/mhd1d/(?:eta_)?decomp/[A-Za-z0-9]+(?:_w\d+)?_N256(?:_eta2)?)_(?P<ax>[xyz])_np(?P<np>\d+)',
     '1', 'x', 'analyze_stage14 (TOL_DEC1D)'),
    (r'(?P<fam>stage4/alfven/(?:2d_k101_a1_N32_default_p800|3d_k112_a1_N16_default_p400))(?:_np(?P<np>\d+))?',
     None, None, 'analyze_stage4 (TOL_MPI)'),
    (r'(?P<fam>stage14/ot/dec_N256)_np(?P<np>\d+)(?P<suf>_n10)?', '16', None, 'analyze_stage4 (ot_decomp)'),
    (r'(?P<fam>stage6q/track/dec_N64)_np(?P<np>\d+)', '64', None, 'analyze_stage6 (TDEC_REF)'),
    (r'(?P<fam>stage15/pre/dec_turb1_N64_[HL]_r1)_np(?P<np>\d+)', '64', None, 'analyze_stage15 (P0B_NP)'),
]

# Cases that Calliope stopped before t_end (dt collapsed to 1e-8 .. 3e-13, or stopped for cost): t of
# restart/time.dat of the last Calliope run (ps1, calliope_isoth/tests/MHD_COMP_ISOTH/runs) and the reason
# (exit_status of that run). `stop` sets t_end of calliope.in and the file t_end to this t, so that Gallope
# covers the same interval; the nominal files are kept as *.nominal.
CALLIOPE_STOP = {
    'stage4/turb/turb4_N64_a':    (0.12725048142419873, 'dt collapse at t = 0.1273 (dt 1e-8), REPORT_stage4 8.3'),
    'stage4/turb/turb4b_N64':     (0.18346570621299313, 'incomplete (dt 2e-8)'),
    'stage5/res/turb4_N128_a':    (0.11111480706190269, 'incomplete (dt 4e-8)'),
    'stage5/res/turb4b_N128':     (0.21627417443002028, 'incomplete (dt 4e-7)'),
    'stage5/ms/turb3b_N64':       (0.42811413524584641, 'dt 3e-6 = dt0/67, REPORT_stage5 12'),
    'stage5/ms/turb2b_N64':       (0.96222620829256733, 'dt 3e-6 = dt0/171, REPORT_stage5 12'),
    'stage13d/nomass/e0.5_C2':    (1.4642941010717285, 'incomplete (dt 2.8e-13)'),
    'stage15/b/turb1_N128_ESf':   (2.0518317232605416, 'stopped for cost (dt 7.7e-12)'),
    'stage15/a/turb1_N64_ESf':    (2.9972186100924301, 'stopped for cost (dt 3.5e-11)'),
}


# ----------------------------------------------------------------------------------------------- helpers
class Case:
    def __init__(self, tpath, run, src):
        self.tpath = tpath      # path in tests/MHD_COMP_ISOTH, with the Folder
        self.run = run          # path in the run root, without the Folder
        self.src = src          # absolute path of the Calliope case
        self.folder = tpath.split('/')[0]


def tests_dir(calliope):
    return os.path.join(os.path.abspath(calliope), 'tests', 'MHD_COMP_ISOTH')


def calliope_cases(calliope, only=None):
    """The cases of the per-test folders of Calliope (directories with calliope.in), sorted by run path."""
    T = tests_dir(calliope)
    cases = []
    for F in FOLDERS:
        top = os.path.join(T, F)
        if not os.path.isdir(top):
            sys.exit('%s: no folder %s' % (T, F))
        for root, dirs, files in os.walk(top):
            dirs.sort()
            if 'calliope.in' in files:
                dirs[:] = []
                cases.append(Case(os.path.relpath(root, T), os.path.relpath(root, top), root))
    seen = collections.Counter(c.run for c in cases)
    dup = sorted(r for r, k in seen.items() if k > 1)
    if dup:
        sys.exit('run paths of several folders coincide: %s' % ' '.join(dup))
    cases.sort(key=lambda c: c.run)
    return [c for c in cases if selected(c.run, only, c.tpath)]


def selected(run, only, tpath=None):
    if not only:
        return True
    return any(fnmatch.fnmatchcase(run, g) or (tpath and fnmatch.fnmatchcase(tpath, g)) for g in only)


def run_dirs(runroot, only=None):
    """Run directories under the run root, as paths relative to it: directories with calliope.in and symbolic
    links to such directories (the links that `folders` leaves). Other links are not followed."""
    R = os.path.abspath(runroot)
    out = []
    for root, dirs, files in os.walk(R):
        dirs.sort()
        if 'calliope.in' in files:
            dirs[:] = []
            found = [root]
        else:
            dirs[:] = [d for d in dirs if not d.startswith('chunk_')]
            found = [os.path.join(root, d) for d in dirs
                     if os.path.islink(os.path.join(root, d)) and os.path.isfile(os.path.join(root, d, 'calliope.in'))]
        for d in found:
            run = os.path.relpath(d, R)
            if selected(run, only, read(os.path.join(d, 'calliope_case'))):
                out.append(run)
    return sorted(out)


def read(path, default=''):
    try:
        with open(path) as f:
            return f.read().strip()
    except OSError:
        return default


def nml_value(path, key):
    """First value of 'key = value' in a namelist file (Fortran d exponent read as e), or None."""
    m = re.search(r'^\s*%s\s*=\s*([^\s,!/]+)' % key, read(path), re.M | re.I)
    return m.group(1).replace('d', 'e').replace('D', 'e') if m else None


def nml_int(path, key, default=0):
    v = nml_value(path, key)
    return int(float(v)) if v is not None else default


def gpu_cap(nlx, nly):
    """Largest power of two <= MAX_GPU that divides nlx and nly."""
    g, p = math.gcd(nlx, nly), 1
    while 2*p <= MAX_GPU and g % (2*p) == 0:
        p *= 2
    return p


# ------------------------------------------------------------------------------------------------ ngpu rule
def ngpu_rule(cases):
    """{run path: (ngpu, Calliope np, reason)} of the cases."""
    info = {c.run: (int(read(os.path.join(c.src, 'nproc'), '1')),
                    nml_int(os.path.join(c.src, 'calliope.in'), 'nlx', 1),
                    nml_int(os.path.join(c.src, 'calliope.in'), 'nly', 1)) for c in cases}
    out = {r: (1, info[r][0], 'single GPU' + (' (Calliope np > 1 is not a decomposition test)'
                                              if info[r][0] > 1 else '')) for r in info}
    fams = collections.OrderedDict()
    for c in cases:
        for pat, ref_np, ref_ax, where in DECOMP:
            m = re.fullmatch(pat, c.run)
            if m:
                key = m.group('fam') + (m.groupdict().get('suf') or '')
                ref = m.group('np') == ref_np and (ref_ax is None or m.group('ax') == ref_ax)
                fams.setdefault(key, (where, []))[1].append((c.run, ref))
                break
    for fam, (where, members) in fams.items():
        refs = [r for r, is_ref in members if is_ref]
        others = [r for r, is_ref in members if not is_ref]
        if len(refs) != 1 or not others:
            sys.exit('decomposition family %s: %d references, %d other members' % (fam, len(refs), len(others)))
        out[refs[0]] = (1, info[refs[0]][0], 'reference of decomposition test %s, %s' % (fam, where))
        nps = sorted({info[r][0] for r in others})
        for r in others:
            np_, nlx, nly = info[r]
            j = nps.index(np_) + 1
            n = min(2**j, gpu_cap(nlx, nly))
            why = 'decomposition test %s: np %d is no. %d of np %s -> 2^%d' % (
                fam, np_, j, '/'.join(map(str, nps)), j)
            if n < 2**j:
                why += (', capped at MAX_GPU = %d' % n if n == MAX_GPU else
                        ', capped at %d (the ranks divide nlx = %d and nly = %d)' % (n, nlx, nly))
            out[r] = (n, np_, why)
    return out


# ------------------------------------------------------------------------------------------------------ setup
def place(src, dst, report):
    """Copy src to dst unless dst exists; a dst that differs is a conflict."""
    if os.path.exists(dst):
        if os.path.exists(dst + '.nominal'):        # edited by `stop`
            dst += '.nominal'
        if not filecmp.cmp(src, dst, shallow=False):
            report['conflict'].append('%s differs from %s' % (dst, src))
        else:
            report['same'] += 1
        return
    shutil.copyfile(src, dst)
    report['copied'] += 1


def put(text, dst, report, keep=False):
    """Write text to dst unless dst exists (keep: an existing dst wins, it is no conflict)."""
    if os.path.exists(dst):
        if read(dst) != text.strip():
            if keep:
                report['kept'].append('%s: %s (rule: %s)' % (dst, read(dst), text.strip()))
            else:
                report['conflict'].append('%s differs from %r' % (dst, text.strip()))
        else:
            report['same'] += 1
        return
    with open(dst, 'w') as f:
        f.write(text)
    report['copied'] += 1


def out3d_pending(d, st):
    """True if the out3d/ conversion after the last run of d (exit_status st) is not finished, as in out3d_pending
    of job-tsubame.sh: a backup f.dat.gallope or an f.dat not listed in out3d/.calliope_layout."""
    if not st or st.startswith('exit layout'):
        return False
    o3 = os.path.join(d, 'out3d')
    done = read(os.path.join(o3, '.calliope_layout')).split()
    return any(os.path.exists(os.path.join(o3, f + '.dat.gallope'))
               or (os.path.exists(os.path.join(o3, f + '.dat')) and f not in done) for f in OUT3D_FIELDS)


def never_ran(d):
    return not any(os.path.exists(os.path.join(d, f)) for f in ('exit_status', 'calliope.out.std', 'chunk_1'))


def cmd_setup(a):
    allc = calliope_cases(a.calliope)
    rule = ngpu_rule(allc)
    cases = [c for c in allc if selected(c.run, a.only, c.tpath)]
    R = os.path.abspath(a.runroot)
    report = {'copied': 0, 'same': 0, 'conflict': [], 'kept': [], 'new': 0, 'terminate': []}
    for c in cases:
        d = os.path.join(R, c.run)
        if not os.path.isdir(d):
            os.makedirs(d)
            report['new'] += 1
        for f in INPUTS:
            if os.path.exists(os.path.join(c.src, f)):
                place(os.path.join(c.src, f), os.path.join(d, f), report)
        extra = []
        for f in EXTRA:
            if os.path.exists(os.path.join(c.src, f)):
                body = read(os.path.join(c.src, f))
                extra.append(f + (' ' + body if body else ''))
        if extra:
            put('\n'.join(extra) + '\n', os.path.join(d, 'calliope_extra'), report)
        put(c.tpath + '\n', os.path.join(d, 'calliope_case'), report)
        put('%d\n' % rule[c.run][0], os.path.join(d, 'ngpu'), report, keep=True)
        if c.run == CHAIN_CASE:
            t = os.path.join(d, 'terminate')
            if never_ran(d) and not os.path.exists(t):
                open(t, 'w').close()
                report['terminate'].append(t)
    print('setup: %d cases (%s) in %s: %d new directories, %d files written, %d files already identical'
          % (len(cases), ', '.join('%s %d' % (F, sum(c.folder == F for c in cases)) for F in FOLDERS),
             R, report['new'], report['copied'], report['same']))
    for t in report['terminate']:
        print('terminate created: %s' % t)
    for k in report['kept']:
        print('ngpu file kept: %s' % k)
    for k in report['conflict']:
        print('CONFLICT: %s' % k)
    return 1 if report['conflict'] else 0


# --------------------------------------------------------------------------------------------------------- ic
def cmd_ic(a):
    S = os.path.join(tests_dir(a.calliope), 'stages')
    R = os.path.abspath(a.runroot)
    runs = [r for r in run_dirs(R, a.only) if os.path.exists(os.path.join(R, r, 'generated.md5'))]
    pat = re.compile(r'^(.*): (wrote|kept|nothing to write), (md5 OK|missing or different: (.*))$')
    res, errors = {}, []
    for i in range(0, len(runs), a.batch):
        batch = [os.path.join(R, r) for r in runs[i:i + a.batch]]
        # -B: no __pycache__ written into the Calliope tree
        p = subprocess.Popen([a.python, '-B', os.path.join(S, 'regen_inputs.py')] + batch, cwd=os.path.abspath(a.calliope),
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        out, err = p.communicate()
        for ln in out.splitlines():
            m = pat.match(ln)
            if m:
                res[m.group(1)] = (m.group(2), m.group(4).split() if m.group(4) else [])
                print(os.path.relpath(m.group(1), R) + ': ' + ln[len(m.group(1)) + 2:], flush=True)
        missing = [b for b in batch if b not in res]
        if missing:
            tail = [ln for ln in err.splitlines() if ln.strip()][-3:]
            errors += ['%s: no result (regen_inputs.py exit %d: %s)' % (os.path.relpath(b, R), p.returncode,
                                                                          ' | '.join(tail)) for b in missing]
    bad = sorted(os.path.relpath(c, R) for c, (_, names) in res.items() if names)
    what = collections.Counter(w for w, _ in res.values())
    print('ic: %d cases with generated.md5: %s; md5 OK %d, different %d, errors %d'
          % (len(runs), ', '.join('%s %d' % kv for kv in sorted(what.items())), len(res) - len(bad), len(bad),
             len(errors)))
    for b in bad:
        print('DIFFERENT: %s: %s' % (b, ' '.join(res[os.path.join(R, b)][1])))
    for e in errors:
        print('ERROR: %s' % e)
    return 1 if bad or errors else 0


# ------------------------------------------------------------------------------------------------------ check
def strip_comment(line):
    """Code part of a namelist line (up to the first ! outside quotes)."""
    q = None
    for i, ch in enumerate(line):
        if q:
            if ch == q:
                q = None
        elif ch in '\'"':
            q = ch
        elif ch == '!':
            return line[:i]
    return line


def gallope_namelists(gallope):
    """{group: [(file:line, set of keys)]} of the namelists of the MHD_COMP_ISOTH build of Gallope."""
    files = sorted(glob.glob(os.path.join(gallope, 'src', '*.F90'))
                   + glob.glob(os.path.join(gallope, 'src', 'model', 'MHD_COMP_ISOTH', '*.F90')))
    out = collections.defaultdict(list)
    for f in files:
        lines = open(f, errors='replace').read().split('\n')
        i = 0
        while i < len(lines):
            m = re.match(r'\s*namelist\s*/(\w+)/\s*(.*)$', strip_comment(lines[i]), re.I)
            if m:
                lineno, body = i + 1, m.group(2).rstrip()
                while body.endswith('&') and i + 1 < len(lines):
                    i += 1
                    body = body[:-1] + ' ' + strip_comment(lines[i]).strip().lstrip('&').rstrip()
                keys = {k.strip().lower() for k in body.split(',') if k.strip()}
                out[m.group(1).lower()].append(('%s:%d' % (os.path.relpath(f, gallope), lineno), keys))
            i += 1
    return out


def input_groups(text):
    """[(group, set of keys, {key: value text})] of a namelist file."""
    out, grp = [], None
    for line in text.split('\n'):
        code = strip_comment(line).strip()
        if grp is None and code.startswith('&'):
            grp = (code[1:].split() or [''])[0].lower()
            out.append((grp, set(), {}))
            continue
        if grp is not None and code.startswith('/'):
            grp = None
            continue
        if grp is not None:
            for m in re.finditer(r'(\w+)\s*(\(\s*[\d:,\s]*\))?\s*=\s*(\'[^\']*\'|"[^"]*"|[^\s,/]+)?', code):
                out[-1][1].add(m.group(1).lower())
                out[-1][2][m.group(1).lower()] = m.group(3)
    return out


def cmd_check(a):
    G = os.path.abspath(a.gallope)
    job = os.path.join(G, 'tests', 'MHD_COMP_ISOTH', 'job-tsubame.sh')
    R = os.path.abspath(a.runroot)
    decl = gallope_namelists(G)
    runs = run_dirs(R, a.only)
    use = collections.defaultdict(collections.Counter)
    bad_group, bad_key, mixed, other = (collections.defaultdict(list) for _ in range(4))
    for r in runs:
        p = subprocess.run(['bash', job, 'translate', os.path.join(R, r, 'calliope.in')],
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        if p.returncode or p.stderr.strip():
            other['translate failed: ' + p.stderr.strip()[:200]].append(r)
            continue
        groups = input_groups(p.stdout)
        names = [g for g, _, _ in groups]
        for g in set(names):
            if names.count(g) > 1:
                other['group &%s appears %d times' % (g, names.count(g))].append(r)
        if 'mpi_settings' in names:
            other['&mpi_settings left in'].append(r)
        rd = [v.get('restart_dir') for g, _, v in groups if g == 'restart_parameters']
        if rd != ["'./restart_g/'"]:
            other['restart_dir is %s' % rd].append(r)
        if not any('max_wall_time' in k for g, k, _ in groups if g == 'terminate'):
            other['no max_wall_time in &terminate'].append(r)
        for g, keys, _ in groups:
            use[g]['__group__'] += 1
            if 'pruned' in keys:
                other['pruned left in &%s' % g].append(r)
            if g not in decl:
                bad_group[g].append(r)
                continue
            union = set().union(*(k for _, k in decl[g]))
            for k in sorted(keys - union):
                bad_key[(g, k)].append(r)
            if keys <= union and not any(keys <= k for _, k in decl[g]):
                mixed[g].append(r)
    print('check: %d run directories, Gallope namelists of %s (src/*.F90, src/model/MHD_COMP_ISOTH/*.F90)'
          % (len(runs), G))
    print('groups of the translated inputs (inputs using it):')
    for g in sorted(use):
        where = '; '.join('%s' % w for w, _ in decl.get(g, [])) or 'NOT DECLARED'
        print('  &%-24s %4d   %s' % (g, use[g]['__group__'], where))

    def show(title, d):
        for k in sorted(d):
            ex = d[k]
            print('%s %s: %d inputs, e.g. %s' % (title, k if isinstance(k, str) else '&%s %s' % k, len(ex),
                                                 ', '.join(ex[:3])))
    show('UNDECLARED GROUP', bad_group)
    show('UNDECLARED KEY', bad_key)
    show('KEYS OF SEVERAL DECLARATIONS', mixed)
    show('TRANSLATION', other)
    n = len(bad_group) + len(bad_key) + len(mixed) + len(other)
    print('check: %s' % ('PASS (every group and key is declared)' if n == 0 else 'FAIL (%d problems)' % n))
    return 1 if n else 0


# ------------------------------------------------------------------------------------------------------- list
def cmd_list(a):
    R = os.path.abspath(a.runroot)
    rows = []
    for r in run_dirs(R):
        d = os.path.join(R, r)
        folder = read(os.path.join(d, 'calliope_case')).split('/')[0] or '?'
        if a.folder and folder not in a.folder:
            continue
        st = read(os.path.join(d, 'exit_status'))
        pending = out3d_pending(d, st)
        if not any(want == 'all' or (want == 'todo' and (st in ('', 'incomplete') or pending))
                   or (want == 'done' and st not in ('', 'incomplete'))
                   or (want not in ('all', 'todo', 'done') and st.startswith(want)) for want in a.status):
            continue
        ngpu = int(read(os.path.join(d, 'ngpu'), '1'))
        if a.max_ngpu and ngpu > a.max_ngpu:
            continue
        inp = os.path.join(d, 'calliope.in')
        cost = 1
        for k in ('nlx', 'nly', 'nlz', 'nstep'):
            cost *= nml_int(inp, k, 1)
        rows.append((-cost, r, ngpu, (st or '-') + (' (out3d conversion pending)' if pending else ''), folder))
    rows.sort()
    for cost, r, ngpu, st, folder in rows:
        print('%s  %d  %.3e  %s  %s' % (r, ngpu, -cost, folder, st) if a.long else r)
    return 0


# ------------------------------------------------------------------------------------------------------- stop
def with_t_end(text, t):
    """calliope.in text with t_end = t in &time_parameters (an existing t_end replaced)."""
    m = re.search(r'^\s*&time_parameters\b.*?^\s*/', text, re.M | re.S | re.I)
    if not m:
        sys.exit('no &time_parameters')
    group = m.group(0)
    line = '  t_end = %rd0' % t
    if re.search(r'^\s*t_end\s*=', group, re.M | re.I):
        new = re.sub(r'^\s*t_end\s*=.*$', line, group, flags=re.M | re.I)
    else:
        new = re.sub(r'^(\s*/)\Z', line + r'\n\1', group, flags=re.M)
    return text[:m.start()] + new + text[m.end():]


def cmd_stop(a):
    R = os.path.abspath(a.runroot)
    bad = 0
    for run, (t, why) in CALLIOPE_STOP.items():
        d = os.path.join(R, run)
        inp, tf = os.path.join(d, 'calliope.in'), os.path.join(d, 't_end')
        if not os.path.exists(inp):
            print('%s: no calliope.in' % run)
            bad = 1
            continue
        if os.path.exists(inp + '.nominal'):
            ok = nml_value(inp, 't_end') == repr(t) + 'e0' and read(tf) == repr(t)
            print('%s: %s' % (run, 'already stopped at t = %r' % t if ok else 'CONFLICT: *.nominal exists, '
                                   'but t_end differs from %r' % t))
            bad |= not ok
            continue
        if not never_ran(d):
            print('%s: has run, not changed' % run)
            bad = 1
            continue
        text = read(inp) + '\n'
        shutil.copyfile(inp, inp + '.nominal')
        if os.path.exists(tf):
            shutil.copyfile(tf, tf + '.nominal')
        with open(inp, 'w') as f:
            f.write(with_t_end(text, t))
        with open(tf, 'w') as f:
            f.write('%r\n' % t)
        print('%s: t_end %s -> %r (%s)' % (run, read(tf + '.nominal', '-'), t, why))
    return bad


# ---------------------------------------------------------------------------------------------------- folders
ROLLBACK = '''#!/bin/sh
# Undo `make_cases.py folders --apply` of this directory (dry run unless DRY=0). Paths are relative to T.
# moves.txt "src dst rel", last line first: unlink the link T/src if it reads rel, move T/dst back to T/src.
# links.txt "path target": unlink T/path if it reads target. Then remove the empty directories of the folders.
T=%(T)s
A=$(cd "$(dirname "$0")" && pwd)
cd "$T" || exit 1
run() { if [ "${DRY:-1}" = 0 ]; then "$@" || exit 1; else echo "dry: $*"; fi; }
while read -r p t; do
  if [ -L "$p" ] && [ "$(readlink "$p")" = "$t" ]; then run unlink "$p"; else echo "SKIP link $p"; fi
done < "$A/links.txt"
reverse() { awk '{l[NR] = $0} END {for (i = NR; i > 0; i--) print l[i]}' "$1"; }
reverse "$A/moves.txt" | while read -r src dst rel; do
  if [ -L "$src" ]; then
    [ "$(readlink "$src")" = "$rel" ] || { echo "SKIP $src: link reads $(readlink "$src")"; continue; }
    run unlink "$src"
  elif [ -e "$src" ]; then
    echo "SKIP $src: exists and is no link"; continue
  fi
  run mv "$dst" "$src"
done || exit 1
for F in %(FOLDERS)s; do [ -d "$F" ] && run find "$F" -depth -type d -empty -delete; done
echo "rollback: $( [ "${DRY:-1}" = 0 ] && echo done || echo 'dry run (DRY=0 to apply)')"
'''


def cmd_folders(a):
    """Move the cases of the run root R (= T/runs) to T/<Folder>/<run> and leave relative links in R, as on ps1."""
    R = os.path.abspath(a.runroot)
    T = os.path.dirname(R)
    if os.path.basename(R) != 'runs' or os.path.islink(R):
        sys.exit('folders: the run root must be a directory named runs (the folders go next to it): %s' % R)
    moves, linked, bad = [], 0, []
    for root, dirs, files in os.walk(R):
        dirs.sort()
        if 'calliope.in' in files:
            dirs[:] = []
            run = os.path.relpath(root, R)
            tpath = read(os.path.join(root, 'calliope_case'))
            if tpath.split('/')[0] not in FOLDERS or tpath.split('/', 1)[-1] != run:
                bad.append('%s: calliope_case %r is not <Folder>/%s' % (run, tpath, run))
            elif os.path.lexists(os.path.join(T, tpath)):
                bad.append('%s: %s exists' % (run, tpath))
            else:
                moves.append(('runs/' + run, tpath, os.path.relpath(os.path.join(T, tpath), os.path.dirname(root))))
        else:
            dirs[:] = [d for d in dirs if not d.startswith('chunk_')]
            linked += sum(os.path.islink(os.path.join(root, d)) and os.path.isfile(os.path.join(root, d, 'calliope.in'))
                          for d in dirs)
    links, same = [], 0
    C = tests_dir(a.calliope)
    for F in FOLDERS:
        for root, dirs, files in os.walk(os.path.join(C, F)):
            for n in sorted(dirs + files):
                if not os.path.islink(os.path.join(root, n)):
                    continue
                p, t = os.path.relpath(os.path.join(root, n), C), os.readlink(os.path.join(root, n))
                clash = [d for _, d, _ in moves if p == d or p.startswith(d + '/') or d.startswith(p + '/')]
                if clash:
                    bad.append('%s overlaps the case moved to %s' % (p, clash[0]))
                elif not os.path.lexists(os.path.join(T, p)):
                    links.append((p, t))
                elif os.path.islink(os.path.join(T, p)) and os.readlink(os.path.join(T, p)) == t:
                    same += 1
                else:
                    bad.append('%s exists and is not a link to %s' % (p, t))
    count = collections.Counter(d.split('/')[0] for _, d, _ in moves)
    print('folders: %s: %d cases to move (%s), %d already moved; %d summary links to make, %d already there'
          % (T, len(moves), ', '.join('%s %d' % (F, count[F]) for F in FOLDERS), linked, len(links), same))
    for m in moves[:3]:
        print('  %s -> %s (link %s)' % m)
    for b in bad:
        print('CONFLICT: %s' % b)
    if bad or not a.apply or not (moves or links):
        if not bad and not a.apply and (moves or links):
            print('dry run; --apply to move')
        return 1 if bad else 0
    A = os.path.join(T, 'attic', 'folders_' + time.strftime('%Y%m%d-%H%M%S'))
    os.makedirs(A)
    with open(os.path.join(A, 'rollback.sh'), 'w') as f:
        f.write(ROLLBACK % {'T': shlex.quote(T), 'FOLDERS': ' '.join(FOLDERS)})
    os.chmod(os.path.join(A, 'rollback.sh'), 0o755)
    with open(os.path.join(A, 'moves.txt'), 'w') as log:
        for src, dst, rel in moves:
            os.makedirs(os.path.dirname(os.path.join(T, dst)), exist_ok=True)
            os.rename(os.path.join(T, src), os.path.join(T, dst))   # EXDEV across file systems: no copy
            log.write('%s %s %s\n' % (src, dst, rel))
            log.flush()
            os.symlink(rel, os.path.join(T, src))
    with open(os.path.join(A, 'links.txt'), 'w') as log:
        for p, t in links:
            os.makedirs(os.path.dirname(os.path.join(T, p)), exist_ok=True)
            os.symlink(t, os.path.join(T, p))
            log.write('%s %s\n' % (p, t))
    broken = [p for p, _ in links if not os.path.exists(os.path.join(T, p))]
    for p in broken:
        print('BROKEN summary link: %s' % p)
    print('folders: moved %d cases, made %d summary links; record and rollback.sh in %s' % (len(moves), len(links), A))
    return 1 if broken else 0


# --------------------------------------------------------------------------------------------------- ngpu-map
def cmd_ngpu_map(a):
    cases = calliope_cases(a.calliope)
    rule = ngpu_rule(cases)
    R = os.path.abspath(a.runroot) if a.runroot else None
    final = {}
    for c in cases:
        n, np_, why = rule[c.run]
        f = read(os.path.join(R, c.run, 'ngpu')) if R else ''
        if f and int(f) != n:
            n, why = int(f), 'ngpu file of the run directory (rule: %d, %s)' % (n, why)
        final[c.run] = (n, np_, why)
    count = collections.Counter(n for n, _, _ in final.values())
    out = ['# ngpu of the MHD_COMP_ISOTH cases on TSUBAME4.0', '',
           'Written by `make_cases.py ngpu-map` (rule `ngpu_rule`); a file `ngpu` in a run directory overrides it.',
           'A case runs on 1 GPU unless it is a member of a decomposition test (`DECOMP`): the reference of a family',
           'runs on 1 GPU and the other members on 2, 4, ... GPUs in the order of their Calliope np (the j-th',
           'distinct np gets 2^j), capped at %d and at the largest power of two that divides nlx and nly.' % MAX_GPU,
           '', 'Cases: %d; ngpu %s.' % (len(final), ', '.join('%d: %d cases' % kv for kv in sorted(count.items()))),
           '']
    for F in FOLDERS:
        out += ['## %s' % F, '', '| case | Calliope np | ngpu | reason |', '|---|---:|---:|---|']
        for c in cases:
            if c.folder == F:
                n, np_, why = final[c.run]
                out.append('| %s | %d | %d | %s |' % (c.run, np_, n, why))
        out.append('')
    text = '\n'.join(out)
    if a.output:
        with open(a.output, 'w') as f:
            f.write(text)
    else:
        print(text)
    print('ngpu: %s' % ', '.join('%d -> %d cases' % kv for kv in sorted(count.items())), file=sys.stderr)
    return 0


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = p.add_subparsers(dest='cmd')
    s = sp.add_parser('setup')
    s.add_argument('--calliope', required=True)
    s.add_argument('--runroot', required=True)
    s.add_argument('--only', action='append', metavar='GLOB')
    s = sp.add_parser('ic')
    s.add_argument('--calliope', required=True)
    s.add_argument('--runroot', required=True)
    s.add_argument('--python', default=sys.executable)
    s.add_argument('--batch', type=int, default=20, help='cases per call of regen_inputs.py')
    s.add_argument('--only', action='append', metavar='GLOB')
    s = sp.add_parser('check')
    s.add_argument('--gallope', default=os.path.normpath(os.path.join(HERE, '..', '..')))
    s.add_argument('--runroot', required=True)
    s.add_argument('--only', action='append', metavar='GLOB')
    s = sp.add_parser('list')
    s.add_argument('--runroot', required=True)
    s.add_argument('--folder', action='append', choices=FOLDERS)
    s.add_argument('--status', action='append')
    s.add_argument('--max-ngpu', type=int, default=0)
    s.add_argument('--long', action='store_true', help='also print ngpu, cost, folder and exit_status')
    s = sp.add_parser('stop')
    s.add_argument('--runroot', required=True)
    s = sp.add_parser('folders')
    s.add_argument('--calliope', required=True)
    s.add_argument('--runroot', required=True)
    s.add_argument('--apply', action='store_true', help='move (default: dry run)')
    s = sp.add_parser('ngpu-map')
    s.add_argument('--calliope', required=True)
    s.add_argument('--runroot')
    s.add_argument('-o', '--output')
    a = p.parse_args()
    if a.cmd == 'list' and not a.status:
        a.status = ['all']
    cmds = {'setup': cmd_setup, 'ic': cmd_ic, 'check': cmd_check, 'list': cmd_list, 'stop': cmd_stop,
            'folders': cmd_folders, 'ngpu-map': cmd_ngpu_map}
    if a.cmd not in cmds:
        p.print_help()
        return 2
    return cmds[a.cmd](a)


if __name__ == '__main__':
    sys.exit(main())

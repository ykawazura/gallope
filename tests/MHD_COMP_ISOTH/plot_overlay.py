#!/usr/bin/env python3
"""
Calliope's paper figures of the MHD_COMP_ISOTH tests with the Gallope results drawn over them.

  python3 plot_overlay.py list  FOLDER
  python3 plot_overlay.py fetch FOLDER                  inputs of the figures from TSUBAME into --gallope-runs
  python3 plot_overlay.py plot  FOLDER [name ...]       -> fig/<FOLDER>/<name>.pdf (default: all figures)

FOLDER is Euler1D, MHD1D, Smooth_wave, OT2D or Turb3D. The figures are those of
<calliope>/<FOLDER>/plot_<folder>.py (--calliope: tests/MHD_COMP_ISOTH of a Calliope worktree at 1d924a9), drawn
twice on the same axes: from --calliope-runs (Calliope's runs/), then from --gallope-runs. The second pass is
restyled: solid lines dotted, dashed lines dash-dotted, markers small and white inside. Its curves identical to
those of the first pass (exact solutions, references, guides, legend proxies), its texts and its legends are
removed, and a key is added below the axes. The colour maps (figures drawn with figstyle.maps: ot_rho, ot_J,
turb_slice) are written as <name>_calliope.pdf, <name>_gallope.pdf and <name>_diff.pdf, the last one with
(Gallope - Calliope)/max|Calliope| of every map on one colour scale.

fetch first draws the figures from --calliope-runs and records the files they open there (an audit hook). It then
copies the same files from RUNROOT on TSUBAME (--remote is CALLIOPE_ROOT there, whose tests/MHD_COMP_ISOTH/runs
is RUNROOT) with figdata.fetch_copy, and runs the Extract items of the plot script (arrays of large outputs) there
with the Python of the venv (--remote-python). Nothing is written on TSUBAME and no local file is overwritten.
An out1d record that RUNROOT lacks under Calliope's name is copied from the Gallope record of the same t
(fetch_by_time).
"""
import argparse
import importlib.util
import os
import re
import shlex
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = {'Euler1D': 'plot_euler1d.py', 'MHD1D': 'plot_mhd1d.py', 'Smooth_wave': 'plot_smooth_wave.py',
           'OT2D': 'plot_ot2d.py', 'Turb3D': 'plot_turb3d.py'}
CALLIOPE = os.path.expanduser('~/Dropbox/Codes/calliope_dev/tests/MHD_COMP_ISOTH')
KEY_H = 0.8                                  # inches added below the axes for the key
# Arguments of the Gallope pass. turb1_N128_a and turb1_N256 stop 0.00003 and 0.00006 tau0 before 3 tau0
# (tt >= t_end - dt/2, t_end = 3 pi); analyze_stage4 takes a record at 3 tau0 only if 3 tau0 <= t_last, so the
# Gallope summary has no '3.00' spectrum of N = 128 and 256 (Calliope's runs, older executables, stopped past 3 tau0).
GALLOPE_KW = {('Turb3D', 'turb_spec'): dict(keys=('1.00',))}
OPENED = None                                # files opened while recording (fetch), else None


def audit(event, args):
    if OPENED is not None and event == 'open' and isinstance(args[0], (str, bytes, os.PathLike)):
        OPENED.add(os.path.realpath(os.fsdecode(args[0])))


def load(calliope, folder):
    """The plot script as a module (its __main__ block does not run); fs and fd are its figstyle and figdata."""
    sys.dont_write_bytecode = True
    path = os.path.join(calliope, folder, SCRIPTS[folder])
    spec = importlib.util.spec_from_file_location(os.path.splitext(SCRIPTS[folder])[0], path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class Capture:
    """figstyle.panels/maps/save replaced: save keeps the figure, panels hands out the axes of the first pass again."""

    def __init__(self, fs):
        self.fs, self.orig = fs, {k: getattr(fs, k) for k in ('panels', 'maps', 'save')}
        fs.panels, fs.maps, fs.save = self.panels, self.maps, self.save
        self.reuse = self.kind = self.fig = self.maps_args = None

    def panels(self, *a, **k):
        self.kind = 'panels'
        return self.reuse if self.reuse is not None else self.orig['panels'](*a, **k)

    def maps(self, *a, **k):
        self.kind, self.maps_args = 'maps', (a, k)
        return self.orig['maps'](*a, **k)

    def save(self, fig, test_dir, name):
        self.fig = fig

    def draw(self, f, runs, reuse=None):
        self.reuse, self.kind, self.fig = reuse, None, None
        f(runs)
        if self.fig is None:
            raise RuntimeError('the figure was not saved')
        return self.kind, self.fig


def state(fig):
    ax = {a: dict(lines=list(a.lines), coll=list(a.collections), patches=list(a.patches), texts=list(a.texts),
                  artists=list(a.artists), legend=a.get_legend()) for a in fig.axes}
    return ax, dict(texts=list(fig.texts), legends=list(fig.legends))


def same_line(l, o):
    from matplotlib.colors import to_rgba
    return (np.array_equal(l.get_xydata(), o.get_xydata(), equal_nan=True) and l.get_linestyle() == o.get_linestyle()
            and l.get_marker() == o.get_marker() and to_rgba(l.get_color()) == to_rgba(o.get_color()))


def same_shape(p, o):
    a, b = p.get_path().vertices, o.get_path().vertices
    return type(p) is type(o) and a.shape == b.shape and np.allclose(p.get_transform().transform(a),
                                                                      o.get_transform().transform(b))


def same_coll(c, o):
    pc, po = c.get_paths(), o.get_paths()
    return (type(c) is type(o) and len(pc) == len(po)
            and all(np.array_equal(a.vertices, b.vertices) for a, b in zip(pc, po))
            and np.array_equal(c.get_offsets(), o.get_offsets()))


def restyle(l):
    from matplotlib.colors import to_rgba
    ls = l.get_linestyle()
    l.set_linestyle({'-': ':', '--': '-.', '-.': (0, (3, 1, 1, 1, 1, 1)), ':': (0, (1, 3))}.get(ls, ls))
    if l.get_marker() not in (None, 'None', '', ' '):
        edge = l.get_markeredgecolor()
        l.set_markeredgecolor(l.get_color() if to_rgba(edge)[3] == 0 else edge)
        l.set_markerfacecolor('white')
        l.set_markersize(0.45 * l.get_markersize())
        l.set_markeredgewidth(0.6 * l.get_markeredgewidth())
    l.set_zorder(max(l.get_zorder(), 4))


def merge(fig, before):
    """Restyle what the second pass added to fig and remove its duplicates, texts and legends; the number of
    restyled and of removed curves."""
    ax0, f0 = before
    kept = dropped = 0
    for a, s in ax0.items():
        for l in [l for l in a.lines if l not in s['lines']]:
            if any(same_line(l, o) for o in s['lines']):
                l.remove()
                dropped += 1
            else:
                restyle(l)
                kept += 1
        for c in [c for c in a.collections if c not in s['coll']]:
            if any(same_coll(c, o) for o in s['coll']):
                c.remove()
            else:
                c.set_facecolor('none')
                c.set_linestyle(':')
        for p in [p for p in a.patches if p not in s['patches']]:
            if any(same_shape(p, o) for o in s['patches']):
                p.remove()
        for t in [t for t in a.texts if t not in s['texts']]:
            t.remove()
        for x in [x for x in a.artists if x not in s['artists']]:
            x.remove()
        a.legend_ = s['legend']
    for t in [t for t in fig.texts if t not in f0['texts']]:
        t.remove()
    fig.legends[:] = f0['legends']
    return kept, dropped


def add_key(fs, fig):
    """Room of KEY_H inches below the axes (the figure grows, the axes keep their size) and the key in it."""
    import matplotlib.pyplot as plt
    W, H = fig.get_size_inches()
    s, d = H / (H + KEY_H), KEY_H / (H + KEY_H)
    for a in fig.axes:
        p = a.get_position(original=True)
        a.set_position([p.x0, p.y0 * s + d, p.width, p.height * s])
    for t in fig.texts:
        if t.get_transform() == fig.transFigure:
            x, y = t.get_position()
            t.set_position((x, y * s + d))
    fig.set_size_inches(W, H + KEY_H)
    h = [plt.Line2D([], [], color='k', lw=fs.LW, ls='-', marker='o', ms=fs.MS * 0.7, mew=fs.MEW, label='Calliope'),
         plt.Line2D([], [], color='k', lw=fs.LW, ls=':', marker='o', ms=fs.MS * 0.7 * 0.45, mew=fs.MEW * 0.6,
                    mfc='white', label='Gallope')]
    fig.legend(handles=h, loc='lower center', ncol=2, frameon=False, fontsize=fs.LEGEND_FONTSIZE,
               bbox_to_anchor=(0.5, 0.0), handlelength=2.4, columnspacing=2.0)


def save(fig, out, name):
    import matplotlib.pyplot as plt
    os.makedirs(out, exist_ok=True)
    path = os.path.join(out, name + '.pdf')
    fig.savefig(path)
    plt.close(fig)
    print('wrote', os.path.relpath(path))


def diff_maps(fs, cap, fc, fg, out, name):
    """(Gallope - Calliope)/max|Calliope| of every map (axes with an image), on the axes of figstyle.maps."""
    pc = [a for a in fc.axes if a.images]
    pg = [a for a in fg.axes if a.images]
    a, k = cap.maps_args
    fig, ax = cap.orig['maps'](*a, **k)
    vmax, ims = 0.0, []
    for x, ac, ag in zip(ax.flat, pc, pg):
        c = np.ma.filled(ac.images[0].get_array().astype(np.float64), 0.0)
        g = np.ma.filled(ag.images[0].get_array().astype(np.float64), 0.0)
        r = (g - c) / np.abs(c).max()
        vmax = max(vmax, np.abs(r).max())
        print('%s %s: max|G - C|/max|C| = %.2e' % (name, ac.get_title(), np.abs(r).max()))
        ims.append(fs.image(x, r.T, cmap='RdBu_r'))
        x.set_title(ac.get_title(), fontsize=fs.LEGEND_FONTSIZE + 2)
    for im in ims:
        im.set_clim(-vmax, vmax)
    fs.colorbar(fig, ax[:, -1], ims[0], r'$\Delta$')
    save(fig, out, name + '_diff')


def plot(mod, opt):
    import matplotlib.pyplot as plt
    cap = Capture(mod.fs)
    out = os.path.join(opt.out, opt.folder)
    for name in opt.names or list(mod.FIGS):
        f = mod.FIGS[name]
        print('== %s (Calliope)' % name)
        kind, fc = cap.draw(f, opt.calliope_runs)
        before = state(fc)
        print('== %s (Gallope)' % name)
        try:
            if kind == 'maps':
                _, fg = cap.draw(f, opt.gallope_runs)
            else:
                for a in fc.axes:
                    a.set_prop_cycle(None)       # the default colours start again
                kw = GALLOPE_KW.get((opt.folder, name))
                g = (lambda runs: f(runs, **kw)) if kw else f
                _, fg = cap.draw(g, opt.gallope_runs, reuse=(fc, panel_axes(fc)))
        except FileNotFoundError as e:
            print('skipped %s: missing %s' % (name, e.filename or e))
            plt.close('all')
            continue
        except (TypeError, ValueError, KeyError) as e:
            # e.g. cp_diag formats gamma of runs that grew in Calliope but have no measurable growth in Gallope
            print('skipped %s: the plot script fails on the Gallope data (%s: %s)' % (name, type(e).__name__, e))
            plt.close('all')
            continue
        if kind == 'maps':
            diff_maps(mod.fs, cap, fc, fg, out, name)
            save(fc, out, name + '_calliope')
            save(fg, out, name + '_gallope')
        else:
            kept, dropped = merge(fc, before)
            print('%s: %d Gallope curves, %d identical to Calliope removed' % (name, kept, dropped))
            add_key(mod.fs, fc)
            save(fc, out, name)


def panel_axes(fig):
    """The ax array that figstyle.panels returned for fig: its axes at their places on the grid."""
    gs = fig.axes[0].get_subplotspec().get_gridspec()
    nrow, ncol = gs.get_geometry()
    ax = np.empty((nrow, ncol), dtype=object)
    for a in fig.axes:
        ss = a.get_subplotspec()
        if ss is not None and ss.get_gridspec() is gs:
            ax[ss.rowspan.start, ss.colspan.start] = a
    return ax


OUT1D = re.compile(r'(?P<case>.+)/out1d/line_\d{5}\.dat$')


def record_time(head):
    m = re.match(r'# t = *(\S+),', head)
    return float(m.group(1)) if m else None


def fetch_by_time(opt, fd, o, paths):
    """Copy, for every out1d/line_NNNNN.dat of paths, the Gallope record of the same t to that path.

    The record numbers of the two codes need not agree. nline1d starts at 0 in every run, so a run that the job
    budget split into chunk_* numbers its records from 0 again. write_due compares floor((tt + dt)/intvl) with
    floor(tt/intvl), and Calliope's build (Intel, fp-model fast) may divide by intvl as a product with its
    reciprocal, which adds a record one step before a multiple of intvl when tt + dt falls just below it.
    The record is looked for in the case directory, then in its chunk_*; its t must equal that of Calliope's
    file to 1e-12 relative."""
    rdir = '%s/%s' % (fd.remote_dir(o), fd.REMOTE_RUNS)
    for p in paths:
        m = OUT1D.match(p)
        if not m:
            continue
        with open(os.path.join(opt.calliope_runs, p)) as f:
            t = record_time(f.readline())
        c = shlex.quote(m.group('case'))
        cmd = 'cd %s && grep -m1 -H "^# t" -- %s/out1d/line_*.dat %s/chunk_*/out1d/line_*.dat 2>/dev/null' % (
            rdir, c, c)
        hits = []
        for ln in subprocess.run(fd.ssh_cmd(o, cmd), capture_output=True, text=True).stdout.splitlines():
            src, _, head = ln.partition(':')
            tg = record_time(head)
            if t is not None and tg is not None and abs(tg - t) <= 1e-12*max(1.0, abs(t)):
                hits.append(src)
        if not hits:
            print('%s: no Gallope record at t = %r' % (p, t))
            continue
        src = min(hits, key=lambda s: ('/chunk_' in s, s))
        r = subprocess.run(fd.ssh_cmd(o, 'cd %s && cat -- %s' % (rdir, shlex.quote(src))), capture_output=True)
        if r.returncode != 0 or not r.stdout:
            print('%s: copy of %s failed' % (p, src))
            continue
        dst = os.path.join(opt.gallope_runs, p)
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        with open(dst + '.part', 'wb') as f:
            f.write(r.stdout)
        os.replace(dst + '.part', dst)
        print('copied %s from %s (t = %r)' % (dst, src, t))


def fetch(mod, opt):
    global OPENED
    import matplotlib.pyplot as plt
    fd = mod.fd
    Capture(mod.fs)
    runs = os.path.realpath(opt.calliope_runs) + os.sep
    sys.addaudithook(audit)
    OPENED = set()
    for name, f in mod.FIGS.items():
        try:
            f(opt.calliope_runs)
        except FileNotFoundError as e:
            print('%s: missing in %s: %s' % (name, opt.calliope_runs, e.filename or e))
        plt.close('all')
    paths, OPENED = sorted(p[len(runs):] for p in OPENED if p.startswith(runs)), None
    extract = [i for i in getattr(mod, 'FETCH', []) if isinstance(i, fd.Extract)]
    skip = {t for i in extract for t in i.targets()}
    copy = [p for p in paths if p not in skip]
    print('%d files to copy, %d to extract' % (len(copy), len(extract)))
    orig = fd.ssh_cmd
    fd.ssh_cmd = lambda o, cmd: orig(o, cmd.replace('&& python3 -B -', '&& %s -B -' % opt.remote_python))
    o = argparse.Namespace(runs=opt.gallope_runs, host=opt.host, remote=opt.remote)
    try:
        fd.fetch_copy(o, fd.Copy(*copy))
    except SystemExit as e:
        print(e)
        fetch_by_time(opt, fd, o, [p for p in copy if not os.path.exists(os.path.join(opt.gallope_runs, p))])
    for i in extract:
        fd.fetch_extract(o, i)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('cmd', choices=('list', 'fetch', 'plot'))
    p.add_argument('folder', choices=sorted(SCRIPTS))
    p.add_argument('names', nargs='*')
    p.add_argument('--calliope', default=CALLIOPE)
    p.add_argument('--calliope-runs')
    p.add_argument('--gallope-runs', default=os.path.join(HERE, 'runs_gallope'))
    p.add_argument('--out', default=os.path.join(HERE, 'fig'))
    p.add_argument('--host', default='tsubame')
    p.add_argument('--remote', default='~/work/calliope_1d924a9')
    p.add_argument('--remote-python', default='~/work/venv-isoth/bin/python')
    opt = p.parse_args()
    opt.calliope_runs = opt.calliope_runs or os.path.join(opt.calliope, 'runs')
    mod = load(opt.calliope, opt.folder)
    unknown = [n for n in opt.names if n not in mod.FIGS]
    if unknown:
        p.error('unknown figures: %s' % ', '.join(unknown))
    if opt.cmd == 'list':
        for name, f in mod.FIGS.items():
            print('%-18s %s' % (name, ' '.join((f.__doc__ or '').split())[:100]))
    elif opt.cmd == 'fetch':
        fetch(mod, opt)
    else:
        plot(mod, opt)


if __name__ == '__main__':
    main()

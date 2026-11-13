# -*- coding: utf-8 -*-
'''
Archive and merge the output of a run that was restarted one or more times.

Gallope reopens every output file from scratch when it restarts -- there is no
append branch in io.F90 -- so job N+1 destroys job N's time history unless it is
archived first.  A run split across several jobs is therefore kept as a numbered
sequence: what the Nth job wrote becomes gallope.out.ncN, cfl.datN, out2dN/, and
so on, beside the live files the next job is about to overwrite.

  [$] python misc/backup_output.py
      archives the live output under the next free number.  The job scripts call
      this themselves once the solver has stopped writing, so every job that
      ends normally -- the last one included -- is already archived; run it by
      hand only after a job that was killed outright.

  [$] python misc/backup_output.py 3
      the same, with the number given explicitly.

  [$] python misc/backup_output.py move
      move instead of copy, which is what production sizes want.  restart/ is
      copied whatever this says: the run may still be continued from it, and the
      archived copy is what lets an earlier job be replayed on its own.

  [$] python misc/backup_output.py all
      merges gallope.out.nc1, gallope.out.nc2, ... into gallope.out.nc_all, and
      out2d1/, out2d2/ ... into out2d_all/.  Point the diagnostics at the whole
      history with

          restart_num = '_all'               # load.py

      Each backup is cut at the first time the next one also carries, so a job
      that wrote diagnostics past the checkpoint it was restarted from leaves no
      duplicate or non-monotonic records.  restart/ is a state rather than a
      history and is not merged -- the highest-numbered restartN is already the
      checkpoint the run ended on.

Merging netCDF needs either the netCDF4 module or NCO (ncrcat/ncks); everything
else -- out2d/, out3d/, cfl.dat, rms.dat -- merges with numpy alone, so the text
and binary halves still work on machines that have neither.

There is no MODEL argument any more.  Every model writes some subset of the same
files, so the list below is their union and whatever a run did not write is
simply skipped -- which is also why a model that never had a branch here
(RMHD, EMHD, MHD_COMP_*) no longer fails.  A MODEL name is still accepted in the
leading position, and ignored, so that existing job scripts keep working.
'''
import warnings
warnings.filterwarnings('ignore')

import sys
import os.path
import os
import re
import shutil
import glob
import subprocess

# numpy is imported by the merging branch only.  Archiving is plain file
# copying, and job scripts call it from compute nodes that may carry nothing but
# the standard library.

argv = sys.argv

# 'move' is a flag rather than a position, so that it can follow the number.
args = [a for a in argv[1:] if a != 'move']
move = len(args) != len(argv) - 1

# The leading argument used to be a MODEL name choosing which files to archive.
# It chooses nothing now, but is still swallowed so that job scripts and tooling
# written against the old signature keep working.  Recognised by shape -- upper
# case, as every MODEL is -- rather than by a list, so that a model this script
# has never heard of is accepted too.
if args and re.match(r'^[A-Z][A-Z0-9_]*$', args[0]):
  args = args[1:]

restart_num = args[0] if args else None

if restart_num is not None and restart_num != 'all' and not restart_num.isdigit():
  raise SystemExit("restart_num must be a number or 'all', not %r -- and the "
                   "MODEL it used to take\nis gone, so a lower-case one is read "
                   "as the number." % restart_num)


# Every model in both codes writes some subset of these -- checked against the
# runname-derived names in src/model/*/io.F90 and src/diagnostics_common.F90 --
# so one list serves them all and a file the run did not write is skipped.
filenames = [
              # Some job scripts redirect stdout to out.std rather than
              # gallope.out.std; whichever is absent is skipped.
              'gallope.out.std',
              'out.std',
              'gallope.in',
              'cfl.dat',
              'rms.dat',
              'gallope.out.nc',
              'gallope.modes.out',
              'gallope.out.SF2.nc',
              'gallope.out.kpar.nc',
              'gallope.out.nltrans.nc',
            ]
dirnames  = ['out2d', 'out3d', 'restart']



#==============================================================================
# helpers
#==============================================================================
ALL   = '_all'
EPS   = 1e-9
CHUNK = 1 << 24


def numbers(name):
  '''The numbers under which `name` has been archived, in numeric order.

  Matched by regexp rather than by the glob alone, so that out2d_all/ and any
  other decorated sibling cannot be mistaken for a backup.'''
  pat = re.compile(re.escape(name) + r'(\d+)$')
  return sorted(int(pat.match(p).group(1))
                for p in glob.glob(name + '[0-9]*') if pat.match(p))


def backups(name):
  return [name + str(n) for n in numbers(name)]


def next_number():
  archived = [n for name in filenames + dirnames for n in numbers(name)]
  return str(max(archived) + 1 if archived else 1)


def first_tt_txt(path):
  '''First value in column 1, skipping blank and non-numeric header lines.'''
  with open(path) as f:
    for line in f:
      fields = line.split(None, 1)
      if fields:
        try:
          return float(fields[0])
        except ValueError:
          continue
  return None


def read_tt_nc(path):
  from scipy.io import netcdf_file
  f = netcdf_file(path, 'r', mmap=False)
  try:
    return np.array(f.variables['tt'][:], dtype=float)
  finally:
    f.close()


def cut_indices(tts):
  '''Records of backup i to keep: everything before backup i+1's first time.'''
  return [len(tt) if i == len(tts) - 1
          else int(np.searchsorted(tt, tts[i + 1][0] - EPS))
          for i, tt in enumerate(tts)]


def copy_n_bytes(src, dst, n):
  while n > 0:
    buf = src.read(min(CHUNK, n))
    if not buf:
      break
    dst.write(buf)
    n -= len(buf)


#==============================================================================
# merging one family of numbered backups at a time
#==============================================================================
def merge_nc(srcs, dst):
  '''Concatenate along the tt (unlimited) dimension with the overlap trimmed.

  Uses netCDF4 when it is installed and NCO otherwise.  Writing netCDF3 with
  scipy is not an option here: its writer insists the unlimited dimension come
  first, which reorders the header into something the netCDF C library refuses
  to open.'''
  tts = [read_tt_nc(p) for p in srcs]
  # A diagnostic that was never switched on leaves a header and no records.
  # Such a backup holds no history, and trimming the one before it against a
  # time it does not have is meaningless, so it is dropped.
  live = [i for i, tt in enumerate(tts) if tt.size]
  if not live:
    # Nothing to concatenate -- but the merged view still has to carry the
    # stream.  load.py picks a history by suffix, and a reader set to '_all'
    # has to find the same set of files it finds live: the header alone is
    # what defines kpbin_log and friends, and the analysis scripts use them
    # whether or not any record was ever written.  Pass the first one through.
    shutil.copy2(srcs[0], dst)
    print ('%-28s <- header only, no records' % dst)
    return None
  srcs = [srcs[i] for i in live]
  tts  = [tts [i] for i in live]
  keep = cut_indices(tts)

  try:
    from netCDF4 import Dataset
  except ImportError:
    return merge_nc_nco(srcs, dst, keep)

  parts = [Dataset(p) for p in srcs]
  ref   = parts[0]
  out   = Dataset(dst, 'w', format='NETCDF3_CLASSIC')
  for a in ref.ncattrs():
    out.setncattr(a, ref.getncattr(a))
  for n, dim in ref.dimensions.items():
    out.createDimension(n, None if dim.isunlimited() else len(dim))
  for n, v in ref.variables.items():
    o = out.createVariable(n, v.dtype, v.dimensions)
    for a in v.ncattrs():
      o.setncattr(a, v.getncattr(a))
    if 'tt' in v.dimensions:
      o[:] = np.concatenate([np.asarray(part.variables[n][0:k])
                             for part, k in zip(parts, keep)], axis=0)
    else:
      o[...] = ref.variables[n][...]
  out.close()
  for part in parts:
    part.close()
  return keep


def merge_nc_nco(srcs, dst, keep):
  '''Fallback for machines without netCDF4.  NCO is only tried, never probed:
  an ncrcat that is on PATH but cannot load its libraries aborts rather than
  reporting an error, so the failure has to be caught around the real call.'''
  parts = []
  try:
    for i, (src, k) in enumerate(zip(srcs, keep)):
      if k == len(read_tt_nc(src)):
        parts.append(src)
      else:
        tmp = dst + '.part%d.tmp' % i
        subprocess.run(['ncks', '-O', '-h', '-d', 'tt,0,%d' % (k - 1), src, tmp],
                       check=True)
        parts.append(tmp)
    subprocess.run(['ncrcat', '-O', '-h'] + parts + [dst], check=True)
  except (OSError, subprocess.CalledProcessError) as err:
    print ('  SKIPPED %s -- no netCDF4 module, and NCO failed: %s\n'
           '          Merge the netCDF files where netCDF4 or a working NCO is\n'
           '          installed; everything else is already complete.'
           % (dst, err))
    keep = None
  for p in parts:
    if p.endswith('.tmp') and os.path.exists(p):
      os.remove(p)
  return keep


def merge_table(srcs, dst):
  '''cfl.dat / rms.dat: whitespace tables whose first column is tt.'''
  bounds = [first_tt_txt(s) for s in srcs[1:]] + [None]
  kept = []
  with open(dst, 'w') as o:
    for src, tmax in zip(srcs, bounds):
      n = 0
      with open(src) as f:
        for line in f:
          if tmax is not None:
            fields = line.split(None, 1)
            if fields:
              try:
                if float(fields[0]) >= tmax - EPS:
                  break
              except ValueError:
                pass
          o.write(line)
          n += 1
      kept.append(n)
  return kept


def merge_text(srcs, dst):
  '''Logs with no time column: plain concatenation.'''
  with open(dst, 'wb') as o:
    for src in srcs:
      with open(src, 'rb') as f:
        shutil.copyfileobj(f, o, CHUNK)
  return None


def merge_bindir(srcs, dst):
  '''out2d / out3d: fixed-size records indexed by the directory's time.dat.

  Only *.dat is merged.  MPI-IO leaves an 8-byte <name>.dat-<id>.lock beside
  every stream it writes; those hold no records, and a merge that took them for
  data would both mis-count and produce garbage.'''
  times = [os.path.join(s, 'time.dat') for s in srcs]
  missing = [t for t in times if not os.path.exists(t)]
  if missing:
    print ('  SKIPPED ' + dst + ' -- no ' + missing[0])
    return None

  lines = [[l for l in open(t).read().splitlines() if l.strip()] for t in times]

  # A job that was killed before its first snapshot leaves the directory behind
  # with an empty time.dat; it contributes nothing and must not divide by zero.
  live = [i for i, ls in enumerate(lines) if ls]
  if not live:
    print ('  SKIPPED %s -- %s hold no records' % (dst, ', '.join(srcs)))
    return None
  if len(live) < len(srcs):
    print ('NOTE: %s hold no records -- merged from the rest'
           % ', '.join(srcs[i] for i in range(len(srcs)) if i not in live))
  srcs  = [srcs[i]  for i in live]
  lines = [lines[i] for i in live]

  tts  = [np.array([float(l.split()[0]) for l in ls]) for ls in lines]
  keep = cut_indices(tts)

  if not os.path.isdir(dst):
    os.makedirs(dst)
  names = sorted(set(os.path.basename(p)
                     for s in srcs for p in glob.glob(os.path.join(s, '*'))
                     if p.endswith('.dat')))
  for name in names:
    with open(os.path.join(dst, name), 'wb') as o:
      for src, ls, k in zip(srcs, lines, keep):
        path = os.path.join(src, name)
        if not os.path.exists(path):
          continue
        if name == 'time.dat':
          o.write(('\n'.join(ls[:k]) + '\n').encode())
        else:
          size = os.path.getsize(path)
          if size % len(ls):
            raise SystemExit('%s: %d bytes is not a whole number of %d records'
                             % (path, size, len(ls)))
          with open(path, 'rb') as f:
            copy_n_bytes(f, o, size // len(ls) * k)
  return keep


#==============================================================================
# merge the numbered backups
#==============================================================================
if restart_num == 'all':
  import numpy as np

  present = sorted(set(n for name in filenames + dirnames for n in numbers(name)))
  if not present:
    raise SystemExit('nothing to merge -- no numbered backups in this directory')
  print ('merging %d backups: %s\n'
         % (len(present), ', '.join('*%d' % n for n in present)))

  for filename in filenames:
    if not (filename.endswith('.nc') and os.path.exists(filename)):
      continue
    last = backups(filename)
    if not last or os.path.getsize(last[-1]) != os.path.getsize(filename):
      print ('CAUTION: %s in the run directory is not the same size as %s.\n'
             '         The last job looks unarchived -- run "%s" first,\n'
             '         or it will be left out of the merge.\n'
             % (filename, last[-1] if last else 'any numbered backup',
                os.path.basename(argv[0])))
      break

  for filename in filenames:
    srcs = backups(filename)
    if not srcs:
      continue
    absent = [n for n in present if n not in numbers(filename)]
    if absent:
      print ('NOTE: %s is missing from %s -- merged from the rest'
             % (filename, ', '.join('*%d' % n for n in absent)))
    dst = filename + ALL
    if os.path.exists(dst):
      os.remove(dst)

    if filename.endswith('.in'):
      # the input of the first job is the one that describes the run
      shutil.copy2(srcs[0], dst)
      print ('%-28s <- %s' % (dst, srcs[0]))
    elif filename.endswith('.nc'):
      keep = merge_nc(srcs, dst)
      if keep is not None:
        print ('%-28s <- %s records' % (dst, ' + '.join(map(str, keep))))
    elif filename.endswith('.dat'):
      keep = merge_table(srcs, dst)
      print ('%-28s <- %s rows' % (dst, ' + '.join(map(str, keep))))
    else:
      merge_text(srcs, dst)
      print ('%-28s <- concatenated' % dst)

  for dirname in dirnames:
    if dirname == 'restart':
      # a checkpoint is a state, not a history: the highest-numbered restartN
      # is already the one the run ended on, so there is nothing to merge.
      continue
    srcs = backups(dirname)
    if not srcs:
      continue
    absent = [n for n in present if n not in numbers(dirname)]
    if absent:
      print ('NOTE: %s/ is missing from %s -- merged from the rest'
             % (dirname, ', '.join('*%d' % n for n in absent)))
    dst = dirname + ALL
    if os.path.isdir(dst):
      shutil.rmtree(dst)
    keep = merge_bindir(srcs, dst)
    if keep is not None:
      print ('%-28s <- %s records' % (dst + '/', ' + '.join(map(str, keep))))

  print ("\nmerged into *%s -- set restart_num to '%s' in load.py" % (ALL, ALL))


#==============================================================================
# archive the live output under a number
#==============================================================================
else:
  if restart_num is None:
    restart_num = next_number()

  clash = ([f for f in filenames if os.path.exists(f + restart_num)] +
           [d for d in dirnames  if os.path.exists(d + restart_num)])
  if clash:
    raise SystemExit('There is already a backup with restart_num = %s (%s)'
                     % (restart_num, ', '.join(c + restart_num for c in clash)))

  for filename in filenames:
    if not os.path.exists(filename):
      print ('skipping ' + filename + ' (not written by this run)')
      continue
    print (('moving ' if move else 'copying ') + filename)
    if move:
      shutil.move(filename, filename + restart_num)
    else:
      shutil.copy2(filename, filename + restart_num)

  for dirname in dirnames:
    if not os.path.isdir(dirname):
      print ('skipping ' + dirname + '/ (not written by this run)')
      continue
    # restart/ is always copied: the run this belongs to may still be continued
    # from it, and the archived copy is what makes this job replayable.
    moving = move and dirname != 'restart'
    print (('moving ' if moving else 'copying ') + dirname + '/')
    if moving:
      shutil.move(dirname, dirname + restart_num)
      os.makedirs(dirname)
    else:
      shutil.copytree(dirname, dirname + restart_num)

  print ('\narchived as *' + restart_num)

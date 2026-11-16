# MHD_COMP_ISOTH: Calliope's paper tests run with Gallope on TSUBAME4.0

These scripts rerun the 589 cases of Calliope's `tests/MHD_COMP_ISOTH` (folders `Euler1D`, `MHD1D`,
`Smooth_wave`, `OT2D`, `Turb3D`) with `bin/gallope-MHD_COMP_ISOTH`. The run directories have the layout that
Calliope's analysis scripts read, so the scripts work on the Gallope runs unchanged.

| file | role |
|---|---|
| `make_cases.py` | sets up the run directories, writes the generated inputs, checks the translated namelists, lists the cases |
| `job-tsubame.sh` | Grid Engine job: runs a list of cases concurrently on the GPUs of one node (port of Calliope's `stages/job.pbs-ps-fes-dev`) |
| `gallope_layout.py` | converts restart files and `out3d/` between the Gallope and the Calliope layout |
| `ngpu_map.md` | GPUs of every case (`make_cases.py ngpu-map`) |
| `judge-tsubame.sh` | Grid Engine job on CPUs: Calliope's judgments of the chosen folders and the comparison with the s13d outputs |
| `check_mhd1d.py` | Calliope's checks D1, D4 and CV1 of `MHD1D` (called by `judge-tsubame.sh`) |
| `compare_s13d.py`, `compare_s13d_all.py` | compare a Gallope run with Calliope's s13d output of the same case; run on a list and tabulate |
| `plot_overlay.py` | Calliope's paper figures with the Gallope results drawn over them (`fig/`, not committed) |
| `fft_roundtrip.F90`, `Makefile.fft`, `setup_g1.sh` | round trip of cuFFTMp with a dimension of size 2 (gate G1 of `validation_tsubame.md`) |
| `validation_tsubame.md` | the results on TSUBAME4.0 |

Each script documents its options in its header or docstring.

## Layout

The runs have the layout of Calliope's `tests/MHD_COMP_ISOTH` on ps1:

```
T/<Folder>/<stage>/<batch>/[<sub>/]<case>/    the case, e.g. T/Euler1D/stage1_es/shock/Sod_N256/
T/<Folder>/...                                symbolic links to the summaries in runs/ (as in Calliope's folders)
T/runs/<stage>/<batch>/[<sub>/]<case>         relative symbolic link to the case, e.g. ../../../Euler1D/stage1_es/shock/Sod_N256
T/runs/...                                    summaries of the batches, judge/, job logs, case lists
```

`T/runs` is the run root (RUNROOT). Its paths are those under `tests/MHD_COMP_ISOTH/runs/` that Calliope's
analysis scripts use (`stages/paths.py`). The cases are first set up and run under RUNROOT and are moved to
their folders with `make_cases.py folders`; every command of `make_cases.py` and the jobs take the links for
cases. On TSUBAME:

```
GALLOPE=~/gallope_isoth                    # built: bin/gallope-MHD_COMP_ISOTH
CALLIOPE_ROOT=~/work/calliope_1d924a9
T=/gs/bs/tgh-26IDS/ua07691/MHD_COMP_ISOTH
RUNROOT=$T/runs                            # also /gs/bs/tgh-26IDS/ua07691/isoth_runs (a symbolic link, the path in the records)
PY=~/work/venv-isoth/bin/python            # numpy, and what Calliope's stages/regen_inputs.py imports
MC=$GALLOPE/tests/MHD_COMP_ISOTH/make_cases.py
ln -s $RUNROOT $CALLIOPE_ROOT/tests/MHD_COMP_ISOTH/runs     # if not there yet
```

`$GALLOPE/tests/MHD_COMP_ISOTH/{Euler1D,MHD1D,Smooth_wave,OT2D,Turb3D,runs}` are symbolic links to `$T/...`,
so the runs also appear next to the scripts, as on ps1. The data stay in the group storage `/gs/bs/tgh-26IDS`
(Lustre project 103190): a rename into another project (`~/work` is project 7691) or another file system
fails, and `folders` stops there instead of copying.

## Setup (once)

```
$PY $MC setup --calliope $CALLIOPE_ROOT --runroot $RUNROOT
$PY $MC ic    --calliope $CALLIOPE_ROOT --runroot $RUNROOT --python $PY
$PY $MC check --runroot $RUNROOT
$PY $MC folders --calliope $CALLIOPE_ROOT --runroot $RUNROOT            # dry run; add --apply to move
```

- `setup` copies `calliope.in`, `nproc`, `case.json`, `generated.md5`, `t_end` and writes `ngpu` (kept if it
  exists, so a hand-edited value survives), `calliope_case` and `calliope_extra` (`exe`/`fftw_estimate` of
  Calliope, which the job ignores), and `terminate` in the restart-chain test
  `stage1_es/smoke/collide_U2_N512_est_chain` before its first run. It never overwrites a file; a different
  file is reported as a conflict. Both commands can be repeated.
- `ic` writes `profile1d.dat` or `ic/` with Calliope's `stages/regen_inputs.py` and compares them with
  `generated.md5` (the md5 of the files of the Calliope runs). A different md5 is reported (exit status 1)
  but does not stop the job; it means that numpy of the venv rounds differently from the one on ps1.
- `check` translates every `calliope.in` as the job does and checks every group and key against the
  namelists declared in Gallope's sources. Gallope reads namelists with `iostat`, so an unknown key would
  otherwise leave its whole group at the defaults without a message. It must print `check: PASS`.
- `folders` moves every case that is a directory in RUNROOT to `T/<calliope_case>` and leaves the relative
  link in its place, then makes the links of Calliope's folders to the summaries. It moves nothing if a
  target path exists (a conflict, exit status 1) and skips what is already moved, so it can be repeated, e.g.
  after more cases were set up. `--apply` records the moves in `T/attic/folders_<time>/` with a
  `rollback.sh` that undoes them (`DRY=0 sh rollback.sh`; without `DRY=0` it prints what it would do). On
  TSUBAME it was applied after all runs had ended (`validation_tsubame.md` still names the paths under
  `isoth_runs`, which resolve through the links).

## Running

Submit from `RUNROOT` (the job runs in the submission directory and the case list is relative to it). The
script submits itself; do not `qsub` it.

```
cd $RUNROOT
$PY $MC list --runroot . --folder Euler1D --status todo > euler1d.txt
CASEFILE=euler1d.txt RESOURCE=node_q H_RT=4:00:00 PYTHON=$PY $GALLOPE/tests/MHD_COMP_ISOTH/job-tsubame.sh
```

`list` orders the cases by decreasing cost (nlx nly nlz nstep), and the job starts them in that order on the
first free GPUs: a case with `ngpu` = n waits for n GPUs that each run fewer than `OVERSUB` cases. A case with
`ngpu` above the GPUs of the job is skipped (`list --max-ngpu` leaves it out). 546 cases run on 1 GPU, 41 on
2 and 2 on 4 (`ngpu_map.md`).

Run one case per GPU (`OVERSUB` = 1, the default). Users cannot start the MPS daemon on TSUBAME4.0, so
cases that share a GPU are time-sliced, and a step of these small grids waits on about 2,500 stream
synchronizations (cuFFTMp and the device reductions). Measured on one H100 (job 9404161, Euler1D
`U2_N1024_fs_nt6`, 200 steps): 0.026 s/step alone, 0.87 s/step for each of 4 cases at once, i.e. 8.4 times
less throughput with `OVERSUB` = 4. For many 1-GPU cases, split the list into disjoint lists (two jobs must
never share an unfinished case) and submit one `node_q` job per list; `node_q` jobs start much sooner than
`node_f` jobs. The cases with `ngpu` = 2 or 4 go to a `node_h` or `node_f` job:

```
$PY $MC list --runroot . --folder MHD1D --status todo --max-ngpu 1 > m1.txt
for k in 1 2 3 4 5 6; do awk -v k=$k 'NR % 6 == k % 6' m1.txt > m1_q$k.txt; done    # round robin
for k in 1 2 3 4 5 6; do CASEFILE=m1_q$k.txt RESOURCE=node_q H_RT=4:30:00 NAME=mq$k PYTHON=$PY \
    $GALLOPE/tests/MHD_COMP_ISOTH/job-tsubame.sh; done
```

| cases | RESOURCE | s/step on 1 GPU |
|---|---|---|
| `Euler1D`, `MHD1D`, `Smooth_wave` (1 GPU) | node_q | 0.026 (1D, N = 1024) |
| `OT2D` | node_q, node_h | 0.031 (256^2) |
| `Turb3D` | node_q .. node_f | up to 256^3 on 1 GPU; set `BUDGET` well below `H_RT` (see below) |

`BUDGET` (minutes, default `H_RT` - 3) is when the runs stop: each run gets `max_wall_time` = the time left,
writes its restart files and ends as `incomplete`; the job continues it from `restart/` while at least 3
minutes are left, and the next submission of the same cases continues it again (the earlier output goes to
`chunk_<k>/`, which the analysis scripts concatenate). After each run the job converts `restart_g/` and
`out3d/` (below); `BUDGET` must leave time for that.

Resubmit with a fresh `list --status todo`, which lists the cases with no `exit_status`, the `incomplete`
ones and the finished ones whose `out3d/` conversion was cut by the end of a job (the job finishes it on the
CPUs). Other cases are skipped; to rerun one, remove `exit_status`, `chunk_*`, `restart`, `restart_g`
and the outputs.

`exit_status` of a case:

| value | meaning |
|---|---|
| (no file) | not run yet |
| `incomplete` | stopped by `BUDGET` or by the file `terminate`, or a `t_end` case that ended before `t_end`; continued |
| `exit 0` | done |
| `exit N` | mpirun returned N |
| `dt_changed` / `dt_stop` | `reset_dt` changed dt in a case without `t_end` / stopped the run (`reset_method = 'stop'`) |
| `exit t_end: ...` | the last row of `cfl.dat` does not end at t of `restart/time.dat` (`t_end` cases) |
| `exit chunk_1 exists` | no `exit_status` but a `chunk_1/`: a broken chain |
| `exit layout: ...` | a layout conversion failed after the run; its message is in `layout.log` |

`DRYRUN=1` prints what the job would do for each case (GPUs, start or continuation, conversions, and the
difference of `calliope_run.in` from `calliope.in`) without submitting, loading modules or changing a file;
it runs on any machine:

```
DRYRUN=1 CASEFILE=euler1d.txt RESOURCE=node_f bash $GALLOPE/tests/MHD_COMP_ISOTH/job-tsubame.sh
bash $GALLOPE/tests/MHD_COMP_ISOTH/job-tsubame.sh translate stage1_es/shock/Sod_N256/calliope.in
```

## What differs from a Calliope run

- Input: `calliope_run.in` is `calliope.in` with `&mpi_settings` and every `pruned = ...` removed,
  `restart_dir = './restart_g/'`, and the job's edits of `max_wall_time`, `init_type`, `nstep` and `dt` as
  in the PBS job. Comments are kept.
- Dealiasing: Gallope has only the 2/3 rule; Calliope's `pruned = T` (the 3/2 rule, padded transforms of
  every field) does not exist. The two cases `stage1_es/scan/{Sod,U2}_N1024_pruned` therefore run the 2/3
  rule and repeat `*_N1024_fs_nt6` (same results to every digit), and Calliope's comparison of the two
  rules (`REPORT_stage1.md` section 9) is not reproduced.
- Restart files: Gallope reads and writes `restart_g/` in its own layout (r2c along z). A case with `ic/`
  starts from `ic/` converted to `restart_g/`; after every run `restart_g/` is converted to `restart/`,
  so `restart/` always holds the Calliope layout (r2c along y) that the chain and the analysis read.
- `out3d/`: converted in place to the Calliope layout after every run (`gallope_layout.py
  out3d-to-calliope`); `out3d/time.dat` is unchanged. The conversion only moves and conjugates modes. It
  needs about one extra field file of disk space while it runs and is resumable: an interrupted conversion
  is redone from the backup `f.dat.gallope` on the next submission. It ran at about 0.7 GB/s on a 128^3
  test file on a laptop; the largest `out3d/` (about 28 GB, `stage5/res/turb4_N128_a`, `turb4b_N128`) may
  take minutes on Lustre (not measured).
- GPUs: `ngpu` replaces Calliope's `nproc`; `exe` and `fftw_estimate` of Calliope are not used.

## Analysis

Calliope's analysis and figure scripts run unchanged from `CALLIOPE_ROOT`, as for Calliope's runs; they find
the runs through the `runs` symlink. `judge-tsubame.sh` runs the judgments of `validation_tsubame.md` as one
CPU job (the analyses read tens of GB, and a login node stops a process after 10 minutes):

```
cd $RUNROOT
FOLDERS="Euler1D MHD1D Smooth_wave" CALLIOPE_ROOT=$CALLIOPE_ROOT PYTHON=$PY REF=/gs/bs/tgh-26IDS/ua07691/s13d \
    H_RT=2:00:00 $GALLOPE/tests/MHD_COMP_ISOTH/judge-tsubame.sh
```

Each folder writes `judge/<folder>.log` (every command, its output and its exit status). `REF` is a copy of
Calliope's s13d outputs, `REF/<folder>/<case>`, without `out3d/`, the NetCDF files and the executables; it
replaces the checks of Calliope that compare two Calliope versions (`validation_tsubame.md`).

The figures are drawn on a laptop with a Calliope worktree at the same commit:

```
python3 plot_overlay.py fetch OT2D       # copies the inputs of the figures from TSUBAME to runs_gallope/
python3 plot_overlay.py plot OT2D        # -> fig/OT2D/*.pdf
```

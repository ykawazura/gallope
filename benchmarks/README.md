# benchmarks/ — the timings behind Section 4 of COMPHY-D-26-00110

Figures 1-5 of the paper are timings: of the distributed FFT alone, and of a
whole Gallope or Calliope step.  What is here is the raw stdout every one of
those runs wrote, the inputs and job scripts that produced it, and the script
that turns it into the PDF the manuscript includes.  Nothing is summarised on
the way in — a point on a curve is a number parsed out of a file in this tree.

The comparison runs of Section 5 are **not** here.  They are pairs, Gallope
against Calliope on the same initial condition, and they live together in their
own repository:

> **[ykawazura/gallope-calliope_comparison](https://github.com/ykawazura/gallope-calliope_comparison)**
> — Figures 6 and 7 and Table 2 (deterministic comparison), Figure 8 (nonlinear MRI).

## Figure by figure

| figure | data | plotting script | output |
|---|---|---|---|
| Fig. 1 — FFT strong scaling, Miyabi-G | `p3dfft_vs_cufftmp/miyabi/strong scaling/202607/` | `plot_combined_202607.py` | `cuFFTMp-scaling.pdf` |
| Fig. 2 — FFT strong scaling, TSUBAME 4.0 | `p3dfft_vs_cufftmp/tsubame/` | `plot_combined_tsubame.py` | `cuFFTMp-scaling-tsubame.pdf` |
| Fig. 3 — step time, Miyabi-G | `gallope/miyabi/strong scaling/202607/` | `plot_202607.py` | `1024-202607.pdf` |
| Fig. 4 — step time, TSUBAME 4.0 | `gallope/tsubame/strong scaling/202607/` | `plot_tsubame.py` | `gallope-scaling-tsubame.pdf` |
| Fig. 5 — weak scaling, both machines | `p3dfft_vs_cufftmp/weak-scaling-combined/` | `plot_weak_combined.py` | `weak-scaling-combined.pdf` |

Each output PDF is committed beside its script and is byte-identical to the
file the manuscript `\includegraphics`-es (`cuFFTMp-scaling.pdf`,
`cuFFTMp-scaling-tsubame.pdf`, `gallope-scaling.pdf`,
`gallope-scaling-tsubame.pdf`, `weak-scaling.pdf` respectively), so a re-plot
can be checked against the figure as printed.

`plot_weak_combined.py` is the only script that reaches outside its own
directory: it reads the Miyabi half from `../miyabi/weak scaling/202607/` and
the TSUBAME half from `../tsubame/`.  The directory names are therefore load
bearing and were kept exactly as they were on the machine that produced them.

## 202601 and 202607

Two campaigns are kept, not one.

* **`202601`** is the measurement as **submitted**: one submission per
  configuration, run directories named `node16`, `node32`, ....
* **`202607`** is the **re-measurement** for R1, under the protocol Section 4.1
  now states: three *independent* job submissions per configuration, named
  `node16-1`, `node16-2`, `node16-3`, with the figures plotting the median and
  whiskers spanning the observed minimum and maximum.  Figures 1 and 3 of the
  revised manuscript are the 202607 versions; the 202601 ones are kept so that
  the difference between the two is inspectable rather than asserted.

The TSUBAME 4.0 measurements are new in R1 and exist only as 202607.

## Protocol, and where it is implemented

Section 4.1 (*Benchmarking protocol*) describes what these runs do; the
correspondence is:

| the paper says | here |
|---|---|
| double precision; forward R2C then inverse C2R on the same field | `p3dfft_vs_cufftmp/_src/{cufftmp-slab,cufftmp-pencil,p3dfft}/` |
| 3 untimed warm-ups, then 10 timed transforms, max over ranks of the per-transform mean | same sources; the per-iteration times are printed in every stdout file |
| step time averaged over steps after the first, excluding initialization and I/O | Gallope/Calliope's own timer block, at the foot of each `*.out.std` |
| three independent submissions, median plotted, min-max whiskers | `_miyabi_src/bench-202607/scripts/campaign.sh`, `_tsubame_bench/campaign.sh` |
| NVHPC 25.9, `nvfortran` + OpenACC, `-O3 -fast`, cuFFTMp; CPU baseline same toolchain with HPC-X and P3DFFT | `_miyabi_src/bench-202607/bin/arch-*.in`, build logs in `logs/` |
| 1 rank/node on Miyabi-G, 72 ranks/node for the CPU baseline; 1 rank/GPU and 160 ranks/node on TSUBAME 4.0 | the generated `job.pbs` / `job.sh-tsubame`, see `scripts/setup_runs.sh` |

`_miyabi_src/bench-202607/` is the campaign as it stood on Miyabi: `scripts/`
builds every binary (`build_all.sh`), lays out every run directory with its
input and job script (`setup_runs.sh`, `setup_fft.sh`), submits them and
resubmits what died (`campaign.sh`, `recover.sh`, `watchdog.sh`), `inputs/`
holds the two frozen input files, `lists/` the run lists fed to the driver, and
`logs/` the build and campaign logs.  `_tsubame_bench/` is the TSUBAME
counterpart — one driver plus `verify_tsubame.py`.  Note that
`_tsubame_bench/campaign.sh` labels its targets `fig1`/`fig2`/`fig3` after the
figure order of an earlier draft; in the manuscript as it now stands these are
Fig. 2, Fig. 4 and Fig. 5(b).

Four Calliope build configurations (`bin/arch-calliope-{A,B,C,D}.in`: with and
without `-O3 -fast`, against a stock and an `-O3` P3DFFT, with and without
`-Mnocache_align`) were screened at 64 nodes before the production runs, so
that the CPU baseline would not be handicapped by its build.  The four builds
and their build logs are here; the screening job's own stdout was left on
Miyabi and is not.

## The Gallope these timings are of

Figures 3 and 4 are not timings of the tip of `main`.  Both campaigns build from
one frozen tree, archived here entire:

    _gallope_src/           src/ arch/ batch/ inputs/ Makefile Makefile.in

`_miyabi_src/bench-202607/scripts/build_all.sh:9` and `setup_runs.sh:7` both set
`G=$ROOT/gallope_performance`, and `_tsubame_bench/campaign.sh` calls that same
tree's `batch/job.sh-tsubame`, so Miyabi-G and TSUBAME 4.0 ran **the same
solver source**.  That was the point: TSUBAME exists in this paper to separate
"twice the GPUs" from "the all-to-all left the node", and the separation is only
clean if nothing else moved.  `_gallope_src/batch/job.sh-tsubame:22-24` says so
in its own words.

The tree is the January 2026 state of the code.  Against the released `main` it
differs by about 3300 lines across 23 source files — `mp.F90` (+236 lines),
`mpiio.F90` (+76), `model/MHD_INCOMP/fields.F90` (+106) and
`model/MHD_INCOMP/advance.F90` are the larger ones — because the intervening
work added the process-grid generalisation, the shearing-box pressure fix that
Section 5.2 depends on, and diagnostics.  None of it was in flight when Section
4 was measured, and re-measuring the whole of Figures 1-5 against a moving code
would have cost more machine time than it bought.  So the honest record is both:
`main` at its tag for the code the paper describes, and `_gallope_src/` for the
exact tree its Section 4 numbers came out of.

`_gallope_src/src/` carries `model/MHD_INCOMP/` only.  The benchmark solves
incompressible MHD and nothing else, and the RMHD model was not built for it.

Two `arch/` files sit beside it because Fig. 2 is a two-build figure:
`tsubame-sep.in` and `tsubame-uni.in` differ by `-gpu=mem:unified`.
`_tsubame_bench/campaign.sh` explains why only the separated build was run there.

## Re-plotting

Each script is run from its own directory and needs only Python 3, NumPy and
Matplotlib, plus a LaTeX installation (`text.usetex` is on, via
`latex_preamble.py`):

```bash
cd "gallope/miyabi/strong scaling/202607"
python3 plot_202607.py
```

Every script prints, per configuration, the median, min, max and spread it
derived, and prints every file it *skipped* and why — a header that does not
match the grid or the rank count the directory name claims is dropped rather
than plotted.  The printed table is the check that a re-plot read the same data
as the published figure.

## Runs on disk that are deliberately not plotted

The campaigns reach further than the figures do, and the surplus is kept:

* **TSUBAME 4.0, 32 nodes (128 GPUs).** Both the solver and the transform turn
  over there, and the weak-scaling panel would otherwise quote a different
  reach from the strong-scaling ones, so every TSUBAME figure stops at 16 nodes
  (64 GPUs).  The 32-node runs stay on disk; each script prints the skip.
* **`384^3-202607/` on TSUBAME** — the first TSUBAME campaign, which bought
  3- and 12-GPU points with a grid that is not a power of two, together with
  the rank-placement control (`gpu4-node2`, `gpu4-node4`, `gpu8-node4`,
  `gpu8-node8`: the same GPU count spread over different node counts).  The
  published Fig. 4 is the `512^3` campaign; `plot_tsubame.py` says why.
* **`768^3` and `1536^3` FFT runs on TSUBAME** — the same idea at the transform
  level, and the measurement of why it was abandoned;
  `plot_combined_tsubame.py` carries the arithmetic.
* **`p3dfft_vs_cufftmp/miyabi/weak scaling/202601/`** — superseded by 202607
  and by the combined two-panel figure.

## What is not here

The compiled binaries.  `p3dfft_vs_cufftmp/_src/` holds the sources, Makefiles
and job scripts of the three FFT micro-benchmarks for all three machines they
were built on, but the one binary that had been sitting in that tree (106 MB)
is excluded by its own `.gitignore`.  The Gallope and Calliope binaries were
never in this tree either — `_gallope_src/` is the Gallope source they were
built from, and Calliope is
[ykawazura/calliope](https://github.com/ykawazura/calliope).

No field data.  These runs write nothing but stdout — the input files set the
step count low and leave diagnostics and I/O off — which is why the whole of
Section 4 fits in 11 MB.

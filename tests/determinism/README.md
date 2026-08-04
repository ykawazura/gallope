# Deterministic Gallope ↔ Calliope comparison

Answers Comment 1.6(a) of the COMPHY-D-26-00110 review: *"a deterministic
comparison from identical initial fields"*. The published MRI comparison is
statistical — two independent realizations of a chaotic flow, agreeing in the
saturated statistics but not pointwise. This test is the deterministic
counterpart.

**Status: run on Miyabi-G. The test found a real bug** — see *What it found*
below — which is now fixed in both codes (`calliope 92f4047`, `gallope
4eea7cf`). All numbers quoted here are post-fix.

## What it found

With the shearing box on, both codes converged at **first** order in `dt`
rather than third, and the two codes differed by `1e-6` after a single step
instead of by `1e-16`. The cause was the pressure in `get_ext_terms`.

In shearing coordinates the constraint is `kxt·u = 0` with
`kxt = kx + q·tsc·ky`, so `kxt` itself moves, `d(kxt)/dt = q·ky`, and
differentiating the constraint gives

```
kxt·dux/dt + ky·duy/dt + kz·duz/dt = -q·ky·ux
```

whose solution is `p = -i/k²(t) · [ K + 2·kxt·uy - (2-2q)·ky·ux ]` with
`K = kxt·nl_x + ky·nl_y + kz·nl_z`. The coded expression carried a spurious
`i` on the shear term and the wrong coefficient on `ky·ux`, so the explicit
right-hand side was not divergence free.

It survived undetected because the `div_free` projection at the end of every
step removed the violation: the solution stayed solenoidal (`is_div_free`
printed `0.000E+00` throughout) and the energy budget still closed, so every
diagnostic the codes carry looked healthy. What was lost was the temporal
order — an O(1) constraint violation cleaned up once per step costs exactly
one order in `dt` per step, which is why third order became first.

The lesson generalizes: **a per-step projection can hide an O(1) error in the
right-hand side, and the only diagnostic that sees it is self-convergence
under `dt` refinement.** Code-to-code agreement would not have found it
either — both codes had the same bug, and reproduced each other's wrong
answer to seven significant figures.

## Results after the fix

Combined relative field difference over the six evolved fields, 27 runs, gate
met in every one:

| steps | t | no shear | control | shear q=3/2 | control | shear, before the fix |
|---|---|---|---|---|---|---|
| 1    | 1e-3 | 1.958366e-16 | 0 | 1.958408e-16 | 0 | 1.000016e-06 |
| 10   | 1e-2 | 1.977629e-16 | 0 | 1.978066e-16 | 0 | 9.999484e-06 |
| 100  | 1e-1 | 2.644458e-16 | 0 | 2.605155e-16 | 0 | 1.000074e-04 |
| 1000 | 1    | 1.551102e-15 | 0 | 1.402888e-15 | 0 | 1.085778e-03 |
| 2000 | 2    | —            | — | 2.678303e-14 | 0 | 4.750540e-03 |

`control` is Gallope on 1 GPU against Gallope on 4, and is exactly zero in all
nine configurations.

Self-convergence at fixed `t = 0.1` over `dt = 1e-3, 5e-4, 2.5e-4`:

| series | gallope | calliope | ratio |
|---|---|---|---|
| shear off             | 2.272099e-10 → 2.840388e-11 | identical to 7 digits | 8.00 |
| shear on, nonlinear on| 2.471032e-10 → 3.088888e-11 | 3.088887e-11          | 8.00 |
| shear on, nonlinear off| 9.582121e-12 → 1.198013e-12 | identical to 7 digits | 8.00 |

Before the fix the sheared ratio was 2.00. The code-to-code difference in the
sheared series is now flat in `dt` — 2.61e-16, 2.59e-16, 2.69e-16 — i.e. pure
roundoff, no longer a discretization difference.

Solenoidal error on the reconstructed cube, evaluated on `kx(t)`: ≤1.53e-24
(no shear), ≤6.59e-24 (shear). Energy-balance residual, identical between the
two codes to every printed digit: 0.19 % of the largest budget term without
shear and 1.4e-4 % with it; against `<D>`, 0.57 % and 0.44 %.

## Why not a shared restart file

The obvious approach — run one code, hand its restart to the other — does not
work, and the reason is worth recording because it is invisible from the file
names.

The two codes take the real transform along **different axes**:

| | real transform | global Fortran array |
|---|---|---|
| Gallope  | `z` | `(nkz, nky, nkx)` = `(nlz/2+1, nly, nlx)` |
| Calliope | `y` | `(nkx, nkz, nky)` = `(nlx, nlz, nly/2+1)` |

(`gallope_dev/src/grid.F90:76-78`, `calliope_dev/src/grid.F90:64-66`.) These
are **not** transposes of one another: each holds a different half of the same
Hermitian cube. They use the same `mpiio_write_one` helper and the same
`MPI_ORDER_FORTRAN` layout, so the files look interchangeable and are not.

So instead both codes are started from the *same analytic initial condition*,
and the comparison is done in post-processing on the reconstructed full cube.

## Protocol

**Initial condition — 3D Orszag–Tang (`init_type = 'OT3'`)**

```
ux = -sin(y+z)    uy = sin(x+z)    uz = 0
bx = -sin(y+z)    by = sin(2x+z)   bz = 0
```

Exactly solenoidal, since no component depends on its own coordinate, and
fully three-dimensional.

Not `OT2`: that vortex is independent of `z`, so it exercises no `kz ≠ 0`
dynamics at all and a code could reproduce it exactly while being wrong in the
third direction. Calliope already provides `OT3`
(`src/model/MHD_INCOMP/fields.F90`); the matching `init_OT3` was added to
Gallope for this test.

**Conditions**

- 128³, box `2π × (1, 1, 1)`, `dt = 1e-3` held fixed, `driven = F`,
  `nonlinear = T`, `eSSPIFRK3`
- `nu = eta = 0`, `nu_h = eta_h = 1e2`, `nu_h_exp = eta_h_exp = 4`,
  identical in both inputs
- `cfl = 2.0`, which makes `dt_cfl ≈ 2·(2π/128)/|u|max ≈ 5e-2 ≫ dt`, so the
  adaptive reset never fires
- Calliope needs **`pruned = .false.`**. With the default `pruned = .true.`
  its spectral grid is `int(2·nl/3)+2` on a side
  (`calliope_dev/src/grid.F90:301-309`) and the arrays are not mode-for-mode
  comparable. `compare.py`'s file-size check rejects this case explicitly.
- Two series: `shear = .false.` (shared solver core) and
  `shear = .true., q = 1.5` (adds the remap; `tremap = ly/(lx·q) = 0.667`, so
  `nstep = 2000` crosses three remaps)
- Terminal restart saved at `nstep ∈ {1, 10, 100, 1000}` (and 2000 with shear)
- Gallope on 1 and on 4 GPUs

**Dealiasing is provably identical** and needs no numerical check: both codes
set the cut-off to `2π(n_j/2)/L_j` (`gallope grid.F90:172` `kx_max`,
`calliope grid.F90:124` `kx_max_noprune`) and both `dealias.F90` zero every
mode with `|k| ≥ (2/3)k_max`. Same threshold, same rule, same retained set.

## Gate

`restart/time.dat` must give `tt = nstep·dt` **exactly** in both codes. If it
does not, the adaptive step fired somewhere and the two runs are at different
times, which makes the field comparison meaningless. Both job scripts print
`time.dat` after every run for this reason. `cfl.dat` records
`tt, dt_cfl, max_vel_{x,y,z}` every step if a closer look is needed.

## What is reported

`compare.py` reconstructs the full `nlx × nly × nlz` complex cube from each
code by Hermitian expansion and compares mode by mode — no inverse transform
and no normalization convention enters, so the comparison is exact. It prints:

1. the restart-time gate;
2. the solenoidal error `Σ|k·f| / (Σ|f|·Σ|k|)` on the reconstructed cube — an
   axis-mapping self-check, which must be at roundoff;
3. the relative field difference `‖f_A − f_B‖₂/‖f_B‖₂` per field and combined;
4. the spectral energy `W = ½Σ_k(|u_k|²+|b_k|²)`.

Reported **as a function of step count**, not at one time. In a chaotic system
the meaningful statement is not that the difference stays small — it cannot —
but that it starts at machine precision and then grows at the rate set by the
flow. The 1-GPU vs 4-GPU pair is the control: it bounds what is attributable
to summation order alone.

The solenoidal error quoted in the response letter comes from each code's own
`is_div_free`, printed to `out.std`, not from `compare.py` — the point being
that it is the codes' own diagnostic and not a post-processing invention. The
two codes sum over their own half spaces, so the two printed values are
comparable in magnitude but are not the same sum.

## Running it

```sh
# locally, no GPU and no run needed: validates compare.py itself
python3 selftest.py

# on Miyabi-G, from a directory holding the binary and the input
qsub job.pbs-gallope        # 18 runs, 128^3, minutes
qsub job.pbs-calliope       #  9 runs

# then
./compare_all.sh <gallope run dir> <calliope run dir> compare_all.out
```

The three probe scripts refine `dt` at fixed final time `t = 0.1`, which is
what measures the temporal order and is the diagnostic that found the bug
above. They are separate jobs because each isolates a different part of the
right-hand side:

```sh
qsub job.pbs-shearprobe    # shear on,  nonlinear on   -> probe_<dt>_<code>/
qsub job.pbs-orderprobe    # shear off, nonlinear on   -> order_<dt>_<code>/
qsub job.pbs-linorder      # shear on,  nonlinear off  -> lin_<dt>_<code>/

./probe_all.sh probe_all.out
```

`probe_all.sh` reports, for each series, each code's **self-convergence**
(the same code at `dt` and `dt/2`) alongside the code-to-code difference at
each `dt`. Self-convergence is the sharper of the two: it is a property of
one code alone and does not depend on the two codes agreeing, so it catches
an error both codes share.

**`--q` is not optional in the sheared series.** The solenoidal check is a
check on `kx(t) = kx + q·tsc·ky`, not on `kx`; `compare_all.sh` passes
`--q 1.5` for the `shear` series and `--q 0.0` for `noshear` for that reason.
Run with the static `kx` the very same fields return `1e-8` instead of
`1e-24` — small enough to pass an eyeball test, and the wrong quantity.

`selftest.py` builds an exactly-solenoidal Hermitian field on a deliberately
unequal 16×12×8 grid with an unequal box aspect — so a transposed axis cannot
hide behind a cubic grid — writes it in both layouts, and checks that the
round trip is exact, that a single mode perturbed by 1e-12 rises 65× above the
roundoff floor, and that a shifted axis and a swapped vector component are
both caught. It requires no MPI, no GPU and no Fortran.

## `rmhd/` — the same test for the RMHD model

Answers Comment 1.7: *"Gallope is presented as both an incompressible-MHD and
RMHD solver, but only the former is tested."* Same protocol, same
reconstruction, same gate; only the model differs. Everything in `rmhd/` is
self-contained except `compare.py` and `selftest.py`, which live one level up
and are shared.

The two models share only the transform layer, the dealiasing and the time
integrator, so this is not a re-run of the test above. RMHD evolves `omg` and
`psi`; its nonlinearity is a pair of perpendicular Poisson brackets rather
than the divergence of a stress tensor; the guide-field coupling enters
through `nabla_par`; there is no pressure and no projection; and the
dissipation is split perpendicular/parallel with four independent
coefficients.

**`compare.py --model RMHD`** switches the field list to `phi, omg, psi` and
nothing else — the Hermitian reconstruction and the layout handling are shared
verbatim, because both codes' RMHD uses the same `grid.F90` as MHD_INCOMP.

**The solenoidal check is replaced.** `u_perp = z x grad_perp(phi)` is
divergence free by construction, so checking it would verify nothing. In its
place `compare.py` checks `omg + kprp^2*phi` on the reconstructed cube, which
the program itself must satisfy since `omg` and `psi` are the evolved fields
and `phi` is stored alongside them. It is the right substitute because
`kprp^2 = kx^2 + ky^2` mixes exactly the two axes the two codes lay out
differently, so a crossed mapping cannot pass it — `selftest.py`'s negative
control shifts x by one index and the residual goes from `0` to `3.0e-01`.

**Conditions**: 128³, box `2π × (1, 1, 1)`, `init_type = 'OT3'`, `dt = 1e-3`
fixed, `nonlinear = T`, `driven = F`, `eSSPIFRK3`,
`nupe_x = etape_x = 1e2`, `nupe_z = etape_z = 1e1`, all four exponents 4.
`nstep ∈ {1, 10, 100, 1000}`; Gallope on 1 and 4 GPUs, Calliope on 16 Grace
cores; probe at fixed `t = 0.1` over `dt = 1e-3, 5e-4, 2.5e-4`. RMHD has no
shearing box, so there is one series rather than two and `--q` never enters.

### Results

| steps | t | phi | omg | psi | combined | control |
|---|---|---|---|---|---|---|
| 1    | 1e-3 | 2.53e-16 | 2.21e-13 | 2.26e-16 | 1.37e-13 | 0 |
| 10   | 1e-2 | 2.34e-16 | 2.16e-13 | 2.89e-16 | 1.33e-13 | 0 |
| 100  | 1e-1 | 2.33e-16 | 1.86e-13 | 1.22e-15 | 1.18e-13 | 0 |
| 1000 | 1    | 5.10e-15 | 8.78e-14 | 6.86e-15 | 7.20e-14 | 0 |

`omg` is three decades larger than the potentials, and that is a property of
the norm, not of the solution: `omg = -kprp^2*phi`, so `‖δomg‖/‖omg‖` is a
`kprp^2`-weighted relative error of `phi`, and `kprp^2` puts nearly all the
weight on the grid-scale modes, where `phi` sits ~13 decades below its peak
and is therefore at the floor set by the rounding of the transforms — which is
absolute and scales with the largest mode, not the local one. It *decreases*
with step count as the hyperdiffusion removes those modes, which is the
opposite of a divergence between the codes.

Self-convergence at fixed `t = 0.1`:

| | dt 1e-3 → 5e-4 | 5e-4 → 2.5e-4 | ratio |
|---|---|---|---|
| gallope  | 5.866494e-10 | 7.334202e-11 | 8.00 |
| calliope | 5.866492e-10 | 7.334204e-11 | 8.00 |

Third order, and the two codes' truncation errors agree to seven significant
figures. Over the same refinement the code-to-code difference is
`1.176018e-13, 1.176148e-13, 1.176533e-13` — flat to four digits across a
factor of four in `dt`, so it carries no discretization component.

`omg + kprp^2*phi` ≤ `2.5e-18` for either code in every configuration. Energy
`W = ½Σ_k kprp²(|phi|²+|psi|²)` is identical to the last bit after one step and
differs by `2.2e-16, 7.8e-16, 3.2e-15` after 10, 100, 1000 steps.

Unlike the MHD_INCOMP test this one found no solver error.

### Two things worth recording

**`init_OT2` in Gallope's RMHD was broken.** It set `phi` and `psi` but never
`omg`, then copied the zero `omg` into `omg_old`. Since `omg` and `psi` are
the evolved fields (`iomg = 1, ipsi = 2` in `advance.F90`), it started from a
state of no flow. Found by reading the code while adding `init_OT3`, not by
the test. It is used in no published result — the scaling runs use
`init_type = 'random'` with `MODEL = MHD_INCOMP` — and is now fixed.

**Calliope's RMHD `init_OT3` has `yy(j)/x0` where `y0` was evidently
intended**, in the second `psi` term. Gallope's copy transcribes it verbatim
and deliberately: the two codes must start from bit-identical fields or the
comparison means nothing, and correcting one side alone would break that. The
test runs with `lx = ly`, where the two are equal, so nothing here depends on
it. Fixing it is a separate change that must be made in both codes at once.

### Running it

```sh
qsub job.pbs-rmhd-gallope    # debug-g,  4 nodes, ~2 min
qsub job.pbs-rmhd-calliope   # short-g,  1 node,  ~4 min
./compare_rmhd.sh            # -> compare_rmhd.out, 15 comparisons
```

## Files

Tracked: `gallope.in`, `calliope.in`, `job.pbs-gallope`, `job.pbs-calliope`,
`job.pbs-shearprobe`, `job.pbs-orderprobe`, `job.pbs-linorder`, `compare.py`,
`compare_all.sh`, `probe_all.sh`, `selftest.py`, this `README.md`, and
`rmhd/{gallope.in, calliope.in, job.pbs-rmhd-gallope, job.pbs-rmhd-calliope,
compare_rmhd.sh}`.
Untracked: `run_*/`, `probe_*/`, `order_*/`, `lin_*/`, `compare_all.out`,
`probe_all.out`, `rmhd/compare_rmhd.out`, PBS logs, `__pycache__/`.

`gallope.in` and `calliope.in` must stay in step, in both directories. The
only intended differences are Calliope's `&mpi_settings` and its
`pruned = .false.`, plus, in `rmhd/`, Gallope's `write_intvl_nltrans`, which
Calliope's RMHD does not accept; `diff` them after editing either.

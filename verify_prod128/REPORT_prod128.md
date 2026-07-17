# prod128 — Driven KRMHD velocity-space cascade (Phase 1-C-ii, Stage 2)

**Status: INCONCLUSIVE — superseded by the balanced rerun (§21 / `echo128bal`).**

> **Why this run does not settle the echo question.** The Alfvénic forcing here is strongly
> **imbalanced** (`xhl_inj = 0.8` → steady σ_c = 0.967, z⁺²/z⁻² = 58.8). The advecting flow is then
> dominated by a single coherent Elsasser field (z⁺), which is **non-stochastic**, so it cannot drive
> the stochastic echo — the "no echo / forward cascade" reading below is an **artifact of the
> imbalanced setup**, not a valid regime test. Three further inconsistencies confirm the run is not
> a clean measurement: (1) χ ≈ 3.4 (strongly nonlinear — should fluidize) yet the flux is forward,
> a contradiction; (2) the k-integrated Γ(m) is **not constant** (declines 0.95 → 0.57 across
> m = 8–100, 6–9 dΓ/dm sign changes), so "0.98 × injection across m = 20–48" over-claims a plateau
> that is not there; (3) no clean Hermite power law (sub-band slopes −0.52 / −0.14 / −0.66), and the
> Alfvén bath is not fully steady (injection 1.0 vs dissipation 0.64). The telescoping certification
> in §3 (top-m Γ = machine zero) is **still valid** — that is a diagnostic-correctness statement,
> independent of the physics regime. **Rerun with balanced forcing** (`xhl_inj = 0`, σ_c ~ 0,
> longer averaging, verify χ ≳ 1) to re-measure ⟨Γ(m)⟩; see `echo128bal.in` and todo.md §21.
>
> The sections below are retained as the record of what was reported at the time; treat every
> physics conclusion in them as **retracted pending the balanced rerun**.

## 1. Result (結論)

At these parameters the driven KRMHD run is **phase-mixing-dominated**: the k⊥+k_z–integrated Hermite
free-energy flux is a **clean forward cascade**, ⟨Γ(m)⟩ ≈ **0.98 × injection** across the inertial
band (m = 20–48; mean 0.98 ± 0.10), decaying through the ν m⁶ range to machine-zero at the top moment. This is the
**opposite** of the stochastic-echo / fluidization signature (⟨Γ(m)⟩ ≈ 0 in the inertial band).
**The stochastic echo is NOT reproduced at these parameters** — essentially all injected compressive
free energy is carried forward in Hermite space to high-m collisional dissipation (Landau/phase mixing).

This is a physically valid, certified result — the counterpoint to Meyrand et al. (2019)'s fluidized
regime — not a numerical artifact (see §3).

## 2. Run configuration

| item | value |
|---|---|
| input | `echo128_restart.in` (chain-continued from `echo128.in`) |
| jobscript | `job.pbs-echo128` (`regular-g`, `select=64:mpiprocs=1`, `walltime=01:00:00`) |
| grid / moments | 128² × 64, nm = 128 |
| process grid | np = 64 → [P_fft=8, P_m=8, P_s=1] (nm_local=16, nly_local=16) |
| g-forcing | driven=T, m=1 low-k shell kmin=(0,0,1)/kmax=(1,1,2), fix_power, **ene_inj_g = 1.0** |
| Alfvén forcing | Elsasser zppe/zmpe, ene_inj = 1.0, xhl_inj = 0.8, freq (0.9, −0.6), stir_seed = 7 |
| α inputs | β_i = τ = Z = 1, alpha_root = +1, v_th = 1 |
| dissipation | ⊥ hyper μ = 10 (nexp_perp=4 → k⊥⁸); Hermite ν = 1 (nexp_m=6 → m⁶) |
| time stepping | dt = 2e-4, cfl = 0.5 (limiter-verified, todo.md §19), eSSPIFRK3 |
| diagnostics | write_intvl = 0.5, write_hermite_flux = T, save_restart_intvl = 2.0 |

**Why g-forcing (todo.md §20):** passive g has no free-energy source — advection {Φ,g} and streaming
(∂_z+{Ψ,·})S_m both conserve W = ½Σ|g_m|² (proven to 2e-9 in 1-B), so the earlier passive run (§19)
merely measured free decay of the initial condition (retracted). Per Meyrand 2019, compressive free
energy must be injected at **m=1** at low wavenumbers with a constant rate — implemented as `ene_inj_g`.

## 3. Diagnostic validation — the Γ_m flux is trustworthy

The k-integrated Hermite flux `Gamma_m_kint` is a **telescoping invariant**: Σ_m src_m(streaming) = 0
identically, so Γ at the top moment must vanish to machine precision. Certified two independent ways:

- **verify_1Cii cheap gates, 6/6 PASS** (jobs 2351320 / 2351971): telescoping
  `TELE_NL_KINT = 4.24e-16` (≪ 1e-11 target), linear-contrast `SIGN_PEAK = +0.137 > 0` and
  `SIGN_MIN_SIG = +7.1e-3 > 0` (forward and positive across the entire significant band), P_m-invariance
  `PMREL = 5.2e-16`.
- **Current g-forcing production binary** (added after the §18 gates): telescoping survives the m=1
  forcing to machine precision at nm=128, in **both** flux segments — `Gamma_m_kint(m=127) = −2.26e-16`
  in `echo128_seg_2358620.nc` (t10–16), and per-record range `[−1.5e-15, +3.1e-16]` across the steady
  segment `echo128_seg_2358914.nc` (t16.5–21.76). The top-moment flux vanishes every diagnostic step.

## 4. Results (steady window t≈19–21.76, later-half time average)

**Figure:** `report_prod128_fig.png` (4 panels; local, not committed — regenerate with `make_report_fig.py`).

- **(A) Hermite flux — forward cascade.** ⟨Γ(m)⟩: rises to +1.0 at m=1–2 (= injection from the m=1
  forcing), stays **≈ injection across the inertial band m=20–48 (mean 0.98 ± 0.10** — a bump to ~1.11
  near m=24, a shallow dip to ~0.83 near m=50), then a monotone decline through the shaded ν m⁶ range to
  machine-zero at m=127. Forward Hermite (phase-mixing) cascade carrying ~all the injection to
  collisional dissipation. **Not the echo** (echo would give ⟨Γ⟩≈0 in this band).
- **(B) Hermite spectrum.** ⟨W_m⟩ power-law slope ≈ **−0.29** over m=3–40 — shallower than both the
  m^{−1/2} phase-mixing and m^{−1} fluidized references, with a mild bump before the ν m⁶ cutoff.
- **(C) Free energy / Gate 2 — steady state reached.** After the IC transient (t=0.5→2.5: 16.3→14.9),
  W_free rises under the drive to a **peak ≈ 17.66 at t≈19.5**, then turns over and settles: the
  later-half mean is **17.59 ± 0.06 (0.3 % variation)**, last-half slope **−0.07/unit** (no longer
  drifting up). The run has **escaped §19's free decay** and reached a **statistical steady state**
  (W bounded, oscillating about a plateau) — Gate 2 met.
- **(D) Turbulence health.** Anisotropic Elsasser perp spectrum, |z⁺|² > |z⁻|² (imbalanced drive),
  both steepening into the ⊥ hyperviscous range; inertial band roughly ~k⊥^{−5/3}.

## 5. Interpretation

- The correct Meyrand fluidization diagnostic is the **k⊥-integrated** Hermite flux Γ(m). Here it is
  ≈ injection → **phase mixing dominates**; the plasma is **not fluidized** at these parameters.
- The project's earlier per-k⊥ `ECHO_RATIO ≈ 0.26` (§19) measured *local* cancellation of src_m at each
  k⊥, not the *net* flux. The net k⊥-integrated flux carries ~all the injection forward, so the low
  per-k ratio does **not** indicate an echo. This report supersedes the "clearer echo" reading of §19.
- Reaching the echo/fluidized regime (⟨Γ⟩→0) would require pushing the nonlinear-to-phase-mixing rate
  ratio χ toward/above the fluidization threshold — stronger Alfvénic drive (larger ene_inj) and/or a
  longer Hermite inertial range (smaller ν) — a separate parameter study (deferred; see §6 caveat).

## 6. Caveats

1. **Steady state reached (was a caveat, now resolved).** W_free peaked at ≈17.66 (t≈19.5) and settled
   to 17.59 ± 0.06 (0.3 %) with a slightly negative last-half slope — a statistical steady state, not a
   monotone drift. Independently, the forward-cascade result is a **flux-balance** statement (injection =
   forward Hermite flux) certified by the telescoping invariant, so it never depended on the W trajectory.
2. The echo regime was **not searched** — this documents the current parameter point only.

## 7. Reproduce

```sh
# on Miyabi (myb), in /work/gr96/o07001/gallope_dev/verify_prod128
/opt/pbs/bin/qsub -v INPUT=echo128_restart.in,FRESH=0 job.pbs-echo128   # continue chain from restart/
# each segment: ~1 h wall, ~6 t-units; NetCDF preserved to echo128_seg_<jobid>.nc
```

Segments: `echo128_seg_2355202` (t0–4), `_2357102` (t4–10, flux off), `_2358620` (t10–16, flux on),
`_2358914` (t16.5–21.76, flux on, steady window). Analysis: `analyze.py`, local `make_report_fig.py`.

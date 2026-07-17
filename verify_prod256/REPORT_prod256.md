# prod256 — Driven KRMHD velocity-space cascade, resolution study (Phase 1-C-ii, Stage 2)

**Status: CLOSED (2026-07-15). The pre-registered 256²/k⁴ decisive test (§6, §8) is complete.
With the UV bottleneck removed (frac k⊥≥30 = 0.9%, cleaner than 128²/k⁴'s 5%), the Hermite flux
Γ(m) DECLINES to ⟨Γ⟩ ≈ 0.58 across the inertial band m=20–48 — the echo direction — confirming
128²/k⁴ (0.64) at doubled resolution and refuting the 256²/k⁸ flat 0.97 as bottleneck-masking.**
**VERDICT: partial velocity-space fluidization at β_i=1, χ~3–4. The stochastic echo suppresses the
forward Hermite flux 34% below the pure-forward value ε_g (0.58 vs 0.89) but does NOT zero it, so
the plasma is a mixed velocity-space cascade — echo-active but not fully fluidized at these
parameters. See §8.**

## 1. Result (結論)

The 256² balanced driven-echo run at the **k⁸** perp-dissipation operator (nexp_perp=4) is a
**forward Hermite cascade**, not the stochastic echo / velocity-space fluidization. The
k⊥-integrated Hermite free-energy flux is **flat** across the inertial band,
⟨Γ(m)⟩ ≈ **0.97 × injection** (band m = 20–48, mean 0.97, CV 1.7%), declining only at m > 64
where it overlaps the grid-scale UV dissipation range. This **reproduces the 128²/k⁸ result**
(flat ⟨Γ⟩ ≈ 0.96) and is clearly distinct from the **128²/k⁴** run (declining ⟨Γ⟩ ≈ 0.64,
echo direction). Doubling k⊥ from 128² to 256² did **not** flip the k⁸ answer toward the echo.

**Caveat that keeps the verdict open (§6):** the plan's premise was that doubling k⊥ would
lengthen the perp inertial range enough to make Γ(m) *operator-insensitive*. It did not — the
k⁸ UV bottleneck did not shrink with resolution (17% ≈ the 128²/k⁸ 15%), so the flat Γ at
256²/k⁸ is still bottleneck-contaminated exactly like 128²/k⁸. The decisive test is a
**256²/k⁴** companion (remove the bottleneck), the analog of the 128²/k⁴ declining-Γ result.

## 2. Run configuration

| item | value |
|---|---|
| input | `echo256_restart.in` (chain-continued from `echo256.in`, FRESH=0) |
| jobscript | `job.pbs-echo256` (`regular-g`, `select=64:mpiprocs=1`, `walltime=06:00:00`; routed to `medium-g`) |
| grid / moments | 256² × 64, nm = 128 |
| process grid | np = 64 → [P_fft=8, P_m=8, P_s=1] (nm_local=16, nly_local=32 = the validated 128² comm pattern, 2× perp slab) |
| g-forcing | driven=T, m=1 low-k shell kmin=(0,0,1)/kmax=(1,1,2), fix_power, **ene_inj_g = 1.0** |
| Alfvén forcing | Elsasser zppe/zmpe, ene_inj = 1.0, **xhl_inj = 0 (balanced)**, freq (0.9, −0.6), stir_seed = 7 |
| α inputs | β_i = τ = Z = 1, alpha_root = +1, v_th = 1 |
| dissipation | ⊥ hyper μ = 10 (**nexp_perp=4 → k⊥⁸**); Hermite ν = 1 (nexp_m=6 → m⁶) |
| time stepping | dt = 4e-5 (seeded at the 256² steady perp-CFL; dt=2e-4 in echo256.in decrements to this), cfl = 0.5, eSSPIFRK3 |
| diagnostics | write_intvl = 0.5, write_hermite_flux = T, save_restart_intvl = 0.5 |

**Single-variable change vs 128²/k⁸ (`echo128bal`):** nlx,nly 128→256 only. nlz=64, nm=128, μ, ν,
forcing identical (dissipation coeffs are normalized to k_max / m=nm ⇒ resolution-agnostic, they
auto-scale to the higher k⊥,max). This isolates perp resolution as the sole variable.

## 3. Diagnostic validation — the Γ_m flux is trustworthy

Same certified `Gamma_m_kint` diagnostic as prod128: a **telescoping invariant** with
Σ_m src_m(streaming) = 0 identically, so Γ at the top moment vanishes to machine precision
(verify_1Cii gates 6/6 PASS; `TELE_NL_KINT = 4.24e-16 ≪ 1e-11`, P_m-invariance `PMREL = 5.2e-16`).
The 256² binary is byte-for-byte the prod128 binary except grid size, so the flux diagnostic is
unchanged and its correctness carries over. Run-to-run determinism (L7) and P_m bit-invariance
(L8, XOR-fold) were established in Phase 1-A/1-B and were relied on here to recover cleanly across
the intermittent cuFFTMp/IB hangs (identical restart ⇒ identical output).

## 4. Results (steady window t ≈ 20–23.5, later-half time average)

**Figure:** `fig_prod256_verdict.png` (6 panels, 2×3; local, not committed — regenerate with
`make_report_fig.py`). Panels (a)–(c),(f) compare 256²/k⁸ (this run) against the two 128²
baselines; panels (d)–(e) show the 256²/k⁸ Alfvénic field spectra.

- **(a) Hermite flux — flat / forward.** ⟨Γ(m)⟩ ≈ injection across the inertial band m = 20–48
  (mean **0.97**, CV **1.7%**; Γ(4)=0.98, Γ(32)=0.99, Γ(48)=0.97), then a decline through
  m = 64 → 127 (Γ(64)=0.87, Γ(96)=0.68) that coincides with the UV pile-up in panel (c) — i.e.
  the dissipation range, **not** an inertial-range echo. The 256²/k⁸ curve overlays the 128²/k⁸
  curve and sits well above 128²/k⁴. **Not the echo** (echo ⇒ ⟨Γ⟩→0 in the band).
- **(b) Hermite spectrum.** ⟨W_m⟩ slope ≈ **−0.30** over m = 3–40 — shallow (near the m^{−1/2}
  phase-mixing reference, far from the m^{−1} fluidized reference), matching 128²/k⁸ (−0.29).
- **(c) g⊥ spectrum — UV bottleneck persists.** A clear grid-scale pile-up
  (tail(k≥30) peak/mid = **2.30**, frac k⊥≥30 = **17%**), essentially unchanged from 128²/k⁸ (15%)
  and much stronger than 128²/k⁴ (5%). Resolution did not separate the injection/dissipation scales.
- **(d) Alfvénic field spectra — E_u(k⊥), E_B(k⊥).** The perp kinetic and magnetic spectra
  (Σ_kz of `upe2_bin`, `bpe2_bin`, later-half average) show a clean **inertial range k⊥ ≈ 3–20
  with slope ≈ −1.4 (E_u) / −1.5 (E_B)** — consistent with the **−3/2** critically-balanced
  RMHD prediction (slightly shallower than the −5/3 guide drawn), rolling over into the
  dissipation range at k⊥ ≳ 40. Magnetic energy modestly dominates, **E_B/E_u = 1.23** (2.56 vs
  2.08). This confirms the Alfvénic bath driving the passive g-cascade is a healthy, resolved
  RMHD turbulence — not a numerical artifact.
- **(e) Elsasser spectra — E_{z⁺}(k⊥), E_{z⁻}(k⊥).** Both counter-propagating Elsasser fields
  (`zppe2_bin`, `zmpe2_bin`) follow the same ≈ −1.46 inertial slope. z⁺ dominates z⁻,
  **z⁺/z⁻ = 1.33** (5.31 vs 3.98), i.e. the residual cross-helicity imbalance **σ_c ≈ 0.14** in
  the verdict window — a mild imbalance consistent with the still-decaying σ_c in panel (f).
- **(f) Approach to steady state.** The Γ band settles to ≈0.96 and the power balance closes
  (injection 1.000 / dissipation **0.997**, ratio 0.997). The imbalance σ_c decays monotonically
  from ≈0.42 to ≈0.09 over t ≈ 10 → 34; the verdict window (t 20–23.5) has σ_c ≈ 0.15,
  the later confirmation window (t 31.5–34) σ_c ≈ 0.09.

Turbulence health: χ = **3.72** (strongly nonlinear), ⟨k⊥⟩_E = 7.0, u⊥,rms = 2.0, ⟨k_z⟩_E = 3.8.
The Alfvénic field spectra (d,e) add: E_B/E_u = 1.23 (mild magnetic excess), inertial slope ≈ −3/2,
σ_c ≈ 0.14 — a well-resolved, mildly imbalanced RMHD bath.

### 4b. g free-energy k-spectrum (full 3D g²(k_z, m, k⊥))

**Figure:** `fig_prod256_gspec.png` (4 panels; local, not committed — regenerate with
`make_gspec_fig.py`). The permanent `g2_bin` diagnostic is the full 3D free-energy spectrum
g²(k_z, m, k⊥); this figure projects it onto k⊥, k_z, and the joint (k⊥, m) plane.

- **(a) g perp spectrum E_g(k⊥)** (total + per-Hermite slices m = 0, 2, 8, 32, 64). Total inertial
  slope ≈ **−1.16** (k⊥ = 3–12), shallower than −5/3; every Hermite moment shares the same perp
  shape (self-similar in m) and the **same grid-scale UV bottleneck upturn** at k⊥ ≳ 60 — i.e. the
  bottleneck lives in **k⊥, at all m** (consistent with panel (c) of the verdict figure), not in
  the Hermite direction.
- **(b) g parallel spectrum E_g(k_z).** Slope ≈ **−1.47** (k_z = 1–8), flattening to a shallow
  plateau at high k_z — the parallel free energy the streaming operator deposits along B₀.
- **(c) Joint spectrum E_g(k⊥, m) = Σ_kz** (log-color map, with the energy-weighted ridge
  ⟨k⊥⟩_E(m)). Free energy is concentrated at low k⊥ / low m (forcing corner). The ridge is
  **U-shaped**: ⟨k⊥⟩_E falls from ≈15 (m=0) to a minimum ≈8 near m ≈ 45, then climbs to ≈30 at
  m ≈ 90.
- **(d) Scale–Hermite correlation ⟨k⊥⟩_E(m), ⟨k_z⟩_E(m).** For **m ≳ 45 both rise with m** — high
  Hermite moments are populated at progressively **smaller parallel *and* perp scales**. This is the
  textbook **phase-mixing / forward-cascade correlation** (m ~ k_z v_th t): the velocity-space
  cascade transports free energy to high m at small scales, corroborating the flat-Γ forward-cascade
  verdict. The m ≲ 45 branch is the injection/inertial region dominated by the outer scale.

**Note on the Hermite dissipation ν (why the W_m spectrum in verdict panel (b) bends).** The g⊥
UV bottleneck is a **k⊥**-space effect (k⁸ operator), present at every m (panel 4b-a). The Hermite
direction shows **no** analogous pile-up: W_m decreases monotonically to the m-cutoff (peak/min =
1.00; fraction m > 0.6·nm = **9.2%**, identical to 128²/k⁸'s 9.3% ⇒ resolution-invariant, not
under-dissipated). The Hermite hypercollision ν(m/nm)⁶ reaches O(0.01) at m ≈ 60 — a dissipation
range cleanly separated *above* the inertial band m = 20–48. The W_m bend is therefore **not**
ν being too small: it is the forcing peak at m = 1 (steep m = 1–4 drop) giving way to the
m^{−1/2} phase-mixing slope (m = 3–10) that flattens further (≈ m^{−0.2}, m = 10–48) — the shallow
forward-cascade signature — then the clean Hermite dissipation steepening (slope ≈ −1.4) at m > 40.
The lever for any future Hermite-bottleneck concern would be the **exponent nexp_m**, not ν; the
k⊥ operator (not ν) is what the decisive 256²/k⁴ test must retune.

## 5. Interpretation

- At the k⁸ operator, doubling k⊥ (128²→256²) leaves the Hermite flux **flat ≈ injection** — the
  same forward-cascade / phase-mixing signature as 128²/k⁸. The plasma is **not fluidized** at
  these parameters and operator.
- The high-m Γ decline is **dissipation-range**, not echo: it tracks the UV pile-up in the g⊥
  spectrum, which the k⁸ hyperviscosity confines to the highest k⊥ but does not remove.
- Because the UV bottleneck did not shrink with resolution, this run **does not settle** the
  128² operator-sensitivity ambiguity (k⁸ flat vs k⁴ declining). It only confirms the k⁸ branch.

## 6. Caveats

1. **Bottleneck-contaminated — verdict not fully closed.** The k⁸ result at 256² is still
   contaminated by a 17% UV pile-up, so "flat Γ" cannot yet be certified as a genuine forward
   cascade rather than an echo masked by the bottleneck. **Decisive next test: 256²/k⁴**
   (nexp_perp=2, μ retuned to remove the bottleneck). If 256²/k⁴ is also flat ⇒ genuine forward
   cascade (echo absent); if it still declines like 128²/k⁴ ⇒ the k⁸ flatness is bottleneck
   masking and the echo survives at high resolution.
2. **Mild residual imbalance.** σ_c ≈ 0.15 in the verdict window (128²/k⁸ baseline had ≈0.04);
   it is still decaying (≈0.09 by t 34). The balanced forcing (xhl_inj=0) drives σ_c→0 slowly.
   The Γ(m) flat-vs-declining contrast (0.97 vs 0.96 vs 0.64) is decisive independent of this.
3. **Intermittent cuFFTMp/IB hangs at 256² scale** (lesson L19): the run was chained across
   hangs by a resident hang-detecting watchdog (lesson L20). Each hang cost ≤0.5 t (dense
   checkpoints); no data corruption (deterministic restart). The run was stopped after the
   verdict (allocation freed); the verdict window is unaffected.

## 7. Reproduce

```sh
# on Miyabi (myb), in /work/gr96/o07001/gallope_dev/verify_prod256
/opt/pbs/bin/qsub -v INPUT=echo256.in,FRESH=1 job.pbs-echo256          # first launch (random IC)
/opt/pbs/bin/qsub -v INPUT=echo256_restart.in,FRESH=0 job.pbs-echo256  # continue chain from restart/
# each segment: ~5.8 h wall (self-terminate), NetCDF in runs_echo256/gallope.out.nc (reset per segment).
# Resident chain/hang watchdog: dev/watchdog_echo256.sh (setsid+nohup on the login node).
```

Baselines used in the figure: 128²/k⁸ = `echo128bal_seg_2359772.nc` (t18.5–25, balanced);
128²/k⁴ = `echo128mu_seg_2362294.nc` (t18–24.3). 256² verdict window:
`echo256_live_t23.nc` (t19–23.5), confirmed by `echo256_final_t25.nc` (t31.5–34).
Figures/analysis: `make_report_fig.py` (verdict, 6-panel) and `make_gspec_fig.py`
(g free-energy k-spectrum, 4-panel) — both local; need the `.nc` files in the scratchpad.

## 8. The 256²/k⁴ decisive test — verdict CLOSED (2026-07-15)

This is the §6 pre-registered companion: the same 256² grid but the **clean k⁴ perp operator**
(nexp_perp=2, μ=500, all k-space coeffs ×5 vs the 64² recipe) that removes the k⁸ UV bottleneck.
It settles the k⁸-vs-k⁴ operator-sensitivity ambiguity that kept the §1 verdict open.

| item | value |
|---|---|
| input | `echo256_restart.in` (chain-continued, FRESH=0; **k⁴**: nexp_perp=2, μ=500) |
| jobscript | `job.pbs-echo256-pf4` (`regular-g`, `select=32:mpiprocs=1`, walltime 06:00:00) |
| process grid | np = 32 → **[P_fft=4, P_m=8, P_s=1]** (P_fft=4 to defer the cuFFTMp/IB hang, L19/L20) |
| grid / moments | 256² × 64, nm = 128; β_i=τ=Z=1, alpha_root=+1, **v_th = √β_i = 1** (now derived, params.F90) |
| dissipation | ⊥ **k⁴** μ=500 (nexp_perp=2), EM k⁸ (perp 2e4 / par 2e3), Hermite ν=3 m⁶ (nexp_m=6) |
| steady window | quasi-steady tt ∈ [16.5, 20.0], later-half mean (8 records); run self-terminated clean at tt=20.12 |
| data | `echo256_live_tt20.nc` (scratchpad); analysis `echo_verdict.py`; fig `fig_echo256_k4_verdict_tt20.png` (local, not committed) |

**Result (3-panel figure `fig_echo256_k4_verdict_tt20.png`):**

- **(a) Hermite flux Γ(m) — DECLINING (echo direction).** Γ rises to the injection/flat-forward
  value **ε_g ≈ 0.89** at m≈2, then **declines** to **⟨Γ⟩ = 0.58** across the inertial band m=20–48
  (band CV 4.7% over the 8 records — settled, no trend), i.e. **34% suppression below the pure
  forward-cascade value**. This lands on the **k⁴ declining branch** (128²/k⁴ ⟨Γ⟩≈0.64) and far
  below the **k⁸ flat branch** (128²/k⁸ 0.96, 256²/k⁸ 0.97). The within-band decline is mild (7%,
  0.63→0.59); the decisive drop is injection→band (0.89→0.58).
- **(b) g⊥ spectrum — bottleneck-free (the crux).** The k⁴ operator leaves **no UV pile-up**:
  frac(k⊥≥30) = **0.9%** (vs 256²/k⁸ 17%, 128²/k⁴ 5%), peak(k≥30)/median(8≤k<20) = 0.23 (monotonic
  decline into the dissipation range). **So the declining Γ in (a) is a genuine measurement, not an
  echo masked by a bottleneck** — exactly the §6 discriminator.
- **(c) Hermite spectrum E_g(m).** Inertial-band slope **−0.32** (shallow, near m^{−1/2}
  phase-mixing), matching the k⁸ run's −0.30 — the spectrum shape is operator-insensitive; the
  **flux** is the fluidization discriminator, and it is not.

**Interpretation / verdict.** The §6 pre-registered logic resolves as: *256²/k⁴ still declines like
128²/k⁴ ⇒ the k⁸ flatness was bottleneck-masking and the echo survives at high resolution.* It does.
The 128²→256² trend on the clean operator is **downward** (0.64→0.58 — more suppression at higher
resolution), the echo-strengthening direction, so the suppression is physical, not a finite-size
artifact. **But Γ settles at 0.58, not 0**: the stochastic echo partially cancels the forward
phase-mixing flux (34% below ε_g) without eliminating it. **Conclusion: partial velocity-space
fluidization at β_i=1, χ~3–4 — a mixed cascade in which the echo is demonstrably active but does
not reach the fully-fluidized Γ→0 limit.** This is the Phase 1-C-ii Stage 2 echo verdict (task #43).

**Caveats (keep this a provisional-but-decisive verdict).** (1) 8 records; W_free drifts −3.3%
over the 2nd half — the *sign and branch* (declining, on the k⁴ side) are robust, the *exact value*
0.58 has ~5% statistical scatter. (2) The Hermite inertial range is short (dissipation cutoff m~55
squeezes the m=20–48 band); the *degree* of suppression, not its existence, is what a longer
m-range / stronger scale separation would refine. (3) χ~3–4 may sit near the fluidization threshold
— whether Γ→0 (full fluidization) emerges at higher χ or different β_i is the natural follow-up,
not required to close the operator-sensitivity question this run was designed to answer.

Note: the v_th=√β_i derivation (params.F90, 2026-07-15) is bit-identical here since β_i=1
(√1=1 exactly), so this verdict is unaffected by that change.

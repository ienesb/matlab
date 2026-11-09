# `eva/` — FIM / CRLB, MCRB, and Pseudo-True Analysis for OFDM vs CP-OTFS

ISAC delay (`tau`) / Doppler (`nu`) estimation bounds for **OFDM** and **CP-OTFS**
(plus **RCP-OTFS** on the FIM side, and precoded **SC/OTFS** for the paper
reproduction), with and without inter-carrier interference (ICI), single- and
multi-target. Four layers of analysis:

1. **FIM / CRLB** — `J` and `J⁻¹` for matched (genie) symbols.
2. **MCRB** — misspecified CRB `A⁻¹ B A⁻¹ (+ BB)` when the estimator uses
   *detected* symbols `Xhat ≠ X`. The mismatch is either an abstract **BSC** of
   rate `p`, or the **realistic CE + demodulation** error pattern of a full
   receiver (`mcrb_ce_ici.m`).
3. **Pseudo-true parameter** — the value `eta_0` the mismatched ML estimator
   converges to: the least-squares projection of the true mean onto the assumed
   means (`pseudo_true.tex`). Adds the parameter-bias floor
   `BB = (eta_bar − eta_0)(eta_bar − eta_0)ᵀ`.
4. **Paper reproduction** — Liu et al. 2026 (`impl/2605.14451v1`): CP-OFDM vs
   frequency-spread waveforms (SC/OTFS) ranging CRB, and its MCRB extension.

**Layout:** all reusable functions live in **`eva/utils/`**; the runnable
**scripts live in `eva/`** and start with
`addpath(fullfile(fileparts(mfilename('fullpath')),'utils'))`. PNG / `.mat`
outputs go to **`eva/results/`**; editable MATLAB **`.fig`** files go to
**`eva/figures/`**. The full catalog of channel models (which script defines
which channel, with parameters and toggles) is in **`CHANNELS.md`**.

---

## Key findings (what the recent runs established)

- **Paper Fig. 1(b) reproduced** (`compute_mcrb_mc_true.m`, `scenario='paper'`):
  under 16-PSK, CP-OFDM sits exactly on the Jensen bound; SC/OTFS are **~2 dB**
  above. OFDM's edge is the *flat* `|Q^H s|²` of an unspread waveform.
- **The waveform advantage is a small-sample / no-ICI effect.** In the realistic
  delay-Doppler + ICI setting (`paperv2`, `N·M = 2048`), the clean-CRB gap
  collapses to **~0.04 dB** and the MCRB gap to **1–5%** — self-averaging drives
  `E[J⁻¹] → (E[J])⁻¹` and ICI spreads *every* waveform. Needs CRN pairing +
  ≥5000 trials to resolve (see `use_crn`).
- **`BB` is waveform-agnostic at equal error rate** (`compute_mcrb_mc.m`
  paperv2): the precoder moves only the covariance term `A⁻¹BA⁻¹` (~1%), not the
  bias `BB` — which is set by the error *rate*, not the waveform.
- **The OFDM BER floor is channel *estimation*, not ICI** (`diag_ofdm_floor.m`):
  swapping in a genie CE drops BER ~33×; `dftdd` CE reaches that genie bound.
- **The residual OFDM-vs-OTFS BER gap is *diversity*, not Doppler**
  (`ber_vs_doppler.m`): at *zero* Doppler OTFS is still ~200× better — uncoded
  OFDM has diversity order 1; no receiver fix crosses that wall.
- **BER and MCRB are not the same currency** (`mcrb_ce_ici.m` + `plot_ber_mcrb`):
  above ~18 dB OTFS has orders-of-magnitude lower BER yet the two MCRBs are
  within **~1%** — sensing integrates the whole frame and is noise-limited once
  the symbol-error rate is small.
- **Same SER ≠ same MCRB** (`do_bsc` arm + `plot_mcrb_real_vs_bsc`): the abstract
  i.i.d. BSC understates the sensing penalty vs the real detector's *clustered*
  errors in the low-SNR / high-SER regime.
- **Delay is ~3× more error-placement-sensitive than Doppler** (2000-trial
  `mcrb_ce_ici_results_2000.mat` + `diag_error_geometry.m`): at low SNR the real
  detector's LB is ~3× the same-count BSC for delay but only ~1.4× for Doppler —
  errors sit in the delay-set frequency fades; the CE tracks the time axis out.
- **OFDM CE reaches the genie bound; the BER floor is *ICI*, not CSI**
  (`ber_ofdm_ce_vs_perfect.m`): DFT+DD CE = perfect CSI to <1%; even perfect CSI
  floors at ~4e-4 (residual ICI), and the floor **vanishes at ν=0**. AWGN
  (`ber_awgn.m`) validates OFDM = OTFS = `Q(√SNR)`.
- **Doppler is *diversity for sensing*** (`mcrb_ce_ici` `doppler_scale=0`,
  `plot_mcrb_nodopp_wf.m`): removing Doppler *raises* the delay CRB ~1.5×
  (distinct per-path Dopplers separate the multipath in the 2-D delay-Doppler
  plane). OFDM and OTFS share the same CRB/MCRB with no Doppler; only BB (BER)
  separates them.
- **ICI is negligible for the bound at nominal Doppler** (`mcrb_ce_ici`
  `use_ici=false`, `doppler_scale=1`): the no-ICI delay CRB matches the full-ICI
  one to <0.3% (νT≈0.03). Doppler *presence* moves the bound, ICI does not.
- **The waveforms split by parameter via error geometry** (`plot_mcrb_nodopp_wf`,
  delay + Doppler figures): OTFS wins **delay** (data delay-guarded from the paths
  → tiny delay BB); OFDM wins **Doppler** at low SNR (frequency-clustered errors
  spread in Doppler → tiny Doppler BB). CRB/MCRB identical; only BB differs.
- **Doppler-only CRB dual of the paper** (`run_fim_avg` `paperdopp`,
  `plot_crb_doppler.m`): CP-OFDM sits ~2.2 dB below CP-OTFS for Doppler-only
  estimation — the Jensen-bound mirror of the paper's ~2 dB delay result.
- **OTFS's BER edge is multipath diversity, and the fading model matters**
  (`ber_rayleigh_flat.m`, `ber_rayleigh_eva.m`): a *single* flat fade → OFDM =
  OTFS (no diversity to harvest); 3-tap Rayleigh → OFDM stays order-1 while OTFS
  harvests delay diversity. The usual fixed-|amp| EVA model *overstates* OTFS's
  gain vs genuine Rayleigh taps. Per-subcarrier `|H|` of the const-|amp| channel
  is Rayleigh-like (diversity order 0.99 — `plot_hist_gain_constmag.m`), which is
  why uncoded OFDM is order-1 regardless of tap count.
- **"BER ≠ sensing" is real but narrow** (`mcrb_sweep_grid.m`,
  `plot_ber_lb_combined.m`): CP-OFDM has the lower LB *despite* worse BER only for
  **Doppler at low SNR (0–7 dB)** (comparable error rates, different error
  geometry). Elsewhere OTFS's lower BER ⇒ lower LB. The OTFS **delay-LB advantage
  is its pilot guard** (a design knob): shrinking it (`ber_confirm_config` /
  `mcrb_sweep_grid` guard param) makes OFDM win both LBs — but floors OTFS's BER,
  so the trade can't be escaped.

**OTFS via (I)SFFT unitarity:** the CP-OTFS FIM/MCRB/pseudo-true all equal the
OFDM ones evaluated at `X = ISFFT(x)`, per realization, because the (I)SFFT is
unitary. So every `*_CP_OTFS_*` function is a thin wrapper that maps the
delay-Doppler grids to time-frequency and delegates to the OFDM builder. (The
*expected* CP-OTFS `B` differs from OFDM only through the bias covariance
`Cov(b)`, since the BSC acts in the delay-Doppler domain — see `cp_otfs_mcrb.tex`
§4.) RCP-OTFS breaks the equivalence and is built directly (FIM only).

---

## Parameter convention

`eta = [alphaR, alphaI, tau, nu]` per target, stacked target-by-target for the
multi-target versions (`eta = [eta_1; …; eta_P]`, a `4P×1` vector).

- **MCRB component matrices:** `B = Jtilde + b·bᵀ`, `A = −Jtilde + A_Δ`,
  `MCRB = A⁻¹ B A⁻¹`. With `Xhat = X` (no errors) `b = 0`, `A = −J`, `B = J`, so
  `MCRB → CRB`.
- **`eta_eval` (evaluation point):** the `getAB_*` realization wrappers take an
  optional trailing `eta_eval` (a `4P×1` stacked vector). `A, B` are then formed
  at `eta_eval` (e.g. the pseudo-true `eta_0`), while the true mean in the
  mismatch uses the ground-truth `(alpha, tau, nu)`. Omit it to evaluate at the
  ground truth (default, backward compatible).
- **BSC symbol model:** `Xhat` = Gray-coded QPSK after a binary symmetric channel
  of crossover `p`; `E[Xhat|X] = (1−2p)X`, `E[E* Xhat] = −2p` (`E = X − Xhat`).
  The expected functions average over i.i.d. unit-energy QPSK with this BSC.

---

## Scripts (run these)

### FIM / CRLB
| Script | What it does |
|--------|--------------|
| **`run_fim_avg.m`** | Main FIM driver. Builds an EVA or paper-mode channel, draws Monte-Carlo QPSK/QAM grids, computes the eight FIMs (CP-OFDM / RCP-OFDM / CP-OTFS / RCP-OTFS × ICI/no-ICI), extracts per-path CRBs via a Schur complement, saves `eva_fim_avg_results.mat`. Selectors `param_mode` / `channel_type` ∈ `{'eva','paper','paperv2','paperdopp'}`. `'paperv2'` = paper's 40-path delay profile plus a per-path Doppler on an M=16-slot frame, estimated jointly. **`'paperdopp'`** = the Doppler-only **dual** of `'paper'`: N=1 (no delay dimension), M=128, 40 targets with τ=0 and distinct Dopplers → Doppler CRB only. The `estimate_delay`/`estimate_doppler` flags (auto-set from the channel) pick the estimated (τ,ν) sub-block in `crb_multi`. RCP-OFDM (`getFIM_RCP_OFDM_multi`) = single frame-level CP; reduces to CP-OFDM at M=1. |
| **`plot_crlb.m`** | Loads `eva_fim_avg_results.mat`, plots CRLB-vs-SNR (trace and path-1, in `s²` or range `m²`), with continuous-time floor and uniform-spectrum reference. If the run estimated Doppler (`paperv2`), also plots the Doppler CRB (trace and path-1, `Hz²`). |
| **`plot_crb_doppler.m`** | Doppler-only (`paperdopp`) trace Doppler CRB vs SNR, CP-OFDM vs CP-OTFS (ICI and no-ICI). The FIM is linear in 1/σ², so CRB(SNR) = CRB(σ²=1)·σ². → `crb_doppler_vs_snr.png`. |

### MCRB
| Script | What it does |
|--------|--------------|
| **`compute_mcrb_mc.m`** | Monte-Carlo **pseudo-true** MCRB (`BB ≠ 0`). `scenario` switch: `'eva'` (original four configs `{OFDM,CP-OTFS}×{no-ICI,ICI}`, QPSK, physical waveforms) / `'paper'` (N=128, M=1, delay-only, 16-PSK, precoded OFDM/SC/OTFS) / `'paperv2'` (N=128, M=16, delay-Doppler + ICI). Paper scenarios: median over trials, `Ω_N` filter (`kappa_thresh`), CRN pairing (`use_crn`), pilot-only CRB, per-trial τ- and ν-traces. Warm-started varpro pseudo-true. Saves `compute_mcrb_mc_{,paper_,paperv2_}results.mat`. |
| **`compute_mcrb_mc_true.m`** | **True-parameter** MCRB (`eta_0 := eta_bar`, `BB = 0`, so `LB = A⁻¹BA⁻¹`). Same `scenario` switch as above. The paper Fig. 1(b) reproducer: 16-PSK precoded OFDM/SC/OTFS, delay (`paper`) or delay+Doppler+ICI (`paperv2`), median stats, pilot-only CRB. Saves `compute_mcrb_mc_true_{paper,paperv2}_results.mat`. |
| **`plot_mcrb.m`** | Loads `compute_mcrb_mc_results.mat` (eva) → two `tr{MCRB_τ}` figures (vs SNR at fixed `p`, vs `p` at fixed SNR), four configs each. |
| **`plot_mcrb_paper.m`** | Loads a `compute_mcrb_mc_true_*` file (`resultfile` toggle). Delay (m²) and, for `paperv2`, Doppler (Hz²) MCRB vs SNR and vs `p`, with the pilot-only line + BER-crossover printout. |
| **`plot_mcrb_paper_pt.m`** | Pseudo-true (`compute_mcrb_mc`) analogue: adds the `BB` bias floor, plus `A⁻¹BA⁻¹`-vs-`BB` decomposition figures and the Doppler figures. |
| **`plot_mcrb_compare_true_pt.m`** | Overlays the true-parameter (`A⁻¹BA⁻¹`) and pseudo-true (`+BB`) bounds on one channel — shows the error moving from the score-bias `F` (at `eta_bar`) to the parameter bias `BB` (at `eta_0`). |
| **`mcrb_avg.m`** | Older single-config driver: OFDM multi-target ICI MCRB **and** CRB, affine `G,F`, sweeping `p ∈ {0,0.1,0.5}`. Saves `results/mcrb_avg_results.mat`. |
| **`compare_ofdm_cpotfs_meanAB.m`** | Closed-form (no MC) check that `E{A}_cpotfs = E{A}_ofdm` and `E{B}_cpotfs − E{B}_ofdm = Cov(b)_cpotfs − Cov(b)_ofdm`. |

### Realistic CE + demodulation (BER → MCRB, EVA channel)
The mismatch here is the error pattern of a **full receiver** (not an abstract
BSC). OFDM: periodic-pilot LS → delay-tap (DFT) or decision-directed (`dftdd`)
CE → 1-tap MMSE. OTFS: embedded-pilot threshold CE → message-passing (MP)
detector. Pilot energy and overhead are matched between waveforms.
| Script | What it does |
|--------|--------------|
| **`ber_ce_channels.m`** | BER vs SNR (no MCRB). Configs = OFDM `{interp, dft, dftdd}` CE × OTFS-MP, on channels `{eva, paperv2comms}`. Data-referenced matched pilots; can reuse OTFS curves from a prior run. → `plot_ber_ce_channels`. |
| **`plot_ber_ce_channels.m`** | Replot of the above from the `.mat` (masks zero-BER points, draws the `1/(nTrials·bits)` measurement floor). |
| **`diag_ofdm_floor.m`** | **Why does OFDM floor?** Decomposes the floor: real CE + ICI vs genie CE + ICI vs genie CE + no-ICI. Verdict: the floor is CE (linear interp chording the phase rotation between pilots), not ICI. |
| **`ber_vs_doppler.m`** | BER vs mobility (velocity → `nu_max`) at fixed SNR, pilot lattice frozen. Isolates diversity from Doppler: at `v=0` (no ICI) OTFS is still ~200× better → the gap is diversity order (1 for uncoded OFDM), not Doppler. |
| **`ber_ofdm_ce_vs_perfect.m(SNR, nu_scale)`** | OFDM BER: 3 CE methods vs **perfect CSI**, reduced EVA, OTFS overlaid. `nu_scale=0` → `_nodopp` file. Shows DFT+DD reaches genie and the ~4e-4 floor is **ICI** (gone at ν=0). Merges/extends over its own SNR ranges. |
| **`ber_awgn.m`** | BER over a flat unit AWGN channel — OFDM (all CE) + CP-OTFS + theory `Q(√SNR)`. Validation: everything collapses onto theory. |
| **`ber_rayleigh_flat.m`** | Single-path FLAT Rayleigh BER (α~CN(0,1), no Dopp/ICI) — OFDM (perfect CSI) vs CP-OTFS (MP) vs Rayleigh theory. Control: no diversity → OFDM = OTFS. |
| **`ber_rayleigh_eva.m(fading, doppler_on)`** | Reduced-EVA 3-tap BER. `fading` = `'rayleigh'` (α_p~CN(0,pow_p)) or `'constmag'` (fixed \|amp\|, random phase); `doppler_on` = false (no Dopp/ICI) or true (~444 Hz Doppler + ICI). OFDM (perfect CSI) vs CP-OTFS (MP). |
| **`plot_ber_rayleigh_vs_constmag.m`** | Overlay of the two tap-fading models (OFDM=orange, OTFS=blue; Rayleigh solid, const-\|amp\| dashed) — shows const-\|amp\| overstates OTFS diversity. |
| **`plot_hist_gain_constmag.m`** | Histogram of per-subcarrier \|H\| for the const-\|amp\| channel vs Rayleigh, + deep-fade-tail slope (diversity order ≈ 0.99). |
| **`mcrb_sweep_grid.m(grid,SNR,nTrials,out,Emax)`** | Low-SNR LB ratio (LB_OTFS/LB_OFDM) for OFDM vs OTFS over a (N,M,Δf,`guard_extra`) grid — the "can OFDM win both LBs?" study. `guard_extra` (grid col 4) shrinks the OTFS delay guard; `Emax` truncates the DD Doppler spread (speed). Per-block checkpoint. |
| **`ber_confirm_config.m(guard_extra,Emax,nTrials,SNR)`** | BER at the B=15.4 MHz config with the same receivers as the LB — full guard (G=13, OTFS cliffs) vs shrunk (G=9, OTFS floors). Checkpointed. |
| **`plot_mcrb_sweep_bw.m` / `plot_mcrb_sweep_guard.m` / `plot_mcrb_confirm.m`** | LB ratio vs bandwidth / vs guard size / G13-vs-G9 before-after. |
| **`plot_lb_g13_absolute.m` / `plot_ber_lb_combined.m`** | Absolute delay & Doppler LB vs SNR; and the one-figure BER+LB "BER ≠ sensing" summary. |
| **`diag_crlb_ici.m`** | Genie CRLB with the ICI kernel modelled vs ignored, swept over M — ICI is second-order for delay, but the *sole* Doppler cue at M=1. |
| **`diag_error_geometry.m`** / **`diag_marginal_bsc.m`** | Where the OFDM detector's errors land (frequency vs time clustering) and whether a marginal-matched BSC reproduces the delay/Doppler bias — the Part-7 placement mechanism. |
| **`mcrb_ce_ici.m`** | Pseudo-true MCRB where the mismatch is the realistic detector's `Xhat`. `ce_mode` (`dftdd` default), matched data-referenced pilots, saves true **BER** *and* SER, MCRB + genie CRB, τ- and ν-blocks. `do_bsc` adds a **BSC-matched control arm**: same per-trial error *count*, re-scattered i.i.d. — isolates error *pattern* from *rate*. Saves `mcrb_ce_ici_results.mat`. |
| **`plot_mcrb_ce_ici.m`** | Delay & Doppler RMSE bounds vs SNR: MCRB (solid) + genie CRB (dashed). `plot_mcrb_ce_ici(matfile, include_BB)` toggles the `BB` term. |
| **`plot_mcrb_ce_ici_bb.m`** | Overlay per parameter: MCRB / without-BB (`N₀G+F`) / genie CRB — shows where `BB` bites and where the bound is noise-limited. |
| **`plot_ber_mcrb.m`** | Three stacked panels sharing the SNR axis: detection BER/SER, delay MCRB+CRB, Doppler MCRB+CRB — makes "BER ≠ MCRB" visible at a glance. |
| **`plot_mcrb_real_vs_bsc.m`** | Overlays the realistic-error MCRB (solid) vs the BSC-matched MCRB (dashed) — tests whether the abstract BSC reproduces the real bound at equal SER. |
| **`plot_bsc_vs_ce.m`** | Compares the abstract-BSC MCRB (`compute_mcrb_mc` `evace`) with the realistic-detector MCRB (`mcrb_ce_ici`) SER-matched at each SNR, for delay & Doppler, plus the detection SER/BER panel. |
| **`plot_mcrb_nodopp_compare.m`** | Per waveform, delay MCRB & LB **with ICI vs without ICI** (nominal vs `nu=0`) — the Doppler-as-sensing-diversity comparison. |
| **`plot_mcrb_nodopp_wf.m`** | **OFDM vs CP-OTFS** LB/MCRB/CRB, **delay AND Doppler** figures, for any `mcrb_ce_ici` run (regime label auto-read from the file: nominal / no-ICI / no-Doppler). Shows CRB/MCRB identical, LB split only by BB — OTFS wins delay, OFDM wins Doppler at low SNR. |

### MMLE (estimator)
| Script | What it does |
|--------|--------------|
| **`compute_mmle_mc.m`** | Monte-Carlo **misspecified MLE** for the four configs `{OFDM, CP-OTFS} × {no-ICI, ICI}` on the reduced EVA channel — the estimator counterpart of `compute_mcrb_mc.m` (**no** `A,B`/MCRB is formed here). Generates the noisy observation from the *transmitted* symbols and runs the MMLE on the *detected* reference; sweeps `p` and SNR. Records per-trial estimates and the empirical `MSE_bar` about `eta_bar`. Saves `results/compute_mmle_mc_results.mat`. |
| **`plot_mmle.m`** | Loads the above and produces two figures of `tr{MSE_τ}` (empirical MMLE delay MSE): vs SNR at fixed `p`, and vs `p` at fixed SNR — four config curves each. Mirrors `plot_mcrb.m` so the curves overlay on the pseudo-true bound. |
| **`plot_mcrb_mmle.m`** | Combined overlay of both: `tr{LB_τ}` (**dashed**, from `compute_mcrb_mc`) and `tr{MSE_τ}` (**solid**, from `compute_mmle_mc`), color = config. Two figures (vs SNR at fixed `p`, vs `p` at fixed SNR), eight curves each — shows the MMLE sitting at/above its MCRB. |

### Pseudo-true parameter
| Script | What it does |
|--------|--------------|
| **`pseudotrue_diff.m`** | Pseudo-true vs true **delay** difference `Σ_p (τ0_p − τ̄_p)²` vs BSC rate `p` at fixed SNR, four curves (OFDM / CP-OTFS × no-ICI / ICI). Saves `results/pseudotrue_diff_compare.png`. Reuses the OFDM cache when present. |
| **`pseudotrue_diff_ofdm.m`** | OFDM-only version of the above (two curves). Produces the cache `results/pseudotrue_diff_ofdm_results.mat`. |

### Verification & tests
| Script | What it checks |
|--------|----------------|
| **`verify_crb.m`** | FIM 4-way equality on reduced EVA: OFDM/OTFS, MC-average vs closed-form expectation, all four cases. |
| **`verify_mcrb.m`** | MCRB `A,B`: realization MC-average vs closed-form expectation (with `Cov(b)`), OFDM and OTFS, all four cases. |
| **`test_pseudo_true.m`** | `pseudoTrue_OFDM`: matched → exact recovery; BSC → gain → `(1−2p)ᾱ`, delay/Doppler ≈ unbiased. |
| **`test_pseudotrue_otfs.m`** | `pseudoTrue_CP_OTFS`: matched, **consistency vs `pseudoTrue_OFDM(ISFFT·)`**, gain attenuation. |
| **`test_getab_pseudotrue.m`** | `getAB_OFDM_pseudoTrue` + the eval/true split in `getAB_OFDM_core` (matched → `A=−B=−J`; backward compat). |
| **`test_getab_eta.m`** | The optional `eta_eval` on all 8 `getAB` wrappers (default = ground truth; `eta_eval` = pseudo-true matches `getAB_OFDM_pseudoTrue`; OFDM↔OTFS consistency). |
| **`check_dmu_dnu.m`** | Finite-difference check of the Doppler derivative `dμ/dν` used by the ICI FIM. |

---

## Functions (`eva/utils/`)

### FIM — realization (genie symbols)
| Function | Description |
|----------|-------------|
| `getFIM_OFDM_{,multi_}{noICI,ICI}.m` | OFDM FIM; the ICI multi variant `J = (2/σ²)·real(Gᴴ G)` is the core everything reduces to. |
| `getFIM_CP_OTFS_{,multi_}{noICI,ICI}.m` | CP-OTFS FIM — wrappers mapping `X = ISFFT(x)` to the OFDM builders. |
| `getFIM_RCP_OTFS_multi.m` (+ `_{,multi_}{ICI,noICI}` wrappers) | **RCP-OTFS** FIM, built directly in the delay-Doppler domain (single frame CP → intra-frame Doppler phase + boundary-wrap phase, useful duration `T`). Matches `rcp_otfs_fim.tex`. |

### FIM — closed-form expected (i.i.d. QPSK)
| Function | Description |
|----------|-------------|
| `getFIM_OFDM_{,multi_}{noICI,ICI}_expected.m` | Analytic `E{J}` (Dirichlet/ICI moment sums; QPSK gives `|X|²=1`). |
| `getFIM_CP_OTFS_*_expected.m` | `E{J}_CP-OTFS = E{J}_OFDM` (unitarity) — delegate to the OFDM expected FIMs. |

### MCRB — realization
| Function | Description |
|----------|-------------|
| **`getAB_OFDM_core.m`** | Core `A, B` builder (all OFDM cases via `numel(alpha)`, `use_ici`). Separates the **evaluation** parameter from the **true** parameter `eta_bar`: `Δ = μ(X, eta_bar) − μ̃(Xhat, eta)`. Returns `[A,B,Jtilde,b]`. |
| `getAB_OFDM_{,multi_}{noICI,ICI}.m` | Thin wrappers over the core; accept the optional `eta_eval`. |
| `getAB_CP_OTFS_{,multi_}{noICI,ICI}.m` | CP-OTFS wrappers: `X = ISFFT(x)` then the OFDM wrapper; pass `eta_eval` through. |
| **`getAB_OFDM_pseudoTrue.m`** | One-shot: solves the pseudo-true `eta_0`, then `A, B` at `eta_0` (true mean from `eta_bar`). |

### MCRB — closed-form expected
| Function | Description |
|----------|-------------|
| `getAB_OFDM_{,multi_}{noICI,ICI}_expected.m` | Analytic `E{A}`, `E{B} = E{Jtilde} + E{b}E{b}ᵀ + Cov(b)`. |
| `getAB_CP_OTFS_*_expected.m` | Reuse OFDM `E{A}`, `E{Jtilde}`, `E{b}`; swap `Cov(b)_OFDM → Cov(b)_CP-OTFS`. |
| `getMCRB_covb_OFDM.m`, `getMCRB_covb_CP_OTFS.m` | Bias covariance `Cov(b)` (4th-moment trace formula; `cp_otfs_mcrb.tex` §4). |

### Pseudo-true parameter
| Function | Description |
|----------|-------------|
| **`pseudoTrue_OFDM.m`** | `eta_0 = argmin_eta ‖μ(X,eta_bar) − μ̃(Xhat,eta)‖` via `lsqnonlin` (scaled τ,ν). All four cases. |
| **`pseudoTrue_CP_OTFS.m`** | DD-domain pseudo-true — `X = ISFFT(x)` then `pseudoTrue_OFDM` (exact by unitarity). |

### MMLE (estimator)
| Function | Description |
|----------|-------------|
| **`mmle_OFDM.m`** | OFDM misspecified MLE (`ofdm_mmle.tex`): the gain-concentrated subspace-projection search `θ̂ = argmax ‖P_C(θ) y‖²`, `α̂ = C†y`, with the steering matrix `C` built from the *detected* `Xhat`. CLEAN (per-target periodogram + refine) initializes a **joint** delay-Doppler refinement; gains by LS. Resolves sub-bin-spaced targets. |
| **`mmle_CP_OTFS.m`** | DD-domain MMLE (`cp_otfs_mmle.tex`): by (I)SFFT unitarity `Y^TF = ISFFT(y^DD)`, `Xhat^TF = ISFFT(xhat^DD)`, then delegate to `mmle_OFDM` (same metric in either grid). |

### Helpers
| Function | Description |
|----------|-------------|
| `get_ici_coeffs.m` | ICI kernel `I(n) = (e^{j2πνT} − 1)/(j2π(n+νT))` (L'Hôpital `=1` at `n+νT=0`). |
| `get_ici_coeffs_derivative.m`, `…_derivative2.m` | `dI/dν` and `d²I/dν²` (with singular-diagonal limits). |
| `sfft_u.m` / `isfft_u.m` | Unitary SFFT / ISFFT (TF ↔ delay-Doppler); mutual inverses. |
| `unpackEta.m` | Stacked `eta` (`4P×1`) → `(alpha, tau, nu)`; inverse of the parameter stacking. |

---

## Typical workflows

```matlab
cd eva

% FIM / CRLB
run_fim_avg;  plot_crlb

% MCRB, reduced-EVA four configs (set scenario='eva' at the top)
compute_mcrb_mc;  plot_mcrb

% Paper reproduction — CP-OFDM vs SC/OTFS ranging CRB/MCRB
%   set scenario='paper'   (delay-only) or 'paperv2' (delay+Doppler+ICI)
compute_mcrb_mc_true;  plot_mcrb_paper       % true-parameter (BB=0)
compute_mcrb_mc;       plot_mcrb_paper_pt    % pseudo-true (+BB)
plot_mcrb_compare_true_pt                    % overlay the two

% Realistic CE + demodulation
ber_ce_channels                              % BER vs SNR (OFDM CE variants vs OTFS-MP)
diag_ofdm_floor                              % why OFDM floors (CE, not ICI)
ber_vs_doppler                               % the gap is diversity, not Doppler
mcrb_ce_ici;  plot_ber_mcrb                  % BER + MCRB from realistic detection
plot_mcrb_ce_ici_bb;  plot_mcrb_real_vs_bsc  % BB split; real vs BSC-matched errors

% Pseudo-true delay difference, OFDM vs CP-OTFS, vs p
pseudotrue_diff

% Sanity checks
verify_crb;  verify_mcrb;  test_getab_eta
```

**Faithful paper Fig. 1** (SC / AFDM / OTFS / OFDM / Jensen, single-block model)
lives in **`impl/2605.14451v1/run_sim.m`** — the `compute_mcrb_mc_true`
`'paper'` scenario is the same physics ported into this repo's precoded-OFDM
FIM machinery (so it also carries the MCRB / detection-error extensions).

Evaluate `A, B` at the pseudo-true parameter:
```matlab
eta0 = pseudoTrue_OFDM(X, Xhat, deltaf, Ts, alpha, tau, nu, use_ici);
[A,B] = getAB_OFDM_multi_ICI(X, Xhat, deltaf, Ts, sigma2, alpha, tau, nu, eta0);  % at eta_0
[A,B] = getAB_OFDM_multi_ICI(X, Xhat, deltaf, Ts, sigma2, alpha, tau, nu);        % at ground truth
% CP-OTFS: eta0 = pseudoTrue_CP_OTFS(x, xhat, ...); getAB_CP_OTFS_multi_ICI(x, xhat, ..., eta0)
```

---

## Reference derivations (`Desktop/latex/*.tex`)

Function headers cite section numbers in these:
- `ofdm_io` / `ofdm_fim` / `ofdm_mcrb` — OFDM I/O, FIM, and MCRB (A, B, Cov(b)).
- `cp_otfs_io` / `cp_otfs_fim` / `cp_otfs_mcrb` — CP-OTFS counterparts (§4 = bias covariance).
- `rcp_otfs_io` / `rcp_otfs_fim` — RCP-OTFS I/O and FIM (the two extra phases).
- `pseudo_true` — pseudo-true parameter (KL = least-squares projection).
- `ofdm_mmle` / `cp_otfs_mmle` — misspecified MLE (matched-filter / subspace-projection estimator on the detected reference).
- `eva_channel` — EVA and reduced-EVA channel; `ct_crlb_benchmark` — continuous-time CRLB floor.
```

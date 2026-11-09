# Experiments Log — OFDM vs OTFS for ISAC (ranging CRB / MCRB / BER)

This document walks through every experiment performed, in the order the
investigation actually unfolded. Each entry states the **question**, the
**script + setup**, the **result**, and the **conclusion** it forced (including
where a later experiment *overturned* an earlier reading).

The arc, in one sentence: *we reproduced a paper's claim that OFDM has the lowest
ranging CRB, extended it to detection-error (MCRB) and realistic-receiver
settings, and found that OFDM's apparent disadvantages are mostly receiver
artifacts and small-sample effects — while OTFS's real advantage is diversity
and shows up in BER, not in the sensing bound.*

Reduced-EVA channel throughout the receiver work: `N=512, M=128, Δf=15 kHz`,
3 taps `[0 1 3]`, `f_c=4 GHz`, `v=120 km/h → ν_max=444 Hz`, 4-QAM.

**Terminology (fixed convention):** `MCRB = A⁻¹BA⁻¹ = N₀·G + F` (noise term +
score-bias floor); `BB = (η̄−η₀)(η̄−η₀)ᵀ` (pseudo-true parameter bias);
`LB = MCRB + BB` (the total bound). Delay/Doppler "RMSE bound" = `√(Σ_p [·]_pp)`
over the τ- or ν-block.

---

## Part 1 — Reproducing the paper's ranging CRB (Fig. 1b)

**Paper:** Liu et al. 2026, *"CP-OFDM Achieves Lower Ranging CRB Than
Frequency-Spread Waveforms in the Large-Sample Regime"* (`impl/2605.14451v1`).
Single block, `N=128`, `L=40` (later `L=10`) delay-only targets. Waveforms are a
unitary precoder `U` on the data: OFDM `U=F_N^H`, SC `U=I`, OTFS `U=F_{32}^H⊗I_4`.
The FIM is `J = (2/σ²)·Re{Hᴴ Diag(|F_N U s|²) H}`; the waveform enters *only*
through the frequency-domain power `|Q^H s|²`.

**Q: can `run_fim_avg`'s physical CP/RCP-OTFS reproduce Fig. 1(b)?**
- **No.** Its OTFS resolves delay over the `N_t` DD rows, so a `32×4` frame gives
  only 32 delay bins → 40 targets alias → rank-deficient FIM (NaN/negative CRBs).
  The paper's `P_τ` is a 128-circulant (receiver `F_128`) — 128 delay bins.
- **Fix** (`compute_mcrb_mc_true.m`, `scenario='paper'`): feed the *precoded*
  subcarrier symbols `X_eff = F_N·U·s` into the 128-length OFDM FIM
  (`getFIM_OFDM_multi_noICI`). This decouples delay resolution (always 128) from
  the precoder, exactly like the paper. `N=128, M=1`, 16-PSK, `L=10`.

**Result (16-PSK, median over trials):**

| waveform | τ-trace CRB @σ²=1 | vs OFDM |
|---|---|---|
| OFDM | 1.925e-17 s² | — (**exactly the Jensen bound**) |
| OTFS | +49.7% | +1.75 dB |
| SC | +63.2% | +2.13 dB |

**Conclusion:** reproduced. Under 16-PSK, `|X|²≡1` for OFDM is *deterministic* →
its FIM equals `E[J]` → it sits exactly on the Jensen bound; spread waveforms
fluctuate, so by Jensen's inequality their averaged CRB is ~2 dB above. Matches
the paper's "~2 dB gain."

---

## Part 2 — MCRB under detection errors (abstract BSC)

The paper is CRB (perfect data). We extended it: the sensing receiver reuses
*detected* symbols `X̂ ≠ X`, so the bound becomes the misspecified CRB
`MCRB = A⁻¹BA⁻¹` at the **true** parameter (`compute_mcrb_mc_true`, `BB=0`), or
`+BB` at the **pseudo-true** parameter (`compute_mcrb_mc`). Errors modeled as a
16-PSK **BSC** of rate `p` (4-bit Gray, `bsc_16psk`).

**Numerical hardening** (the hard part): at `L=40` the α+τ FIM is near-rank-
critical (120 params in 128 obs) and `F ~ A⁻²` is *heavy-tailed*, so the naive
mean over trials is dominated by rare near-singular draws (values swinging
1e4–1e10, non-monotonic). The recipe that fixed it:
1. `L=10` (30 params in 128 obs, well-conditioned);
2. `Ω_N` filter (`kappa_thresh=1e4`) — drop the rare near-singular trials;
3. **median over trials**, not mean — the spread-waveform floor is heavy-tailed
   even when `A` is well-conditioned (one flipped bit corrupts the whole spread
   vector). Store per-trial τ-traces and take the median in the plot.

*(An "average-matrix" bound `A_bar⁻¹ B_bar A_bar⁻¹` is smooth but WRONG — it
collapses all waveforms to Jensen at `p=0`, killing the very gap we're studying.)*

**Result** (`plot_mcrb_paper*.m`): OFDM lowest, OTFS mid, SC highest; the gap
widens with BER (more spreading = more error propagation). Delay-only, `p∈[0,0.1]`.

### 2a. Pilot-only CRB and the data-vs-pilot crossover
Added a pilot-only OFDM CRB: 20% of subcarriers as known comb pilots → **no
mismatch** → pure CRB (no floor), scales `∝σ²`. Data-aided has a floor `F` that
doesn't fall with SNR. **Crossover:** at SNR=20 dB, data-aided OFDM exceeds
pilot-only at `p ≳ 1%` → above that BER, pilots win. At SNR=0 dB data always
wins (energy-limited). The classic pilot-fallback argument.

### 2b. Decomposition `A⁻¹BA⁻¹` vs `BB` (pseudo-true only)
`compute_mcrb_mc` (pseudo-true) splits the floor into the covariance term
`A⁻¹BA⁻¹ = σ²G+F` and the parameter-bias `BB`. **Finding:** at high SNR `BB`
dominates `A⁻¹BA⁻¹` by 2–3 orders → the pseudo-true floor is the *parameter
bias*, not the score bias. (`plot_mcrb_paper_pt.m` figures 3–4.)

### 2c. True vs pseudo-true (`plot_mcrb_compare_true_pt.m`)
Same channel, both bounds overlaid:
- `A⁻¹BA⁻¹` differs ~100× between the scripts — because it's evaluated at `η̄`
  (true, large score bias `b`) vs `η₀` (pseudo-true, `b≈0` by stationarity).
- **Totals:** pseudo-true is 1.2–3.1× **higher** than true-parameter (the extra
  `BB`). The true-parameter script *understates* the ranging floor. The
  pseudo-true bound is the correct one.

---

## Part 3 — paperv2: delay + Doppler + ICI

**Q: with Doppler and ICI, does OFDM still win — for delay AND Doppler?**

**Setup** (`scenario='paperv2'`): same delays + a per-path Doppler
`ν~U[−ν_max,ν_max]`, `ν_max=0.5/(M·Ts)≈39 kHz` (half a Doppler bin). **`M=16`
required** — at `M=1` the Doppler columns are structurally null. Doppler breaks
subcarrier orthogonality → `use_ici=true`. `N·M=2048` cells.

**Result (5000 paired trials, CRN, bootstrap 95% CI):** OFDM has the lower MCRB
for **both** delay and Doppler, but the margin **collapses** from ~2 dB (60%) in
the M=1 case to **1–5%**:

| p=0.1 | SC delay | OTFS delay | SC Doppler | OTFS Doppler |
|---|---|---|---|---|
| gap vs OFDM | +4.3% | +5.3% | +4.0% | +4.0% |

At `p=0` the Doppler gap is `+1.08%`/`+0.88%` with OFDM lower in **100% of
paired trials** (its FIM is deterministic under 16-PSK). At `p>0` the per-trial
win rate is only ~53% (a coin flip) — heavy-tailed `F` dominates per realization.

**Why the collapse:** (1) *self-averaging* — with 2048 cells `J` concentrates so
`E[J⁻¹]→(E[J])⁻¹`, all waveforms hit Jensen (clean-CRB gap only +0.04 dB);
(2) *ICI mixes subcarriers for every waveform*, eroding OFDM's localization edge.

**Pseudo-true paperv2** (`compute_mcrb_mc`): confirms `BB` is **waveform-agnostic**
— the significant, tight gaps are all in `A⁻¹BA⁻¹` (~1.1–1.3%); every `BB` gap's
CI straddles zero. The waveform shapes only the *covariance*, not the *bias*.

**Methodological lesson (important):** at these ~1% margins, **CRN pairing
matters more than trial count**. A 20-trial smoke test gave +27/+40% and a
1000-trial *unpaired* run gave a SC-vs-OTFS ordering that *reversed* under
pairing — both pure noise. Use `use_crn=true` + ≥5000 trials; treat any unpaired
gap < 2% as noise.

---

## Part 4 — Realistic CE + demodulation (BER)

Switched from abstract BSC to a **full receiver** error pattern
(`ber_ce_channels.m`). OFDM: periodic-pilot LS → CE → 1-tap MMSE → hard decision.
OTFS: embedded-pilot threshold CE → 30-iter message-passing (MP) detector. Pilot
energy and overhead **matched** between waveforms (`match_pilot_spacing`).

**Q1: which has lower BER?** OTFS, decisively — and it's an **error floor**, not
an offset. OFDM floors at ~1e-2 (EVA) / 3e-3 (paperv2comms); OTFS drops to 0 by
16 dB. Gap ~10 dB at BER=1e-2, unbounded below it. (Both channels tell the same
story → it's the receiver, not the PDP.)

**Q2 (`diag_ofdm_floor.m`): why does OFDM floor — ICI or CE?** Decomposition at
35 dB:

| OFDM config | BER |
|---|---|
| real CE + ICI | 1.0e-2 |
| **genie CE + ICI** | **5.2e-4** (33× better) |
| genie CE, no ICI | 8.55e-5 (Rayleigh, div-1) |

**Verdict: the floor is channel ESTIMATION, not ICI.** Pilots sit 37 subcarriers
apart; the channel rotates ~78° between them and **linear interpolation chords
that arc** — a large signal-proportional error. ICI floors only at ~5e-4, an
order of magnitude below.

**Q3 (`ber_vs_doppler.m`): is the OTFS gap Doppler robustness?** Swept velocity
0–240 km/h at fixed SNR, pilot lattice frozen. Reducing Doppler to **zero** only
improves OFDM ~2×, and **at v=0 (zero ICI by construction) OTFS is still ~200×
better.** So the gap is **diversity**, not Doppler: uncoded OFDM is diversity
order 1 (`BER~1/(4·SNR)`); OTFS spreads over the DD grid and harvests diversity
≈ the number of resolvable paths (3) → `BER~SNR⁻³`.

### 4a. Improving OFDM's receiver (`ce_mode` in `ber_ce_channels.m`)
Chased the CE bound with better estimators (no coding, no spreading):

| OFDM config (EVA @35 dB) | BER | note |
|---|---|---|
| LS + linear interp (original) | 6.7e-3 | |
| DFT CE (delay-tap LS fit) | 5.9e-4 | 11× |
| **DFT + decision-directed** | **4.7e-4** | 14×, **reaches genie CE** |
| Rayleigh wall (div-1) | 8.6e-5 | unreachable without coding |
| OTFS | ~0 | |

Also fixed the pilot-power convention: it was `∝√N0` (noise-referenced), so
pilots shrank as SNR grew and caused a BER *upturn* past 32 dB. Made it
**data-referenced** and matched for both waveforms → the upturn vanished.

**Conclusion:** with CE fixed, OFDM sits at its diversity-1 Rayleigh wall. The
remaining ≥200× to OTFS is **diversity**, recoverable only by coding/interleaving
(or DFT-spread, i.e. SC-FDMA — which *is* a spread waveform). *Uncoded* OFDM is
genuinely beaten; but uncoded OFDM is a strawman — deployed OFDM always codes.

---

## Part 5 — MCRB from the realistic detector (`mcrb_ce_ici.m`)

Fed the realistic detector's `X̂` (now with the good `dftdd` CE and matched
pilots) into the pseudo-true MCRB. SNR 0:1:25, 50 trials, saves true **BER** and
SER. This is the payoff experiment for the original question — *"OTFS has lower
BER, so does it have lower MCRB?"*

**Result (`plot_ber_mcrb.m`, 3 stacked panels):**

| SNR | OFDM BER | OTFS BER | OFDM delay MCRB | OTFS delay MCRB |
|---|---|---|---|---|
| 10 | 3.5e-2 | 3.7e-3 | 0.726 ns | 0.227 ns |
| 15 | 1.0e-2 | 6.2e-7 | 0.136 ns | 0.115 ns |
| 20 | 3.3e-3 | **0** | 0.0655 ns | 0.0630 ns |
| 25 | 1.3e-3 | **0** | 0.0365 ns | 0.0362 ns |

**Conclusion: lower BER, essentially the SAME MCRB.** Above ~18 dB the two
waveforms are within **~1%** on delay, both sitting on their genie CRB, despite
OTFS having *zero* bit errors vs OFDM's 1e-3. The reason: the MCRB is dominated
by the noise term `N₀·G` (waveform-independent to ~1% per paperv2); the bias
`BB` only bites when SER is large. **BER and MCRB are not the same currency** —
sensing integrates the whole frame and is noise-limited, not error-limited, once
SER ≲ 1e-2.

**Two things the fine grid revealed:**
- **Doppler is inverted at low SNR:** below ~10 dB OTFS is 3–4× *worse* on
  Doppler (MP failure biases `ν` harder than MMSE). Crossover ~10 dB.
- OTFS's sensing advantage lives **only in the ~6–15 dB window** — where its `BB`
  is already ~0 while OFDM's isn't.

**⚠️ This overturns the script's original premise** ("OTFS has far lower BER →
lower MCRB"), which was written when OFDM's CE floored its SER at ~1e-2. With
the CE fixed, OFDM's SER falls far enough that `BB` stops mattering.

`plot_mcrb_ce_ici_bb.m` overlays MCRB / without-BB (`N₀G+F`) / genie CRB: all
three converge above ~18 dB, and the two waveforms' *genie CRBs lie exactly on
top of each other* — any separation of the MCRBs is mismatch, never the waveform.

---

## Part 6 — Does the abstract BSC match the real detector? (`do_bsc` arm)

**Q: is the MCRB a function of `(SNR, error rate)` alone?** i.e. would
`compute_mcrb_mc` (i.i.d. BSC) equal `mcrb_ce_ici` (real detector) at matched
SER? Added a **BSC-matched control arm**: per trial, take the real detector's
error *count* and re-scatter that many errors i.i.d. uniformly (`bsc_match`) —
same channel, data, SNR, SER; **only the placement differs**.

**Result (`plot_mcrb_real_vs_bsc.m`, SER matched per trial):**

| SNR | OFDM real/BSC | regime |
|---|---|---|
| 0–8 dB | **2.5–3.7×** | real ≫ BSC (clustered errors much worse) |
| ~9–12 dB | ~1.0 | crossover |
| 13–25 dB | **0.6–0.8×** | real < BSC |

**Conclusion: same SER does NOT give the same MCRB — off by up to ~3.7×, with a
sign change.** The two regimes have *different* mechanisms (both hinge on the
mismatch being `Δ = E⊙H` = symbol error × channel; note detectors are HARD, so
this is pure error-pattern, not soft/hard):
- **High SNR (real < BSC, ~0.65×):** the 1-tap MMSE detector errs where its
  post-eq SNR `|H|²/N0` is lowest — i.e. at **deep fades**. Verified: at 20 dB
  the mean `|H|²` at real error cells is 0.023 vs 0.999 overall (**43× below
  average**). Since `Δ = E⊙H`, those errors are multiplied by a near-null `H` →
  contribute almost nothing to the sensing bias. The BSC scatters the same
  errors uniformly (`|H|²≈1`) → each biases ~43× more → `BB_real ≪ BB_bsc`.
- **Low SNR (real ≫ BSC, 2.5–3.7×): a NONLINEAR effect, not yet fully pinned
  down.** The *first-order* bias `F = A⁻¹bbᵀA⁻¹` at `η̄` says real ≪ BSC at ALL
  SNR (`F_real/F_bsc ≈ 0.01` at 6 dB — deep-fade + low `|H|` always shrinks the
  linear bias), so the observed real>BSC is *beyond first order*: at SER~0.3 the
  LS projection `η₀` drifts far and the Taylor `η₀−η̄≈A⁻¹b` breaks. A "coherent
  errors" hypothesis was **ruled out** (`coh_real ≈ coh_bsc ≈ 1`, both
  random-walk). Leading (unverified) hypothesis: the noisy channel estimate `Ĥ`
  is delay-biased, real errors are demodulated *through* it, so re-fitting `η₀`
  inherits that bias (a structured pull a big-perturbation LS amplifies); BSC
  errors are channel-agnostic. Test: genie-`Ĥ` detection at low SNR.
- **OTFS geometry differs.** OTFS errors are NOT at deep fades (`|H|²`@err ≈ 1,
  vs OFDM's 0.02–0.32) — MP is a *joint* DD detector, no per-subcarrier
  equalization, so it has no fail-at-the-null discount. The high-SNR
  deep-fade mechanism is **OFDM-specific**. (`diag_bias_pattern.m`.)

**The abstract BSC understates the sensing penalty at low SNR / high SER**
(exactly where OTFS's lower BER matters) and **overstates it at high SNR**. So
`compute_mcrb_mc` (BSC) and `mcrb_ce_ici` (real detector) will *not* agree at
matched SER — the error *placement relative to the channel* is a first-order
effect the BSC discards. *(Diagnostic: `mean |H|² at error cells` check above.)*

---

## Part 7 — Delay vs Doppler: error *placement*, at 2000 trials

**Q: which parameter is more sensitive to *where* the detector's errors land?**
The 50-trial `do_bsc` arm hinted delay ≫ Doppler; confirmed at **2000 trials**
(`mcrb_ce_ici_results_2000.mat`, ~36 h) and cross-checked against the abstract
BSC (`compute_mcrb_mc` `scenario='evace'`, `plot_bsc_vs_ce.m`).

**Result — OFDM real-detector LB vs same-count BSC (same channel, only placement
differs):**

| SNR | delay real/BSC | Doppler real/BSC |
|---|---|---|
| 0 dB | **3.3×** | 1.4× |
| 6 dB | **2.6×** | 1.2× |
| ≥12 dB | ~1.0 | ~0.9 |

**Conclusion: at low SNR the real detector's errors bias the DELAY ~3× more than
random errors of the same rate, but the DOPPLER only ~1.4×.** Mechanism (verified
`diag_error_geometry.m`, `diag_marginal_bsc.m`):
- Real errors sit in the channel's **deep fades** (`|H|²@err ≈ 0.3–0.5×` avg at
  low SNR), and are more clustered in **frequency** than time.
- Delay is read from the phase across **subcarriers**; the fade frequencies are
  set by the **delays** → error placement is delay-signature-correlated →
  strongly biases delay, unlike uniform BSC.
- Doppler is read across **symbols**; the decision-directed CE tracks the (slow)
  time variation, so errors are ~uniform in time → BSC-like → Doppler barely
  differs. Above ~12 dB both ratios drop below 1 (the deep-fade nulling of Part 6).

---

## Part 8 — Channel estimation vs perfect CSI, and the ICI BER floor

**Q: how far is each OFDM CE method from a genie receiver, and what causes the
floor?** (`ber_ofdm_ce_vs_perfect.m`, reduced EVA, 500 trials, SNR −10:50; OTFS
overlaid from the existing run.) Perfect-CSI arm = true diagonal `H` in the same
1-tap MMSE.

| OFDM CE (nominal Doppler) | BER floor (35–50 dB) |
|---|---|
| **Perfect CSI ≈ DFT+DD** | ~4e-4 (coincident) |
| DFT CE | ~5e-4 |
| LS + interp | ~6.4e-3 (~15× worse) |

- **DFT+DD reaches the genie bound** (within 1% above ~10 dB) — it fits the 4
  delay taps from all N subcarriers, so estimation error is negligible. LS+interp
  keeps a permanent ~15× CE-error floor (chords the phase rotation between pilots).
- **Even perfect CSI floors at ~4e-4.** Not a CSI problem — the 1-tap equalizer
  treats **residual ICI** (off-diagonal, from Doppler) as noise; it's SNR-
  independent → the curve flattens. On deep-fade subcarriers the ICI leakage from
  strong neighbours dominates the faded signal, setting the floor height.

**Confirmed by `nu=0` (`ber_ofdm_ce_vs_perfect(-10:50, 0)`):** with no Doppler the
floor **vanishes** — perfect-CSI/DFT+DD descend monotonically to ~1e-6 at 50 dB
(slope −1, order-1 Rayleigh, `BER≈1/(4·SNR)`). So the floor is **ICI**; the
residual ~10× gap to OTFS is **diversity** (survives at `nu=0`). LS+interp also
stops flooring at `nu=0` — its floor too was Doppler-driven.

**AWGN validation (`ber_awgn.m`, flat unit channel, 500 trials):** OFDM (all CE),
CP-OTFS, and the theory `Q(√SNR)` **all collapse onto one curve** (e.g. 3.5e-5 @
12 dB) — a unitary spread over flat AWGN can't change the statistics. Only
LS+interp sits ~2.5 dB worse (noisy per-pilot interpolation). Confirms both
simulation chains are correct. *(Fixed a latent `l_tau=0` divide-by-zero in
`match_pilot_spacing` that this flat channel exposed.)*

---

## Part 9 — Clean CRLB: does *modelling* ICI help? (`diag_crlb_ici.m`)

**Q: with genie symbols, does the ICI kernel raise or lower the delay/Doppler
CRLB?** Toggle `use_ici` on the SAME nominal-Doppler channel, sweep the time
aperture M.

| M | delay CRLB ratio (ici/noici) | Doppler CRLB ratio |
|---|---|---|
| 128 | 1.00 | 1.00 |
| 8 | 1.00–1.01 | 0.99 |
| 1 | 1.00 | **∞ → finite (ICI *enables* it)** |

**Conclusion: modelling ICI is second-order for the CRLB at large apertures.** It
becomes decisive only when the *conjugate* aperture is short: at **M=1 Doppler is
observable ONLY through ICI** (the inter-symbol phase ramp is gone). Delay is
essentially unaffected (≤1%). This is a *different knob* from Part 10 — here the
channel's Doppler is fixed and only the receiver's ICI model changes.

---

## Part 10 — No Doppler: Doppler is *diversity for sensing* (`mcrb_ce_ici` `doppler_scale=0`)

**Q: how do MCRB and LB change with no ICI?** Re-ran `mcrb_ce_ici` at `nu=0`,
2000 trials (`mcrb_ce_ici_nodopp_2000_results.mat`; OFDM 0:25, CP-OTFS 0:15,
per-block checkpointing). Compared to the nominal (ICI) run.

**Surprise: removing Doppler makes the delay bound WORSE, not better.** The genie
CRB itself is **~1.5× higher** without Doppler, by a *constant* factor across SNR
(κ ~ 10² at ≥6 dB, so not numerical):

| SNR | CRB (ICI) | CRB (no ICI) | ratio |
|---|---|---|---|
| 6 dB | 0.324 ns | 0.488 ns | 1.50 |
| 25 dB | 0.036 ns | 0.055 ns | 1.53 |

**Mechanism — path separability, not interference.** With nominal Doppler the 3
paths have *distinct* Dopplers → separated in the **2-D delay-Doppler plane** →
better-conditioned FIM → lower delay CRB. At `nu=0` they collapse onto the delay
axis (taps 0,1,3 — the 0/1 pair barely one bin apart) → more cross-coupling →
~1.5× higher delay bound. Doppler acts as a **second resolution dimension**: it
*hurts comms* (the ICI floor of Part 8) but *helps sensing*. (No contradiction
with Part 9 — that toggled the ICI *model* on a fixed channel; this changes the
channel's Doppler *content*.) The LB gap is even larger at low SNR because poor
geometry also amplifies the detection-error bias BB.

**OFDM vs CP-OTFS, no Doppler (`plot_mcrb_nodopp_wf.m`):**
- **CRB and MCRB are essentially identical** for the two waveforms (same sensing
  model, matched energy) — no intrinsic delay-information advantage.
- **LB differs only through BB:** OTFS's far lower BER → smaller BB → ~5× lower LB
  at 0 dB, converging with OFDM by ~12 dB where both are error-free. Same
  project-wide theme (Part 5): the waveform sets BB via its error rate, never the
  fundamental bound.

---

## Part 11 — No ICI but Doppler kept: ICI is negligible for the bound

**Q: with Doppler present, how much does ICI itself change the sensing bound?**
Ran `mcrb_ce_ici` with `use_ici=false` (no-ICI sensing model + no-ICI OFDM channel)
but `doppler_scale=1` (ν_max≈444 Hz kept), 1000 trials
(`mcrb_ce_ici_noici_1000t_results.mat`).  OTFS keeps its native DD channel;
only its MCRB *scoring* is no-ICI.

**Result — OFDM delay CRB [ns] across all three regimes:**

| SNR | full ICI | **no ICI (Doppler kept)** | no Doppler |
|---|---|---|---|
| 0 dB | 0.6445 | **0.6436** | 0.9826 |
| 12 dB | 0.1623 | **0.1630** | 0.2458 |
| 25 dB | 0.0360 | **0.0359** | 0.0551 |

**Conclusion: turning ICI off changes the bound by <0.3%.** At ν_max≈444 Hz,
νT≈0.03 (~−25 dB ICI), so the no-ICI and full-ICI curves are on top of each other
(confirms Part 9). What moves the bound is the **presence of Doppler**, not ICI:
both Doppler cases sit ~1.5× below the no-Doppler case (Part 10 path-separation).

## Part 12 — Delay vs Doppler: OFDM beats OTFS on Doppler, not delay

**Q (same no-ICI/full-ICI runs, both figures):** why does OFDM have the *lower*
LB for **Doppler** at low SNR but the *higher* LB for **delay**? Since CRB/MCRB
are identical between waveforms, the whole difference is **BB** (detection-error
bias), and the two detectors place errors differently in the delay-Doppler plane:

- **OTFS (DD detection, guard in DELAY):** `guard_rows = 1:(2·l_tau+1)` puts the
  paths+pilot inside the delay guard, so **all data (and its errors) sit at delay
  bins ≥7, disjoint from the paths (taps 0,1,3)** → errors can't pull the path
  delays → **low delay BB (OTFS wins delay)**. But there is **no Doppler guard**,
  so DD-localized errors concentrate at Doppler bins overlapping the paths →
  coherent ν pull → **high Doppler BB (OTFS loses Doppler)**.
- **OFDM (per-subcarrier MMSE):** errors cluster at **deep-fade subcarriers**
  (delay-set, Part 7) → strong delay pull → **high delay BB (OFDM loses delay)**;
  those same frequency-localized errors are **spread across Doppler** → diffuse ν
  pull → **low Doppler BB (OFDM wins Doppler)**.

|  | delay BB | Doppler BB |
|---|---|---|
| OFDM (frequency-clustered errors) | high (loses) | low (**wins**) |
| OTFS (delay-guarded, Doppler-exposed) | low (**wins**) | high (loses) |

Each waveform's error geometry hurts a *different* parameter → OFDM beats OTFS on
Doppler but not delay, only at low SNR (both → identical CRB above ~10–12 dB).
Identical for ICI/no-ICI because none of these mechanisms involve OFDM's ICI
(delay-fade clustering = frequency selectivity; OTFS guard = structural).
Delay side verified (Part 7 + the guard is in code); the OTFS Doppler-concentration
side is the reasoned complement (an OTFS DD error-geometry diagnostic would confirm).

## Part 13 — Doppler-only CRB: the dual of the paper (`run_fim_avg` `paperdopp`)

**Q: does OFDM's ranging-CRB advantage carry over to Doppler-only estimation?**
Added a `paperdopp` mode to `run_fim_avg`: the **dual of `paper`** — single
subcarrier (N=1, no delay dimension), M=128 slots, 40 targets with **τ=0 and
distinct Dopplers** filling `[0,1/T)`.  New `estimate_delay` flag + a Doppler-only
branch in `crb_multi` (Schur over the ν block, τ known).

**Result (trace Doppler CRB, σ²=1, 16-PSK):** CP-OFDM **2.2 dB below** CP-OTFS —
the exact dual of the paper's ~2 dB delay result, same Jensen mechanism (16-PSK →
OFDM power deterministic → sits on Jensen; OTFS fluctuates → above). Holds in both
no-ICI (6.77e3 vs 1.12e4 Hz²) and ICI (4.48e5 vs 7.51e5) models. The ICI pair is
~66× higher because filling the Doppler window pushes top targets to νT≈1 (intra-
symbol smearing) — a Doppler-specific effect delay doesn't have, but the ~2.2 dB
OFDM/OTFS ratio is preserved. (`plot_crb_doppler.m` → `crb_doppler_vs_snr.png`.)

---

## Part 14 — Rayleigh vs constant-magnitude taps: where OTFS's diversity comes from

**Q: a controlled BER study — genuine Rayleigh taps vs the usual EVA
fixed-amplitude/random-phase model, and single-path vs multipath.** All no
Doppler / no ICI, OFDM = perfect CSI + per-subcarrier ML, OTFS = pilot CE + MP.
(`ber_rayleigh_flat.m`, `ber_rayleigh_eva.m` with `fading` ∈ {rayleigh, constmag}.)

**14a. Single-path FLAT Rayleigh (`ber_rayleigh_flat.m`):** OFDM and OTFS land on
the **same** order-1 Rayleigh curve `½(1−√(γ/(1+γ)))`. A single flat fade gives no
diversity to spread over, and OTFS's ISFFT is unitary → identical per-symbol SNR.
**This is the control proving OTFS's BER edge comes entirely from *multiple
resolvable paths*.**

**14b. 3-tap Rayleigh (`ber_rayleigh_eva('rayleigh')`):** **OFDM stays on the
order-1 curve** even with 3 taps (each subcarrier is a sum of Gaussians → still
Rayleigh → diversity 1; tap count is irrelevant to the slope). **CP-OTFS drops
steeply below it** (~40× lower by 20 dB) — it harvests the 3-path delay diversity.

**14c. Rayleigh vs constant-magnitude (`plot_ber_rayleigh_vs_constmag.m`):**
- **OFDM ~model-independent:** const-mag only ~1.5–2× lower (a *coding gain*, not
  a diversity gain); both order-1. Confirmed by `plot_hist_gain_constmag.m`: the
  per-subcarrier `|H[n]|` of the const-mag channel is **Rayleigh-like with fitted
  diversity order 0.99** (finite deep-fade density at 0), just slightly fewer deep
  fades (bounded support ≤ Σ|amp| ≈ 1.73).
- **OTFS strongly model-dependent:** const-mag OTFS plunges to ~5e-7 while
  Rayleigh OTFS **flattens ~5e-5**. With fixed-magnitude taps the paths *never*
  fade, so OTFS gets clean full diversity; with Rayleigh taps all taps occasionally
  fade together → the threshold-CE + MP receiver fails → a tail floor (partly MC
  undersampling at 1000 trials). **The standard const-magnitude EVA model
  overstates OTFS's diversity gain.**

**Receiver asymmetry (matters for interpretation):** OFDM was given **perfect
instantaneous CSI** (exact realized `H` each frame) + per-subcarrier 1-tap MMSE,
which for 4-QAM / no-ICI **equals the ML detector** (MMSE-slice ≡ ZF-slice ≡ ML,
since a positive scaling can't change a QPSK quadrant, and the diagonal channel
makes per-subcarrier jointly optimal). OTFS got **pilot-based CE + suboptimal MP**.
So the comparison *handicaps* OTFS — the Rayleigh-OTFS floor is a CE/MP limit, not
fundamental; genie-CSI OTFS would keep its steep slope.

Toggle `doppler_on` added to `ber_rayleigh_eva.m` (no Doppler ⇒ no ICI; Doppler ⇒
ICI on) to extend this study to the mobile case.

---

## Part 15 — Can CP-OFDM beat CP-OTFS on BOTH delay and Doppler LB at low SNR?

**Q (the payoff question): at low SNR CP-OFDM already has the lower *Doppler* LB
despite worse BER (Part 5/12); can a design regime give it the lower LB on delay
too?** Since the low-SNR LB gap is entirely **BB (error geometry)** — CRB/MCRB are
identical — it is a *design-dependent*, not fundamental, quantity, so this is a
legitimate question. Explored with `mcrb_sweep_grid.m` (low-SNR LB for OFDM &
OTFS over a grid; ratio r = LB_OTFS/LB_OFDM, r>1 ⇒ OFDM wins).

**15a. Grid sweep (N, M, Δf), SNR 0:2:6 (`mcrb_sweep_grid`, `plot_mcrb_sweep_bw`).**
- **M and the N/Δf split don't matter** for delay; only **bandwidth B = N·Δf** does.
- Doppler always OFDM's (r≈2–5). Delay ratio climbs from ~0.7 (B=7.68 MHz) toward
  1 as B grows — but **non-monotone**: it peaks at a *tie* (~1.0) at **B≈15.4 MHz
  (Δf=30 kHz, taps [0,2,6])** then falls again (the delay guard `l_tau` grows with
  B, so beyond the sweet spot OTFS re-wins delay). Bandwidth tuning ⇒ delay tie, not
  a win.

**15b. OTFS delay guard is the lever (`mcrb_sweep_grid` guard knob, `plot_mcrb_sweep_guard`).**
Parameterised the guard `G = l_tau+1+guard_extra` (constraint `G−1 ≥ l_tau` so the
largest delay always fits; `match_pilot_spacing` keeps the OFDM/OTFS data:overhead
ratio matched as G shrinks). At B=15.4 MHz, **shrinking the guard** (G=13→9) lifts
the delay ratio monotonically above 1 while Doppler stays >1 → **OFDM wins BOTH at
SNR 2–4 dB** (confirmed 800 trials, `plot_mcrb_confirm`: G=9 delay r=1.01–1.11,
Doppler 4.8–5.4). *(G=8 is a pilot-lattice discretization outlier — 4 time pilots
break OFDM's CE — excluded.)*

**15c. The honest arc, and the trade (0–25 dB LB + BER companion).**
Extending the G=13 (full-guard) LB to **0–25 dB** (`mcrb_confirm_g13_full`,
`plot_lb_g13_absolute`) overturns the optimistic reading:
- **OFDM's only sensing win is Doppler, and only at 0–7 dB.** Delay is a tie there.
- **Above ~8 dB OTFS wins BOTH** and increasingly so: past ~15 dB OTFS is error-free
  (LB = genie CRB) while OFDM keeps residual errors (SER ~1e-3 at 24 dB → a BB
  **floor**). OFDM's delay LB floors ~0.05 ns and Doppler ~0.7 Hz; OTFS's keep
  dropping. So **the ratios do NOT converge to 1** — OTFS's lower BER *does* buy
  lower LB across most of the range.
- The BER companion at G=13 (`ber_confirm_config`, `Emax=8` DD-truncation +
  checkpointing): **OTFS cliffs to 0 by ~15 dB**, OFDM decays order-1 — OTFS wins
  comms decisively.
- **Trade:** the *decisive* "OFDM wins both LBs" only appears at the **shrunk guard
  G=9**, which **floors OTFS's BER** (its CE is starved) — so there OFDM has the
  *lower* BER too, and the "despite higher BER" premise collapses. You cannot have
  OTFS at full BER strength *and* OFDM winning both sensing LBs.

**Defensible claim:** *At low SNR (0–7 dB), CP-OFDM attains a lower Doppler-
estimation bound than CP-OTFS despite CP-OTFS's far lower BER* — because there the
two have comparable error *rates* but different error *geometry* (Part 12).
Everywhere else, lower BER ⇒ lower sensing bound (OTFS). The combined one-figure
view is `plot_ber_lb_combined` (`ber_lb_combined_g13.fig`).

---

## Overall conclusions

1. **OFDM's ranging-CRB advantage is real but fragile** — it's a small-sample,
   no-ICI, unspread-waveform effect (~2 dB). In a realistic delay-Doppler + ICI
   frame it shrinks to ~1%.
2. **Waveform choice barely affects the sensing bound.** The MCRB is
   noise-limited (`N₀·G`, waveform-independent to ~1%); the mismatch bias `BB`
   is set by the *error rate*, not the precoder.
3. **OTFS's genuine advantage is BER (diversity), not sensing.** It wins BER by
   orders of magnitude, but that buys ~1% of ranging accuracy once SER is small.
   Its BER edge is diversity order, not Doppler robustness.
4. **Receiver artifacts dominated the original picture.** OFDM's BER floor was
   channel estimation (linear interp), fixable to the genie bound; its MCRB
   inflation came from that same CE. Fix the receiver and most of the gap closes.
5. **BER and MCRB are different currencies**, and **the abstract BSC is not a
   faithful stand-in** for a real detector's clustered errors (off by up to
   ~4×, and delay ~3× more placement-sensitive than Doppler — Part 7). Both
   matter for how results should be cited.
6. **Doppler is diversity for *sensing*** (Part 10): distinct per-path Dopplers
   separate the multipath in the 2-D delay-Doppler plane, so removing Doppler
   *raises* the delay CRB ~1.5×. The same Doppler that floors the *comms* BER via
   ICI (Part 8) *improves* the ranging bound. **ICI itself is negligible** for the
   bound at nominal Doppler (Part 11, <0.3%); it's Doppler *presence* that counts.
   Waveform still doesn't matter for the fundamental bound — only BER (hence BB)
   separates OFDM and OTFS.
7. **The waveforms split by parameter via error *geometry*** (Part 12): OTFS wins
   delay (its data is delay-guarded away from the paths → tiny delay BB); OFDM
   wins Doppler at low SNR (its frequency-clustered errors spread in Doppler →
   tiny Doppler BB). And the paper's ~2 dB OFDM ranging edge has a Doppler dual
   (Part 13): OFDM's Jensen advantage is ~2 dB for Doppler-only estimation too.
8. **OTFS's BER edge is multipath diversity, and the channel model matters**
   (Part 14): on a single flat fade OFDM = OTFS; uncoded OFDM is order-1 regardless
   of tap count (each subcarrier is Rayleigh, order 0.99 confirmed); the
   fixed-amplitude EVA model *overstates* OTFS's diversity vs genuine Rayleigh
   taps. And these BER comparisons give OFDM the fairer receiver (perfect CSI + ML
   vs OTFS's pilot-CE + MP).
9. **"BER ≠ sensing" is real but narrow** (Part 15): CP-OFDM's lower-LB-despite-
   worse-BER holds only for **Doppler at low SNR (0–7 dB)**, where the two have
   comparable error *rates* but different error *geometry*. Everywhere else OTFS's
   lower BER *does* give it the lower LB (past ~15 dB OFDM's residual errors floor
   its bound while OTFS = genie CRB). A *decisive* "OFDM wins both LBs" needs the
   OTFS delay guard shrunk — which floors OTFS's BER, so it's no longer "despite
   higher BER." The delay-LB gap is set by the OTFS **guard** (design), not by any
   fundamental waveform property.

**Caveats carried throughout:** `BB`'s mean is heavy-tailed → use the median and
CRN pairing; unpaired gaps < 2% are noise; low-SNR points at 50 trials are the
least reliable. All comparisons are **uncoded** — coding would recover OFDM's
frequency diversity and change the BER story (not the MCRB one).

---

## Scripts and outputs produced

| Script | Output(s) |
|---|---|
| `compute_mcrb_mc_true.m` (paper/paperv2) | `compute_mcrb_mc_true_{paper,paperv2}_results.mat` |
| `compute_mcrb_mc.m` (paper/paperv2) | `compute_mcrb_mc_{paper,paperv2}_results.mat` |
| `plot_mcrb_paper.m` / `_pt.m` | `mcrb_paper_*` / `mcrb_paper_pt_*` PNGs (+ Doppler for paperv2) |
| `plot_mcrb_compare_true_pt.m` | `mcrb_compare_true_pt_{vs_p,vs_snr}.png` |
| `ber_ce_channels.m` + `plot_ber_ce_channels.m` | `ber_ce_channels_results.mat`, `ber_ce_channels.png` |
| `diag_ofdm_floor.m` | console (floor decomposition) |
| `ber_vs_doppler.m` | `ber_vs_doppler_results.mat`, `ber_vs_doppler.png` |
| `mcrb_ce_ici.m` (`ce_mode`, `do_bsc`, `doppler_scale`, `ofdm_only`) | `mcrb_ce_ici_results{,_2000,_nodopp_2000}.mat` (+ `CKPT_*`) |
| `plot_mcrb_ce_ici.m` / `_bb.m` | `mcrb_ce_ici_{delay,doppler,bbsplit_*}.png` |
| `plot_ber_mcrb.m` | `mcrb_ce_ici_ber_mcrb.png` |
| `plot_mcrb_real_vs_bsc.m` | `mcrb_ce_ici_real_vs_bsc.png` |
| `plot_bsc_vs_ce.m` | `mcrb_ce_ici_bsc_vs_ce_{delay,doppler,ser}.png` (evace BSC vs CE) |
| `diag_error_geometry.m`, `diag_marginal_bsc.m` | console (Part 7: error placement) |
| `ber_ofdm_ce_vs_perfect.m(SNR, nu_scale)` | `ber_ofdm_ce_vs_perfect{,_nodopp}.png/.mat` |
| `ber_awgn.m` | `ber_awgn.png`, `ber_awgn_results.mat` (Part 8 AWGN validation) |
| `diag_crlb_ici.m` | console (Part 9: ICI CRLB vs aperture) |
| `plot_mcrb_nodopp_compare.m` / `plot_mcrb_nodopp_wf.m` | `mcrb_nodopp_compare_{ofdm,cp-otfs}.png` / `..._ofdm_vs_otfs_{delay,doppler}.png` (any regime) |
| `mcrb_ce_ici.m` (`use_ici=false`, doppler on) | `mcrb_ce_ici_noici_1000t_results.mat` (Part 11–12) |
| `run_fim_avg.m` (`paperdopp`) + `plot_crb_doppler.m` | `eva_fim_avg_results.mat` (Doppler-only) / `crb_doppler_vs_snr.png` (Part 13) |
| `ber_rayleigh_flat.m` | `ber_rayleigh_flat.png/.mat` (Part 14a: single flat fade, OFDM=OTFS) |
| `ber_rayleigh_eva.m(fading, doppler_on)` | `ber_{rayleigh,constmag}_eva{,_dopp}.png/.mat` (Part 14b–c) |
| `plot_ber_rayleigh_vs_constmag.m` | `ber_rayleigh_vs_constmag.png/.fig` (Part 14c overlay) |
| `plot_hist_gain_constmag.m` | `hist_gain_constmag.png/.fig` (per-subcarrier gain ~ Rayleigh, order 0.99) |
| `mcrb_sweep_grid.m(grid, SNR, nTrials, out, Emax)` | `mcrb_sweep_{grid,bw,guard,confirm}_results.mat`, `mcrb_confirm_g13_{ext,full}_results.mat` (Part 15; guard knob = grid col 4) |
| `plot_mcrb_sweep_bw.m` / `plot_mcrb_sweep_guard.m` | `mcrb_sweep_{bw,guard}.png/.fig` (LB ratio vs B / vs guard) |
| `plot_mcrb_confirm.m` / `plot_lb_g13_absolute.m` | `mcrb_confirm_guard.fig` (G13 vs G9) / `lb_g13_absolute.fig` (abs delay & Doppler LB) |
| `ber_confirm_config.m(guard_extra, Emax, nTrials, SNR)` | `ber_confirm_config_g{13,9}_results.mat` (BER at the config; G=13 OTFS cliffs, G=9 floors) |
| `plot_ber_lb_combined.m` | `ber_lb_combined_g13.fig` (BER + LB-ratio, one figure — the whole thesis) |

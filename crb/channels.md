# Channels used in `matlab/crb/`

Two primary channels are defined by the two papers in `papers/`. Everything in
`matlab/crb/` (CRB, BER, MCRB/LB) runs on one of them or on a documented variant.

| Name | Source | What it is |
|---|---|---|
| **EVA** | P. Raviteja, K. T. Phan, Y. Hong, "Embedded Pilot-Aided Channel Estimation for OTFS in Delay–Doppler Channels," IEEE TVT 68(5), 2019 (`papers/Embedded_Pilot-Aided_…pdf`), Sec. IV; tap table from 3GPP TS 36.104 (its ref. [19]) | 9-path Extended Vehicular A multipath channel with one Jakes Doppler per path, 120 km/h at 4 GHz. Mobile comms channel; delay **and** Doppler. |
| **MT40** | F. Liu, Y. Xiong, Y.-F. Liu, J. Yang, C. Masouros, S. Jin, "CP-OFDM Achieves Lower Ranging CRB Than Frequency-Spread Waveforms in the Large-Sample Regime," arXiv:2605.14451 (`papers/2605.14451v1.pdf`), Sec. V-C | 40 static point targets with distinct random delays in one 128-sample block at 160 MHz; **delay only, no Doppler**. Sensing (ranging) channel. The paper gives it no name; "MT40" = *Multi-Target 40*, after the paper's own phrase "multi-target ranging". |

A "channel" here is the set of paths `(α_p, τ_p, ν_p)` plus the rule for drawing them per trial.
Receiver knobs are **not** channels and live in `README.md`: the BSC error model, the OFDM
channel-estimation mode, the MP damping, the OTFS guard size, the DD spread truncation `E`, the
pilot budget.

## 1. Conventions and the translation table

The code follows the `latex/` derivations: `N` = number of subcarriers = **delay axis = rows**
(`n`, `l`), `M` = number of OFDM symbols = **Doppler axis = columns** (`m`, `k`). Both the
time–frequency grid `X[n,m]` and the delay–Doppler grid `x[l,k]` are `N×M` matrices. Raviteja's
paper uses the opposite letters, so every number lifted from it must be translated:

| Quantity | Raviteja TVT'19 / WCNC'18 | Liu 2026 | This library |
|---|---|---|---|
| subcarriers / delay bins | `M = 512` | `N = 128` | `N` |
| OFDM symbols / Doppler bins | `N = 128` | none (single block) | `M` (128 for EVA; 1 or 16 for MT40) |
| DD grid | `x[k,l]`, Doppler rows | vector `s` | `x[l,k]`, delay rows (transpose of Raviteja) |
| delay tap of path `i` | `l_τi` (integer by assumption) | `τ_l` in samples, continuous | `l_p = round(τ_p N Δf)`; guard bound `l_τ = ceil(τ_max N Δf)` |
| Doppler tap | `k_ν = 4` at 120 km/h; `ν = (k+κ)/(NT)` | — | `k_ν = ceil(ν_max M T_S)`; `k_real = ν M T_S` (fractional kept) |
| full-guard pilot overhead | `(2l_τ+1)·N = 5248` cells | — | `(2l_τ+1)·M = 5248` cells (8.0 %) |
| DD channel kernel (WCNC eq. 13) | Dirichlet `/N`, twist phase `e^{−j2π k_real l/(NM)}` | — | Dirichlet `/M`, **no twist** (the demodulator is referenced to each symbol's own start, so the effective gain equals the physical gain), `T_S` inside the Doppler kernel |
| path gain | `h_i` | `β_l` | `α_p` |
| estimated parameters | — | `[τ; Re β; Im β]` (3L) | `η_p = [α_R; α_I; τ; ν]` stacked to 4P; the `ν` block is masked out when `ν ≡ 0` |
| noise / SNR | `σ² = 1`, signal scaled | `SNR := 1/σ²` | unit-power data symbols, `N₀ = 1/SNR` per bin (monostatic: sensing noise = comm noise) |
| range unit | — | monostatic, `c/(2B) = 0.9375 m` per sample | same for MT40; **bistatic `c·τ` for EVA** (as in `latex/main.tex`). The two differ by 4× in m²; every m² axis states which. |

Timing constants: `T = 1/Δf` (useful symbol duration), `T_s = 1/(NΔf)` (sample period),
`T_S = T + T_cp` for the per-symbol-CP waveforms and `T_S = T` for the reduced-CP (frame-CP)
waveforms. Delay resolution `T_s`; Doppler resolution `1/(M T_S)`.

## 2. EVA (9 paths, delay + Doppler)

### 2.1 Numerology

| Parameter | Value | Note |
|---|---|---|
| carrier `f_c` | 4 GHz | Raviteja Sec. IV. (`latex/eva_channel.tex` shows a 2 GHz / 222 Hz worked example — stale; `matlab/CLAUDE.md`'s "540 Hz" is wrong.) |
| speed `v` | 120 km/h | Raviteja default; 30 and 500 km/h appear only in his Fig. 7/13 |
| `ν_max = v f_c / c` | 444.44 Hz | `c = 3e8` |
| `N` (subcarriers) | 512 | Raviteja's `M` |
| `M` (OFDM symbols) | 128 | Raviteja's `N` |
| `Δf` | 15 kHz | |
| `T`, `T_s`, `B` | 66.667 µs, 130.21 ns, 7.68 MHz | |
| `T_cp` | `l_τ·T_s = 20 T_s = 2.604 µs` | equals Raviteja's OFDM CP (WCNC companion, 2.6 µs); makes the CP-model assumption `τ_max ≤ T_cp` hold. User-confirmed default; `T_cp = 0` reproduces the legacy `old/eva` numbers. |
| `T_S` (CP waveforms) | 69.27 µs | Doppler bin `1/(M T_S) = 112.8 Hz` (117.19 Hz at `T_cp = 0`); frame `8.87 ms` |
| `ν_max T` | 0.0296 | ICI strength: about −25 dB, "small but present" |
| `k_ν` | 4 | `ceil(444.44·128·69.27e−6) = ceil(3.94)` — matches Raviteja |
| `l_τ` | 20 | `ceil(2510e−9 · 7.68e6) = ceil(19.28)`. A guard **bound**, not an occupied tap (`round` gives 19). |
| modulation | Gray 4-QAM, unit power | 16-QAM only for Raviteja's Fig. 12/14 |
| SNR grid | `0:1:25 dB` (`N₀ = 1/SNR`) | Raviteja Fig. 9 check uses `SNR_d = 10:2:18` with his convention `σ² = 1` and scaled signal |

### 2.2 Power-delay profile (3GPP EVA)

| path | delay [ns] | rel. power [dB] | linear | normalized `P_l` | tap `l_p` (round, N=512, Δf=15 kHz) |
|---|---|---|---|---|---|
| 1 | 0 | 0.0 | 1.00000 | 0.24120 | 0 |
| 2 | 30 | −1.5 | 0.70795 | 0.17076 | 0 |
| 3 | 150 | −1.4 | 0.72444 | 0.17473 | 1 |
| 4 | 310 | −3.6 | 0.43652 | 0.10529 | 2 |
| 5 | 370 | −0.6 | 0.87096 | 0.21007 | 3 |
| 6 | 710 | −9.1 | 0.12303 | 0.02967 | 5 |
| 7 | 1090 | −7.0 | 0.19953 | 0.04813 | 8 |
| 8 | 1730 | −12.0 | 0.06310 | 0.01522 | 13 |
| 9 | 2510 | −16.9 | 0.02042 | 0.00492 | 19 |

Sum of linear powers `4.14594` (+6.18 dB); `P_l = linear / 4.14594` so that `Σ P_l = 1`
(mandatory here; the tex calls the normalization optional — with it the per-bin SNR is the SNR
axis). `τ_max = 2.51 µs`, RMS delay spread ≈ 357 ns.

**Tap collision.** At this grid the 0 ns and 30 ns paths both round to tap 0 (and 310/370 ns to
taps 2/3). The two colliding paths are separable only through their different Dopplers, which
makes some channel draws nearly singular for the bounds. `run/diag_cond_eva.m` measures how often;
`chan_eva(...,'merge_colliding',true)` sums the colliding gains into one tap (`P = 8`) as a
documented alternative.

### 2.3 Per-trial draw

- **Gains** (user-confirmed default): `α_l ~ CN(0, P_l)`, i.i.d. across paths, redrawn per trial
  (`gain_model = 'rayleigh'`). Toggle `'constmag'`: `|α_l| = sqrt(P_l)` with a uniform random phase
  — the `old/eva` default, known to overstate OTFS's diversity gain. Within a trial the gains are
  deterministic unknowns, so every bound is conditional on the draw and then trial-averaged; no
  Bayesian bound is used anywhere.
- **Doppler**: one shift per path, `ν_l = ν_max cos θ_l`, `θ_l ~ U[−π, π)` i.i.d. (Jakes). The
  fractional part of `k_real = ν_l M T_S` is kept (Raviteja's "fractional Doppler" case, his Fig.
  9). `doppler_scale` multiplies `ν_max` (0 gives the no-Doppler control).
- **Delays**: snapped to the tap grid, `τ_p = l_p T_s` with `l_p = round(τ_p^{phys} N Δf)`
  (`snap = true`, Raviteja Sec. II-C "nearest sampling point"). This is what the DD receiver
  assumes, and the bounds are evaluated on the same snapped channel so that the detector-error
  MCRB and the BER refer to one realization. `snap = false` (physical delays, off-grid) is
  supported by the CP builders only and is never mixed with `snap = true` in one figure.

### 2.4 Receiver constants recorded with the channel

Pilot at `(l_p, k_p) = (l_τ, M/2)`; full-Doppler guard rows `0 … 2l_τ` (5248 cells, 8.0 %);
threshold CE at `3σ_p`; DD Doppler spread truncated to `|q| ≤ E = 10` (Raviteja Fig. 2, WCNC Fig.
3); OFDM pilot lattice matched to the OTFS pilot+guard cell count within 5 % (exact matching is
impossible at full EVA: nearest attainable lattice has 5472 cells) with the pilot **energy** matched
exactly and data-referenced.

### 2.5 Variants

| Variant | Definition | Used for |
|---|---|---|
| **EVA-reduced** | 3 strongest, well-separated taps, delay-sorted: `[0 150 370] ns`, `[0 −1.4 −0.6] dB`, taps `[0 1 3]`, normalized to unit power; same Doppler/gain rules | well-conditioned control; the bound workhorse if the conditioning gate in `run/diag_cond_eva.m` fails for the 9-path channel (this is what `old/eva` ended up using) |
| **legacy regression preset** | the *power-sorted* pairing of `old/eva/verify_crb.m`: `τ = [0 370 150] ns`, `P = [0 −0.6 −1.4] dB`, fixed Jakes angles `θ = [0.5 2.0 3.5]`, fixed phases `φ = [0.7 1.1 −0.5]` assigned **by index**, `f_c = 2 GHz` (`ν_max = 222.2 Hz`), `N = 16, M = 8`, `T_cp = 0`, un-normalized gains | 1e-12 regression tests against `old/eva/utils`. Never delay-sort it: re-pairing `(τ, θ, φ)` gives a different channel. |
| **EVA-reduced-fast** | EVA-reduced at **240 km/h** (`ν_max = 888.9 Hz`, `k_ν = 8`, `ν_max T = 0.059`); same delay profile, same draw rules. `chan_eva(g,'variant','reduced-fast')`, or `'channel','eva-reduced-fast'` in `run_ber_eva`. | the regime where CP-OFDM's **lower bound beats CP-OTFS's in both the delay and the Doppler block** while the genie CRBs stay equal and OTFS keeps the lower BER (`eva_reduced_fast_results.md`) |
| custom profile | `'taus_ns'`, `'pows_db'`, `'name'`: any delay/power profile with the EVA draw rules | ad-hoc variants (collisions, stretched delay spread) |
| `merge_colliding` | colliding taps summed, `P = 8` | conditioning study |
| `doppler_scale = 0` | `ν ≡ 0`, no ICI | "Doppler as sensing diversity" comparison |
| `constmag` gains | fixed magnitudes, random phases | reproducing old figures |

## 3. MT40 (40 static targets, delay only)

### 3.1 Numerology (Liu Sec. V-C)

| Parameter | Value | Note |
|---|---|---|
| targets `L` | 40 | `L = 10` sub-variant for MCRB at `p > 0` (the 3L = 120-of-128 case is at the identifiability edge) |
| `N` | 128 | subcarriers = block samples = data symbols |
| bandwidth `B` | 160 MHz | `Δf = B/N = 1.25 MHz`, `T_s = 6.25 ns`, `T = 800 ns` (derived; the paper states only `N` and `B`) |
| carrier | not specified | irrelevant: no Doppler |
| delays | continuous, distinct, in `[0, 127]` samples | see 3.2 |
| gains `β_l` | `= 1` for all `l` | Liu requires only `β_l ≠ 0`. The absolute m² levels of his Fig. 1 depend on this repo convention. Toggles: `'randphase'`, `'cn01'`. |
| Doppler | none, `ν ≡ 0` | the Doppler block of the FIM is structurally null in the no-ICI model; estimated set `[α_R, α_I, τ]` |
| modulation | 16-PSK **and** 16-QAM | kurtosis 1 vs 1.32: with PSK, OFDM attains the Jensen bound exactly and the waveform gap vanishes; with 16-QAM the gap is strict. 4-QAM (kurtosis 1) is degenerate like PSK and is labelled so. |
| SNR | `SNR := 1/σ²`, `−10:1:6 dB` | Liu Fig. 1 axis |
| range unit | `c/(2B) = 0.9375 m` per sample, `CRB[m²] = 0.87890625·CRB[samples²]` | monostatic |

The paper's Fig. 1 plots `Tr{T J⁻¹}` **summed over the 40 targets** (≈ 4.85 m² for OFDM at
−10 dB is a sum, not a per-target value). The Jensen bound `Tr{T J̄⁻¹}` with
`J̄ = (2/σ²)Re{HᴴH}` is waveform- and constellation-independent and scales exactly as `1/SNR`.

### 3.2 Geometry (fixed, not redrawn per trial)

```matlab
rng(20260514);                      % MATLAB default generator
usable   = 127 - (L-1);             % = 88 for L = 40
u        = sort(usable * rand(L,1));
tau_samp = u + (0:L-1)';            % >= 1-sample separation, max <= 127
tau      = tau_samp * T_s;          % seconds
beta     = ones(L,1);
```

Bit-identical to `impl/2605.14451v1/run_sim.m` and to the old `'paper'` scenario, so results are
comparable across the three code bases. Liu requires only distinctness; the one-sample separation
is a repo convention. Monte Carlo is over the **data symbols** (and over gains only for the
`'randphase'`/`'cn01'` toggles).

### 3.3 Framings

| Framing | Grid | Purpose |
|---|---|---|
| **MT40-1D** | `M = 1`, one 128-sample block; Liu's unitary precoders `Qᴴ = F_N U` — OFDM `Q = I`, SC `Q = F_N`, OTFS `Q = F_N (F₃₂ᴴ ⊗ I₄)` (`N₁ = 32, N₂ = 4`), AFDM `F_N Λ_{1/16}ᴴ F_Nᴴ Λ_{1/8}ᴴ` (optional) | validation against Liu Fig. 1; CRB and BSC-MCRB at `p = 0` |
| **MT40-2D** | `N = 128 × M = 16` frame; per-symbol CP `T_cp = l_τ T_s` with `l_τ := N − 1 = 127` (delays may span the block), so `T_S = 255 T_s = 1.594 µs` (about 50 % overhead; the CP assumption then holds by construction — Liu's model is circular by fiat instead). `T_cp = 0` reproduces the old `paperv2`-style circular block. | the CP/RCP OFDM/OTFS case matrix of the derivations on a delay-only channel (ν masked). The 1-D and 2-D framings are on different SNR-per-frame axes because `T_S` differs. |
| **MT40-D** | MT40-2D plus per-target Doppler `rng(20260707); ν_l = ν_max (2·rand − 1)`, `ν_max = 0.5/(M T_S)` for the active `T_S`, ICI on | the only MT40 form that exercises the Doppler half of the derivations |
| **MT40-int** | delays rounded to samples | required by the sampled RCP-OTFS branch (derived for sample-spaced delays only) |
| **MT40-10** | `L = 10` | MCRB at `p > 0` (40 targets is too ill-conditioned once symbols are in error) |

MT40 has **no BER simulation** (user decision): its delay spread covers the whole block, so a
Raviteja-style embedded-pilot guard would need more than `N` cells. Its mismatch in Phase 3 comes
from the abstract BSC / symbol-error model only.

## 4. Idealized controls (tests only)

| Control | Definition | Checks |
|---|---|---|
| AWGN | single unit-gain tap, `τ = ν = 0`, random phase | OFDM = OTFS = Gray 4-QAM theory |
| flat Rayleigh | one tap `α ~ CN(0,1)`, `τ = ν = 0` | OFDM = OTFS on the order-1 Rayleigh curve |
| EVA-reduced, `ν = 0` | as above | order-1 slope for OFDM, delay diversity for OTFS |

## 5. What the bounds and the receivers assume

- Gains are deterministic unknowns within a trial (all `η_p` entries are estimated jointly).
- The closed-form expected FIM/MCRB assume unit-modulus i.i.d. data on the **full** grid; embedded
  pilots and guards invalidate them, so the detector-error bounds are computed from the realization
  builders on the actually transmitted grid.
- `E` (DD spread truncation) is a receiver approximation, never part of a bound.
- The DD receiver's tap model is derived from the same io layer that generates the received signal;
  it differs from Raviteja's WCNC eq. (13) by exactly the twist phase and by `T_S` versus `T`.

## 6. Which script uses which channel

| Script | Channel / variant | Seed | Notes |
|---|---|---|---|
| `run/diag_cond_eva.m` | EVA full and reduced, both gain models | 9000 | conditioning gate |
| `run/run_crb_mt40_block.m` | MT40-1D, 16-PSK and 16-QAM | 20260514 (geometry), 7000 (data) | Liu Fig. 1 |
| `run/run_crb_eva.m` | EVA (gated variant) | 9000 | CRB vs SNR, all waveforms |
| `run/run_crb_mt40_frame.m` | MT40-2D (+ MT40-D, MT40-int, MT40-10) | 20260514 / 20260707 / 7000 | frame case matrix |
| `run/run_ber_fig9.m` | EVA full, Rayleigh, Raviteja noise convention | 2019 | Raviteja Fig. 9 check |
| `run/run_ber_eva.m` | EVA (+ controls AWGN, flat Rayleigh, EVA-reduced with `doppler_scale = 0`) | 9100, substream `1000 + 10·SNR_dB` | main BER |
| `run/run_mcrb_bsc_eva.m` | EVA-reduced (default) or full EVA (`'reduced', false`) | 20260707 | MCRB/LB vs p (BSC) |
| `run/run_mcrb_rx_eva.m` | EVA-reduced / full EVA; same seed and substream key as `run_ber_eva` (identical draws at equal SNR); `'legacy', true` = old `mcrb_ce_ici` configuration (constmag, `T_cp = 0`, posterior MP, seed 9000) | 9100 / 9000 | detector-error MCRB/LB |
| `run/run_mcrb_bsc_mt40.m` | MT40-2D (`'framing', '2d'`, default `L = 10`, plus `L = 40` at `p = 0`) and MT40-1D with Liu's OFDM/OTFS/SC precoders (`'framing', '1d'`) | 20260514 (geometry), 7000 (data) | MCRB/LB vs p, vs SNR |

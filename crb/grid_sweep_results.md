# Can the grid (N, M) widen OFDM's delay margin over OTFS?

Short answer: **no.** The delay margin is flat in M and non-monotone in N with a maximum at the
512 baseline. But the reason is subtler than "the grid changes nothing", and one of the parallel
theory predictions turned out to be right where a first reading of the data suggested otherwise.

Setup: channel `EVA-reduced-fast` (three largest-gain EVA paths at 240 km/h), Δf = 15 kHz,
`T_cp` = 20 samples, 4-QAM, the standard receivers, 40 trials per point, SNR 0/4/8/12 dB, paired
per-trial OTFS/OFDM ratios (> 1 = OFDM wins). `run_mcrb_rx_eva` gained `'N'`, `'M'`, `'deltaf'`.

> **Units warning.** The `LB` ratios below are ratios of *traces*, i.e. **variance (MSE) ratios**.
> The RMSE ratio is their square root. A delay MSE ratio of 1.6 is a 1.27× RMSE advantage. The
> same applies to the ratio columns of `eva_reduced_fast_results.md`.

## Measurement (geometric mean over the four SNRs)

| grid | l_τ | k_ν | delay margin | Doppler margin | OFDM anisotropy | OTFS anisotropy | gap |
|---|---|---|---|---|---|---|---|
| 256 × 128 (taps collide) | 2 | 8 | 0.88 | 3.04 | 1.98 | 4.38 | 1.86 |
| 384 × 128 | 3 | 8 | 1.08 | 4.26 | 1.62 | 3.92 | 1.98 |
| 512 × 32 | 3 | 2 | 1.67 | 15.92 | 0.81 | 3.61 | 3.09 |
| 512 × 64 | 3 | 4 | 1.36 | 12.93 | 0.83 | 3.99 | 3.08 |
| **512 × 128 (baseline)** | 3 | 8 | **1.62** | 13.49 | 1.04 | 5.12 | 2.89 |
| 512 × 256 | 3 | 16 | 1.27 | 16.09 | 0.95 | 5.64 | 3.57 |
| 1024 × 128 | 6 | 8 | 1.15 | 11.18 | 1.15 | 6.57 | 3.12 |

Anisotropy = (bias in Doppler bins)/(bias in delay samples); `gap` = sqrt(Doppler margin / delay
margin), how much more decisively OFDM wins Doppler than delay.

## What the grid does and does not control

1. **The delay margin is flat in M**: slope −0.09 per e-fold over an 8× range (1.67, 1.36, 1.62,
   1.27), against +1.0 predicted by the pseudo-true algebra. That is an order-of-magnitude
   disagreement, not a noise-level one, so the flatness itself is safe to quote.
2. **The Doppler margin is also flat in M** (slope +0.01) but its *absolute* level is large
   everywhere (11 to 16), so the imbalance persists on every grid: the gap is 2.9 to 3.6 on all
   four grids where OFDM wins delay at all.
3. **N is non-monotone with a maximum at 512.** Both margins peak there and fall on either side.
   The fall below 512 is caused by OFDM, not OTFS: OFDM's anisotropy breaks from ~1 to 1.62 and
   1.98 at N = 384 and 256, and its Doppler bias jumps 2.8× for a 1.33× grid change, because the
   pilot lattice carries fewer frequency pilots and the decision-directed estimate degrades. At
   N = 256 the taps collide ([0 1 1]) and OFDM loses the delay block outright.
4. **The anisotropy is *not* grid-invariant** — it obeys a clean law,
   **`R_OTFS = 0.17·N^0.34·M^0.23`** (refitted here: +0.229 in M, +0.338 in N). It rises because
   OTFS's delay bias falls faster with the grid (N^−0.98 over 512→1024) than its Doppler bias
   (N^−0.62): more subcarriers separate the taps in delay ([0 1 3] → [0 2 6]), suppressing the
   cross-path leakage that feeds OTFS's delay bias. **Enlarging the grid closes the delay margin
   rather than opening it.**
5. **Why the margins are nonetheless flat.** Every bias here — both waveforms, both blocks — falls
   with the grid at roughly `(NM)^−1/2` (measured exponents in M: −0.60, −0.50 for OFDM delay and
   Doppler, −0.72, −0.49 for OTFS). Each margin is therefore a ratio of two identically-scaling
   laws, i.e. a ratio of *prefactors*, which the grid does not set. Two different mechanisms give
   the same exponent: for the incoherent channels the bias is a random walk of `p·NM` error signs
   against a ramp, giving `sqrt(p/(NM))` with OFDM's SER grid-invariant (0.129 to 0.134 at 8 dB on
   all seven grids); for OTFS's coherent Doppler channel the leftover fractional-Doppler replica
   scales as `1/S` with the retained tap count `S ∝ sqrt(NM·SNR)`.
6. **The delay block has no coherent mis-modelling channel at all**, because the delays are
   snapped and the delay response is therefore exact. OTFS's delay bias sits on the same
   incoherent floor as OFDM's, and the whole delay margin is a structure factor of 0.78 to 1.47.
   That is the real reason it cannot be tuned by the grid — and it also says what *would* work:
   create a coherent delay mis-modelling channel, i.e. unsnap the delays.

## Predictions versus measurement

Three derivations were run in parallel with the sweep.

- **Delay-Doppler support geometry.** Its closed-form geometric law (`anisotropy ∝ N^0.15 M^-0.29`)
  failed, but the empirical law it fitted from its own delay-Doppler measurements,
  `R_OTFS = 0.17 N^0.372 M^0.264`, **survived on both exponents** (measured +0.23 and +0.34) and
  was uniformly 10 to 30 % high in level on all seven grids. It also predicted correctly that
  N = 256 is the one grid where the delay margin drops below parity, and that at 512 × 32 the OTFS
  Doppler bound would be absurd in absolute terms (predicted 55 Hz against a 4.1 Hz genie CRB;
  measured 52.2 against 4.02). Its one failure was the recommendation: it expected the margin to
  grow with `NM`, and the margin falls.
- **Pseudo-true algebra.** Its central law `margin_τ ∝ sqrt(M/N)` failed with the wrong sign and
  a tenfold magnitude error (predicted ×2.83 over the M range, measured ×0.87). Its microscopic
  prediction `bias_τ,OTFS ∝ 1/l_τ ∝ 1/N` was right (measured exponent −0.98).
- **Detector reasoning.** Its headline — both margins flat in the grid — was right, and its
  prediction of the Doppler anisotropy level (6.4, range 2.4 to 18) matched almost exactly
  (measured median 6.22, range 3.2 to 12.5).

## Recommendation

**Keep 512 × 128.** It carries the largest delay margin of any grid that preserves the story
(1.62, statistically tied with 512 × 32's 1.67), and 8 dB on the 512 × 128 and 512 × 256 grids are
the only two points in the whole sweep where OFDM wins the delay block in 80 % of trials
(p = 1.8e-4). The alternatives all cost more and deliver less: 512 × 256 doubles the frame to
17.7 ms for a *smaller* delay margin, 1024 × 128 doubles the bandwidth to 15.4 MHz and lowers both
margins, and 512 × 32 ties on delay but makes OTFS lose on BER as well and produces a 52 Hz
Doppler bound against a 4 Hz CRB.

**To equalise the two advantages the grid is the wrong knob.** The delay margin is a structure
factor near unity that the grid does not move; raising it requires creating a coherent delay
mis-modelling channel. The two untested-in-combination candidates are unsnapped delays
(`chan_eva('snap', false)`, which reached a 3.7× delay margin at 2 dB but broke OFDM's estimator
above 10 dB) and higher `ν_max`, which is what lifted the delay margin from ~1.0 to ~1.6 in the
first place.

## What 40 trials support, and what they do not

Supported: flatness in M (an order-of-magnitude disagreement with the competing prediction); the
Doppler margin's grid-invariance (12 % spread across 8× in M and 2× in N); every Doppler sign test
at 0 to 8 dB (p ≤ 6e-3, usually ≤ 1e-7); OFDM's SER invariance; and the `R_OTFS` exponents, which
are measured on paired within-run ratios.

Not supported: **any single grid-to-grid delay comparison.** Only 7 of 28 delay-block (grid, SNR)
points clear p < 0.05, and 4 of those are at 8 dB; win fractions of 0.65 and 0.60 correspond to
p = 0.08 and p = 0.27. So "512 × 64's 2.70 at 8 dB beats 1024 × 128's 1.46" is *not* a supported
claim — "the delay margin is 1.0 to 1.7 everywhere and does not trend with M" is. Also not
supported: the N axis as a law (four points, one octave at the top, confounded by the pilot
lattice being re-optimised at each N), the 12 dB column anywhere, and any extrapolation outside
256 ≤ N ≤ 1024, 32 ≤ M ≤ 256.

Raw logs: session scratchpad `nm_sweep.log`, `nm_sweep2.log`. Result files
`results/mcrb_rx_eva_reduced_fast_n40_nm*.mat` and `..._nmlo*.mat`.

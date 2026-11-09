# Our MT40 CRB results vs Liu et al., arXiv:2605.14451, Fig. 1

Compared on 2026-09-08. "Paper" values are the exact curve data extracted from the vector
graphics of Fig. 1 (page 8 of `papers/2605.14451v1.pdf`, 17 points per curve, precision
~0.001 m²), not pixel read-offs. "Ours" = `results/crb_mt40_block.mat` (`run_crb_mt40_block`,
500 data draws, mean over draws, `β = 1`, geometry `rng(20260514)`, no draw dropped);
"earlier" = the independent reproduction `impl/2605.14451v1/crb_new_results.mat`. All three use
`CRB = Tr{T E[J⁻¹]}` summed over the 40 targets, `c²/(4B²) = 0.87890625 m²/sample²`,
`SNR = 1/σ²`, `N = 128`, `N₁ = 32`, `N₂ = 4`, `c₁ = 1/16`, `c₂ = 1/8`. Every curve in every data
set is exactly proportional to `1/SNR`, so the −10 dB column carries all the information.

## Absolute levels at −10 dB (m², sum over 40 targets)

| curve | paper | ours (mean) | ours (median) | earlier | ours/paper |
|---|---|---|---|---|---|
| 16-QAM SC | 6.357 | 7.218 | 7.129 | 7.145 | 1.135 |
| 16-QAM AFDM | 6.302 | 7.066 | 6.936 | 7.032 | 1.121 |
| 16-QAM OTFS | 6.062 | 6.909 | 6.843 | 6.892 | 1.140 |
| 16-QAM OFDM | 4.835 | 5.466 | 5.437 | 5.399 | 1.131 |
| 16-PSK SC | 6.313 | 7.190 | 7.072 | 7.143 | 1.139 |
| 16-PSK AFDM | 6.231 | 6.919 | 6.838 | 6.962 | 1.110 |
| 16-PSK OTFS | 5.953 | 6.702 | 6.548 | 6.720 | 1.126 |
| 16-PSK OFDM | 3.877 | 4.330 | 4.330 | 4.330 | 1.117 |
| Jensen (both panels) | 3.877 | 4.330 | 4.330 | 4.330 | 1.117 |

Every curve, the deterministic Jensen line included, is 1.11–1.14× the paper: one common scale
factor (the curve-to-curve spread of the ratio is at the ~1 % Monte-Carlo level). The paper does
not specify the target placement rule or the gains; the Jensen constant is dominated by the closest
target pair and is heavy-tailed over placements (under our rule, 5000 placements: median 0.94 m²,
p95 4.09 m², max > 500 m²; distinct-integer placement gives a median of 557 m²). Both the paper's
3.88 m² and our 4.33 m² sit near the 94th–96th percentile of that distribution, so the 12 % offset
cannot be attributed without the authors' delays. The dB gains below are placement-invariant
(60 placements: ±0.03–0.05 dB, correlation with the Jensen level 0.02).

## OFDM gain, dB = 10·log10(CRB_X / CRB_OFDM), SNR-independent

| | paper | ours (mean, 500) | ours (median) | fresh 4000-draw run | earlier |
|---|---|---|---|---|---|
| 16-QAM vs SC | 1.189 | 1.207 | 1.177 | 1.194 | 1.217 |
| 16-QAM vs AFDM | 1.151 | 1.115 | 1.057 | 1.108 | 1.148 |
| 16-QAM vs OTFS | 0.982 | 1.017 | 0.999 | 1.015 | 1.060 |
| 16-QAM vs cluster mean | 1.108 ("≈ 1 dB") | 1.114 | 1.078 | – | 1.142 |
| 16-QAM OFDM above Jensen | 0.960 | 1.012 | – | 0.986 | 0.958 |
| 16-PSK vs SC | 2.118 | 2.202 | 2.130 | 2.170 | 2.174 |
| 16-PSK vs AFDM | 2.061 | 2.035 | 1.984 | 2.047 | 2.062 |
| 16-PSK vs OTFS | 1.863 | 1.897 | 1.796 | 1.920 | 1.909 |
| 16-PSK vs cluster mean | 2.016 ("≈ 2 dB") | 2.047 | 1.972 | – | 2.050 |
| 16-PSK OFDM above Jensen | 0.000 | 0.000 (1e-12) | 0.000 | 0.000 | 0.000 |

Ours minus paper: 16-QAM +0.02 / −0.04 / +0.04 dB (+0.05 for OFDM above Jensen); 16-PSK
+0.08 / −0.03 / +0.03 dB — 1–2 σ of the Monte-Carlo error (~0.03 dB at 500 draws) plus the
placement scatter (0.03–0.05 dB). The acceptance windows of `run_crb_mt40_block` (0.5–1.8 dB and
1.3–2.7 dB) pass with wide margin.

## Claims of the paper reproduced exactly

- Corollary 1: OFDM with 16-PSK equals the Jensen bound on every draw (max relative deviation
  1e-12); in the paper the OFDM-PSK curve sits on the Jensen markers to 1e-4.
- Proposition 1: the Jensen bound is waveform- and constellation-independent (same value in both
  panels) and exactly ∝ 1/SNR; all curves are ∝ 1/SNR in the paper too.
- Ordering SC > AFDM > OTFS > OFDM ≥ Jensen at every SNR and both constellations, following
  α_SC > α_AFDM > α_OTFS > 0 (spread factors 0.28867 / 0.28641 / 0.27951 reproduced to 7 digits;
  eq. 147 holds with equality for all three bases). Jensen bounds the expectation only: individual
  16-QAM OFDM draws dip up to 1.8 % below it (2 of 500), which the paper does not exclude.
- "SC, AFDM and OTFS show close CRB levels": SC/OTFS spread 4.5 % ours vs 5.4 % paper (16-QAM),
  7.3 % vs 6.5 % (16-PSK). "Only marginal improvement from 16-QAM to 16-PSK" for the spread
  waveforms: paper −0.7 / −1.1 / −1.8 %, ours −0.4 / −2.1 / −3.0 % (SC / AFDM / OTFS), same ordering.
- No singular FIM in 500/500 draws (Sec. V-A: probability 16^−(N−1)); kurtosis 1.32 / 1.00 (Table I).

## Where the shapes differ

- SC–AFDM separation: ours 0.09 dB (16-QAM) / 0.17 dB (16-PSK) vs the paper's 0.04 / 0.06 dB. The
  second-order gap of eq. 145 weights every frequency pair by a geometry-dependent μ_nm, so this
  separation depends on the placement; it is within ~1.5 σ of the placement scatter of the gains.
- QAM/PSK gap ratio: 0.61–0.63 in our data and 0.59–0.63 in the paper's own figure, against the
  0.68 = (2 − κ_QAM)/(2 − κ_PSK) of the leading-order theory (both data sets show the same finite-N
  shortfall; at N = 512 ours gives 0.64–0.67).

## The paper's large-sample formula at N = 128

Evaluating the paper's second-order term `E{Tr(T R₂)} = Tr{W Σ_v(U)}` (eqs. 124–125, 140, 145)
with our geometry (two independent constructions agree to 6 digits) predicts OFDM gains of
1.12 / 1.01 / 0.88 dB (16-QAM) and 1.76 / 1.60 / 1.40 dB (16-PSK), i.e. 0.1–0.5 dB below both
the paper's figure and our numbers, and covers 70–84 % of the measured gap vs OFDM. The reason is
that at N = 128 with 3L = 120 unknowns the resolvent expansion is outside its convergence region:
the spectral radius of `J̄^{-1/2} Δ J̄^{-1/2}` is ≈ 2 for SC/AFDM/OTFS on every draw (0.73 for
OFDM-QAM, 0 for OFDM-PSK), so the event Ω_N of Proposition 2 essentially never occurs, although
the exact identity `J⁻¹ = J̄⁻¹ − R₁ + R₂ − R₃` still holds and R₂ carries most of the gap. At
N = 512 (same relative geometry) the same formula is accurate to 1–7 % and the relative gaps
scale as ≈ 1/N (ratios 3.4 and 4.5 between N = 128 and 512), consistent with the paper's Θ
statements. This is a remark about the regime of the paper's theory, not about the code.

## Does "OFDM beats OTFS" carry over to the EVA channel?

No measurable gap. On the EVA sweep (`results/crb_eva_smoke.mat`: 512 × 128 grid, 9 Rayleigh
paths with Doppler, Gray 4-QAM, 20 channel + data draws) the per-draw ratio of the CP-OTFS to the
CP-OFDM genie CRB trace is 1.0002 (mean) / 0.9998 (median) with range 0.991–1.011 for delay and
0.9999 / 1.0004 (0.991–1.006) for Doppler, with or without ICI; OTFS is worse in exactly half of
the draws. Two reasons, both predicted by Liu's own analysis: (i) 4-QAM has kurtosis 1, so OFDM
attains its Jensen-type bound on every draw (per-draw OFDM CRB = closed-form `(E J)⁻¹` to 3e-5);
(ii) the OTFS excess scales as `(2 − κ) α² / (number of samples)`, and here there are 65 536
samples for 36 unknowns (MT40: 128 samples for 120 unknowns), so the second-order term is buried
under the ±1 % scatter of the ISFFT-spread symbol energies. The 8 % higher Doppler CRB of the
reduced-CP waveforms relative to CP-OFDM is the shorter symbol period `T` vs `T + T_cp`
(ratio 1.079), not a waveform effect. Liu's statement is therefore a small-sample / identifiability
-edge effect and is irrelevant on a communication-sized frame such as EVA.

### Can the OFDM advantage be brought back by shrinking the grid or adding targets?

Yes, in the direction Liu's theory predicts, but the effect is set by the ratio of unknowns to
samples and stays small for anything frame-sized (`scratchpad/eva_edge.m`, 100 draws each,
median of the per-draw ratio of the CP-OTFS to the CP-OFDM CRB trace, [p10, p90] in brackets).

Knob A — shrink the grid at fixed bandwidth 7.68 MHz; EVA with the colliding taps merged
(P = 8, 36 unknowns), Rayleigh gains, Doppler rescaled so that `k_ν ≈ 4` on every grid:

| grid | samples | 4-QAM delay | 4-QAM Doppler | 16-QAM delay | draws with OTFS worse (4-QAM) |
|---|---|---|---|---|---|
| 512 × 128 | 65 536 | 0.999 [0.995, 1.004] | 1.001 | 1.000 | 45 % |
| 128 × 32 | 4 096 | 1.002 [0.989, 1.023] | 1.005 | 1.001 | 57 % |
| 64 × 16 | 1 024 | 1.017 [0.981, 1.056] | 1.015 | 1.010 | 67 % |
| 32 × 8 | 256 | 1.052 [0.985, 1.146] | 1.056 | 1.027 | 81 % |

Knob B — more targets on a 64 × 8 grid (512 samples), delay-only, random continuous delays with
≥ 0.5-sample separation, Rayleigh gains (rows with `cond > 1e18` are pinv-truncated and only
indicative):

| P | unknowns | 4-QAM delay | 16-QAM delay | draws with OTFS worse (4-QAM) | cond (median) |
|---|---|---|---|---|---|
| 4 | 12 | 1.002 [0.960, 1.052] | 1.011 | 55 % | 3e2 |
| 8 | 24 | 1.027 [0.974, 1.083] | 1.016 | 71 % | 2e3 |
| 16 | 48 | 1.057 [0.990, 1.115] | 1.018 | 85 % | 7e5 |
| 32 | 96 | 1.09 | 1.05 | 89 % | 2e18 |
| 100 | 300 | 1.14 [1.10, 1.19] | 1.09 | 100 % | 4e19 |

Reading: the median OTFS penalty is roughly 0.3–0.6 × (unknowns / samples) for constant-modulus
data (36/256 → 5 %, 48/512 → 6 %, Liu's 120/128 → 32 %), smaller for 16-QAM because OFDM itself
then sits above the Jensen bound, and it is a statistical advantage: even at 32 × 8 OTFS is the
better waveform in one draw out of five. Reaching Liu's ~1 dB needs the number of unknowns to
approach the number of samples, i.e. a sensing-only block with tens of targets, not a
communication frame.

## Verdict

The reproduction matches everything Fig. 1 and Sec. V-C assert: both identities hold exactly,
the ordering and the constellation dependence are reproduced, and every OFDM gain agrees with the
figure to 0.01–0.08 dB. The only unexplained difference is a common 12 % scale factor on the
absolute m² levels, which comes from the paper's unspecified target placement/gain convention and
cannot be resolved without the authors' geometry. Scripts of the analysis: session scratchpad
`liu_check.m` … `liu_check4.m`, `verify/extract_paths.py`, `verify/fig1_values.py`.

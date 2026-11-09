# Our EVA BER results vs Raviteja, Phan, Hong, IEEE TVT 2019, Fig. 9

Compared on 2026-09-09. Fig. 9 of the paper (`papers/Embedded_Pilot-Aided_…pdf`, page 8):
BER vs `SNR_d`, fractional Doppler, full guard, 4-QAM, EVA at 120 km/h, curves = embedded-pilot
CE + MP at `SNR_p = 40 / 45 / 50 dB`, OTFS with known channel + MP, OFDM with known channel.

**Paper values** are exact: the figure is vector graphics in the PDF, and the curve vertices were
extracted from the page content stream (precision ~0.1 %). The four OTFS curves have vertices at
10, 15 and 18 dB only; the values at 12/14/16 dB (marked \*) are the log-linear chords the authors
drew, not simulated points. The OFDM curve is the 10/15/20/25 dB curve of the paper's Fig. 13
clipped at 18 dB, so its 18 dB value is also a chord. **Ours** = `results/ber_fig9_dev.mat`
(`run_ber_fig9`, 50 trials, seed 2019, a fresh channel per (SNR, trial), `E = 10`, per-edge MP
0.7 × 20, io-layer CP-OTFS with the exact ICI operator and a 20-sample CP, threshold `3σ`,
`σ² = 1`, data amplitude `√SNR_d`, `|x_p|² = SNR_p`). **impl48** = the earlier independent
reproduction `impl/otfs_ber` (48 trials, seed 2019, Raviteja's DD-circular twist model without
CP, one channel per trial reused for all SNR); **seed 2020** = the same code re-run with seed 2020.
All reproductions use Rayleigh gains `CN(0, P_l)` with `ΣP_l = 1` and the snapped EVA taps
`[0 0 1 2 3 5 8 13 19]`; the paper does not state its gain distribution or normalization, its MP
iteration count/damping, its number of frames, or the pilot position.

Error counts behind our points: 6.55 M bits per point (ideal, OFDM), 6.03 M (CE curves). The
95 % intervals are frame-level (mean ± 1.96 σ/√50 over the 50 per-frame BERs): errors cluster
in frames (one frame carries 22–58 % of the errors at 14–18 dB), so pooled-Poisson intervals
would be 2–80× too narrow.

## Side-by-side (BER; ours/paper in the last column)

OTFS, known channel

| SNR_d | paper | ours (errors / frames with errors) | 95 % interval | impl48 | seed 2020 | ours/paper |
|---|---|---|---|---|---|---|
| 10 | 1.51e-2 | 1.38e-2 (90485 / 50) | [8.9e-3, 1.9e-2] | 1.78e-2 | 1.11e-2 | 0.91 |
| 12\* | 2.96e-3 | 4.20e-3 (27539 / 46) | [1.2e-3, 7.2e-3] | 3.18e-3 | 2.82e-3 | 1.42 |
| 14\* | 5.81e-4 | 3.65e-4 (2395 / 37) | [7.6e-5, 6.6e-4] | 9.69e-4 | 5.46e-4 | 0.63 |
| 16\* | 9.99e-5 | 2.66e-5 (174 / 18) | [6e-6, 4.7e-5] | 2.89e-4 | 1.50e-4 | 0.27 |
| 18 | 1.51e-5 | 2.14e-6 (14 / 5, below the 20×-floor rule) | [0, 4.5e-6] | 1.24e-4 | 1.42e-4 | 0.14 |

OTFS CE, `SNR_p = 40 dB`

| SNR_d | paper | ours | 95 % interval | impl48 | seed 2020 | ours/paper |
|---|---|---|---|---|---|---|
| 10 | 2.37e-2 | 2.67e-2 (160676 / 50) | [1.6e-2, 3.7e-2] | 3.22e-2 | 2.11e-2 | 1.13 |
| 12\* | 9.07e-3 | 9.81e-3 (59150 / 49) | [5.2e-3, 1.4e-2] | 1.61e-2 | 7.12e-3 | 1.08 |
| 14\* | 3.47e-3 | 2.33e-3 (14038 / 47) | [8.8e-4, 3.8e-3] | 1.11e-2 | 2.32e-3 | 0.67 |
| 16\* | 1.35e-3 | 5.57e-4 (3356 / 39) | [2.0e-4, 9.1e-4] | 4.75e-3 | 5.26e-4 | 0.41 |
| 18 | 5.34e-4 | 1.50e-4 (902 / 25) | [3.5e-5, 2.7e-4] | 3.37e-3 | 3.58e-4 | 0.28 |

OTFS CE, `SNR_p = 45 dB`

| SNR_d | paper | ours | 95 % interval | impl48 | seed 2020 | ours/paper |
|---|---|---|---|---|---|---|
| 10 | 1.63e-2 | 1.62e-2 (97635 / 50) | [1.1e-2, 2.1e-2] | 2.51e-2 | 1.51e-2 | 0.99 |
| 12\* | 4.75e-3 | 5.55e-3 (33477 / 46) | [2.2e-3, 8.9e-3] | 1.11e-2 | 3.97e-3 | 1.17 |
| 14\* | 1.38e-3 | 9.53e-4 (5746 / 44) | [2.4e-4, 1.7e-3] | 1.98e-3 | 8.57e-4 | 0.69 |
| 16\* | 3.28e-4 | 1.41e-4 (851 / 27) | [3.4e-5, 2.5e-4] | 9.46e-4 | 1.25e-4 | 0.43 |
| 18 | 6.36e-5 | 1.69e-5 (102 / 11) | [1.4e-6, 3.3e-5] | 3.42e-4 | 2.73e-5 | 0.27 |

OTFS CE, `SNR_p = 50 dB`

| SNR_d | paper | ours | 95 % interval | impl48 | seed 2020 | ours/paper |
|---|---|---|---|---|---|---|
| 10 | 1.34e-2 | 1.34e-2 (80718 / 50) | [8.7e-3, 1.8e-2] | 1.85e-2 | 1.08e-2 | 1.00 |
| 12\* | 3.04e-3 | 4.05e-3 (24385 / 46) | [1.2e-3, 6.9e-3] | 3.28e-3 | 2.78e-3 | 1.33 |
| 14\* | 6.94e-4 | 4.80e-4 (2894 / 39) | [1.0e-4, 8.6e-4] | 1.88e-3 | 4.37e-4 | 0.69 |
| 16\* | 1.41e-4 | 5.62e-5 (339 / 17) | [4e-6, 1.1e-4] | 4.44e-4 | 4.84e-5 | 0.40 |
| 18 | 2.54e-5 | 3.98e-6 (24 / 8, 1.2× the 20×-floor rule) | [0, 8.6e-6] | 9.04e-5 | 1.04e-5 | 0.16 |

OFDM, known channel

| SNR_d | paper | ours | 95 % interval | impl48 | seed 2020 | ours/paper |
|---|---|---|---|---|---|---|
| 10 | 4.21e-2 | 4.18e-2 (273934 / 50) | [3.7e-2, 4.7e-2] | 4.49e-2 | 4.10e-2 | 0.99 |
| 12\* | 2.80e-2 | 2.94e-2 | [2.6e-2, 3.3e-2] | 2.99e-2 | 2.70e-2 | 1.05 |
| 14\* | 1.87e-2 | 2.00e-2 | [1.7e-2, 2.3e-2] | 1.96e-2 | 1.77e-2 | 1.07 |
| 16\* | 1.19e-2 | 1.21e-2 | [1.0e-2, 1.4e-2] | 1.27e-2 | 1.13e-2 | 1.02 |
| 18\* | 7.24e-3 | 8.03e-3 (52642 / 50) | [6.7e-3, 9.4e-3] | 8.10e-3 | 7.26e-3 | 1.11 |

At the paper's exact 15 dB vertex (ours as the geometric mean of the 14 and 16 dB points):
ideal 2.57e-4 vs 0.99e-4 (0.38), CE40 2.15e-3 vs 1.14e-3 (0.53), CE45 7.45e-4 vs 3.67e-4
(0.49), CE50 3.31e-4 vs 1.64e-4 (0.50), OFDM 1.53e-2 vs 1.56e-2 (1.02). Expressed horizontally
(SNR at which our curve reaches the paper's vertex BER, positive = ours better): ideal
+0.15 / +0.73 / +1.55 dB at the 10 / 15 / 18 dB vertices; CE40 −0.24 / +0.89 / +1.94; CE45
+0.02 / +0.74 / +1.25; CE50 0.00 / +0.65 / +1.40; OFDM +0.04 / −0.08 / about −0.5.
The paper's value lies inside our frame-level 95 % interval for every curve at 10, 12 and 14 dB
and for OFDM at every SNR; it lies above the interval for all four OTFS curves at 16 and 18 dB.
(A 200-trial run at exactly 15 dB is reported in the last section.)

## What the paper claims and what we reproduce

- **Ordering `SNR_p` 40 > 45 > 50 dB and "increasing the pilot power improves the error
  performance"**: reproduced at every SNR (ours at 18 dB: 1.5e-4 > 1.7e-5 > 4.0e-6).
- **"At `SNR_p = 50 dB` the embedded CE attains similar performance to the known channel"**:
  reproduced with the same magnitude. Horizontal gap CE50 → ideal: paper 0.31 dB (15 dB) and
  0.55 dB (18 dB); ours 0.22 dB (14), 0.57 dB (16), 0.49 dB (18). The paper's low-SNR crossover
  (CE50 below ideal at 10 dB, with CE50 < ideal < CE45) is reproduced with the same ordering
  (1.339e-2 < 1.381e-2 < 1.620e-2 vs 1.335e-2 < 1.512e-2 < 1.634e-2), but it is at the noise
  level: per frame, CE50 beats ideal in 26 of 50 frames.
- **Pilot-power penalties**: CE45 vs ideal 1.30 dB (15 dB) / 1.52 dB (18 dB) in the paper vs
  1.27 / 1.64 dB here; CE40 2.61 / 3.90 dB vs 2.34 / 3.32 dB. The CE40 curve is the shallowest OTFS
  curve in both (0.21 decades/dB in the paper, 0.27–0.30 here). Mechanism visible in our per-frame
  data: the `3/x_p = 0.030` threshold at 40 dB misses the weakest EVA path (`P = 0.0049`) with
  probability 0.17 per frame and the −12 dB path with 0.057 (0.017 / 0.0095 at 45 / 50 dB).
- **"OTFS with channel estimation significantly outperforms OFDM with known channel"**:
  reproduced; OFDM/CE40 at 18 dB is 13.6× in the paper and 54× here. The OFDM curve itself
  matches the paper within 1–11 % at every SNR (slope 0.09 vs 0.10 decades/dB), although the
  paper's OFDM reference is the 30 km/h ideal-CSI ICI-aware MP receiver of its Fig. 13 (an
  inference from identical vector data) and ours is a 120 km/h genie banded MMSE with `W = 8`.
- Not testable here: the integer-vs-fractional pilot-power statement, reduced guards, 16-QAM,
  rectangular pulses.

## Where the levels differ and why

At 10 dB every curve agrees within 13 % (0.25 dB); OFDM agrees everywhere. The four OTFS curves
diverge progressively: 1.9–2.0× (CE) to 2.6× (ideal) below the paper at its 15 dB vertex and
3.6–7× below at 18 dB, i.e. our OTFS waterfalls are steeper (ideal 0.48 vs 0.38 decades/dB over
10–18 dB; CE50 0.44 vs 0.34). Ranked by evidence:

1. **Rare tap-truncation floor frames dominate every 16–18 dB point at ~50 trials.** The known-
   channel detector uses the exact channel truncated to `|q| ≤ E = 10`; on ~2 % of EVA draws the
   truncation mismatch alone gives a per-frame BER ≈ 7e-3 independent of SNR (paired test on the
   seed-2020 draws: trial 16 has 6.7e-3 at 14, 16 and 18 dB in the impl chain, 1.4–2.0e-2 in
   ours, and is essentially all of that run's 16 and 18 dB errors; the four highest-power draws
   are error-free, so total power does not identify these frames). Our 50 seed-2019 draws contain
   none (largest per-frame ideal BER 4.6e-4 at 16 dB, 5.3e-5 at 18 dB), which happens with
   probability ≈ 0.35 per SNR point. One such frame would move our 18 dB ideal point from 2.1e-6
   to 1.3e-4 and the 16 dB point from 2.7e-5 to 1.6e-4, the exact spread between the three
   reproductions (2.1e-6 / 1.2e-4 / 1.4e-4 at 18 dB). The paper's 1.5e-5 lies between those two
   outcomes, so the 16–18 dB ideal points discriminate nothing at 50 trials. Note that a 2 %
   rate of 7e-3 frames gives an expected mean of 1.4e-4, 9× the paper's value, and that the
   paper's CE50 curve lies *above* its ideal curve at 18 dB whereas in both impl runs the CE curve
   (whose detector uses all thresholded taps, no `E` truncation) falls below the ideal curve; the
   most plausible reading is that the paper generated and detected the known-channel case with the
   same `q ∈ [−10, 10]` model (removing the mismatch), which cannot be verified.
2. **Unstated simulation details** (MP iterations/damping, frame count, pilot position, the
   truncation used when generating the known channel) can still shift the 16–18 dB tail by
   factors of this size without changing the 10 dB level; but the 200-trial exact-vertex check at
   15 dB (last section) removed the apparent "0.7 dB steeper waterfall" of the 50-trial run, so
   there is no evidence left for a systematic slope difference below 16 dB. The paper's own
   figures show a 40 % inconsistency at 15 dB (Fig. 10 at `3σ_p` gives 4.63e-4 where Fig. 9 gives
   3.31e-4 at identical settings) and reuse the same data across Figs. 9/10/11 at 18 dB.
3. **Channel gains** (unstated in the paper; `CN(0, P_l)`, `ΣP_l = 1` here): a different total
   power would be a pure horizontal shift, bounded to ≈ 0.25 dB by the 10 dB agreement and to
   ≈ 0.4 dB by the OFDM agreement, i.e. at most 1.3–1.4× at our tail slope. No evidence for it.
4. **Channel application model** (measured, wrong direction): the paper and impl apply the
   DD-circular twist model without CP; we apply the exact TF io layer with per-symbol CP. On
   identical draws the MP detectors are bit-identical and our chain is 1.15–1.5× *worse* per
   typical frame, so switching models would move us further below the paper.
5. **OFDM receiver**: all three variants agree within 12 % at every SNR; immaterial.
6. **Read-off**: none for the vertices; the 12/14/16 dB chords lie slightly below the true
   (concave) curve, so the true ratios at 14 and 16 dB are somewhat closer to 1 than tabulated.

## Verdict

At 10–15 dB every curve matches the paper within its frame-level 95 % interval (the 15 dB vertex
with 200 draws: ratios 0.84–1.56, see the last section); the OFDM reference matches at all SNRs
(1–11 %); the curve ordering, the `SNR_p = 50 dB` ≈ ideal statement with its 0.3–0.6 dB gap and
low-SNR crossover, the 45/40 dB pilot penalties and the shallow CE40 slope, and the > 1 decade
OTFS-over-OFDM margin are all reproduced. What cannot be tested numerically with a few hundred
frames are the OTFS levels at 16–18 dB, where the 50-trial points are 2.3–7× below the paper:
the population mean there is set by a ~2 %-rate class of tap-truncation floor frames (per-frame
BER ~1e-2 independent of SNR), one occurrence of which changes a 50-trial point by 6–60× and
brackets the paper's value. A decisive tail comparison needs of the order of 1000 frames per
point, or a driver that keeps the channel fixed across SNR so the floor frames can be identified.

## Exact-vertex check at 15 dB with 200 trials

`run_ber_fig9('tier', 200, 'snr_d_db', 15)` → `results/ber_fig9_n200.mat` (200 fresh channel
draws, 26.2 M bits per point for ideal/OFDM, 24.1 M for the CE curves), against the paper's exact
15 dB vertices:

| curve | paper | ours | errors / frames with errors | frame-level 95 % interval | ours/paper | without the single worst frame |
|---|---|---|---|---|---|---|
| OTFS ideal | 2.57e-4 | 4.00e-4 | 10498 / 99 | [6.5e-5, 7.4e-4] | 1.56 | 2.52e-4 (0.98; that frame = 37 % of the errors) |
| CE 40 dB | 2.15e-3 | 1.81e-3 | 43647 / 174 | [1.0e-3, 2.6e-3] | 0.84 | 1.59e-3 (0.74) |
| CE 45 dB | 7.45e-4 | 7.79e-4 | 18785 / 125 | [3.0e-4, 1.3e-3] | 1.05 | 5.94e-4 (0.80) |
| CE 50 dB | 3.31e-4 | 4.80e-4 | 11570 / 101 | [9.4e-5, 8.7e-4] | 1.45 | 3.08e-4 (0.93; 36 %) |
| OFDM ideal | 1.53e-2 | 1.50e-2 | 393257 / 200 | [1.40e-2, 1.60e-2] | 0.98 | 1.48e-2 (0.97) |

Every paper vertex lies inside our frame-level interval, and the ratios are 0.84–1.56 with the
worst single frame included (a 3 % frame of the floor class described above) and 0.74–0.98
without it, bracketing the paper. The "2× below at 15 dB" of the 50-trial run was the same
frame lottery in the other direction (those 50 draws contained no floor frame). So at 15 dB the
paper is reproduced within the draw-to-draw scatter on all five curves; only the 16–18 dB
vertices remain beyond the reach of a few hundred frames.

# What the MCRB / LB results say

Conclusions from the Phase-3 runs, cross-read against the genie CRB (F1.2, `liu_comparison.md`)
and the BER results (F2.1/F2.2, `raviteja_comparison.md`). Written 2026-09-10.

**Definitions.** `MCRB0 = A⁻¹BA⁻¹` evaluated at the numerically solved pseudo-true parameter
`η₀`, which splits as `N0·G + F` (`G = A⁻¹J̃A⁻¹` the noise term, `F = A⁻¹bbᵀA⁻¹` the score-bias
floor, both built at `σ² = 1`). `BB = (η̄−η₀)(η̄−η₀)ᵀ` is the pseudo-true bias. `LB = MCRB0 + BB`.
`CRB` is the genie bound (`X̂ = X`). Traces are variances; RMSE is their square root, so a trace
ratio of 4 is a factor 2 in RMSE.

**Data.** Real detector: `results/mcrb_rx_eva_reduced_n100_n100.mat` (EVA-reduced, 3 paths,
Rayleigh gains, OFDM dftdd CE + one-tap MMSE, OTFS embedded pilot + MP, **100 trials**,
SNR 0:2:24 dB) and the legacy 50-trial regression. Abstract errors:
`mcrb_bsc_eva_reduced_smoke.mat` (20 trials, `p = 0:0.05:0.5`) and the MT40 files (20 draws).
All statistics below are **medians of per-trial ratios** over valid trials; with Rayleigh gains
the trial mean is unusable (mean/median of the per-trial CRB reaches 500 at some SNRs).

---

## 1. The score-bias floor F is zero at the pseudo-true point, so `LB = N0·G + BB`

`η₀` is by definition the least-squares projection of the true mean onto the mismatched model, and
its stationarity condition is exactly `b = Re{G̃ᴴΔ} = 0`. Measured: `max F/(N0·G)` over all trials
is 2.3e-1 in one trial of the 100-trial run and below 1e-11 in the medians (typical values 1e-15 to
1e-23), 1.3e-9 in the legacy run, ≤ 7e-6 for `p ≤ 0.3` in the BSC sweep, ≤ 3e-9 on MT40. So `F` is
a solver-residual witness, not a term. Three consequences used throughout:

- `MCRB0` is exactly proportional to `N0`, so every `MCRB0` ratio is SNR-independent.
- All SNR-floor behaviour of `LB` comes from `BB`.
- `F` is a genuine term **only** in the `MCRB@η̄` arm (Section 5), where it makes up 70–99 % of it.

## 2. How the error rate affects MCRB0: the exact mean-field shrink law

Detection errors shrink the effective symbol, `ρ = E[X*X̂]/E|X|²`, and the noise term inflates as
`1/ρ²`, independently of where the errors fall. Measured against the BSC prediction
`ρ = 1−2p` (per-trial medians, 20 trials, CP-OFDM, delay; the other seven waveform/ICI/block
combinations agree within 1 %):

| p | 0.05 | 0.10 | 0.15 | 0.20 | 0.25 | 0.30 | 0.35 | 0.40 | 0.45 |
|---|---|---|---|---|---|---|---|---|---|
| measured MCRB0/CRB | 1.24 | 1.57 | 2.04 | 2.79 | 3.96 | 6.19 | 10.8 | 24.4 | 105 |
| `1/(1−2p)²` | 1.235 | 1.562 | 2.041 | 2.778 | 4.00 | 6.25 | 11.1 | 25.0 | 100 |

The law holds to 1–3 % up to `p = 0.4`. In RMSE this is a 11 % penalty at `p = 0.05`, 25 % at 0.1,
68 % at 0.2. MT40 extends it to constellation-specific shrinks (16-PSK Gray-bit flips vs uniform
16-QAM symbol errors) within 8 %. The linearised "excess ≈ 4p" is 17–20 % low already at `p = 0.05`
and 2.2× low at `p = 0.2`; use the full law.

## 3. How the error rate affects BB: super-linear, noise-independent, and it dominates

`BB` grows faster than the error rate and does not depend on `N0`. Log-log regression of the
per-trial `BB` on the per-trial SER, pooled over SNR (100-trial run):

| | slope | n | Spearman |
|---|---|---|---|
| OFDM delay | 1.94 | 1251 | 0.90 |
| OFDM Doppler | 1.88 | 1251 | 0.91 |
| OTFS delay | 1.31 | 904 | 0.87 |
| OTFS Doppler | 1.64 | 904 | 0.93 |

The exponent is between 1.3 and 1.9 and is not pinned down better than ±0.4 (the legacy run gives
2.4 for OFDM and 1.3 for OTFS). Because `N0·G` falls with SNR while `BB` does not, at a **fixed**
error rate `LB` flattens into a bias floor: in the BSC sweep at `p = 0.05` the ratio `BB/(N0·G)`
is 0.14 at 0 dB, 1.4 at 10 dB and 14 at 20 dB, and on MT40 at 20 dB `BB` is ≥ 98 % of `LB`
already at `p = 0.025` (its 16-PSK `LB` improves only 1.25× between 0 and 20 dB while the CRB
improves 100×). With a **real** detector the error rate itself falls with SNR, so `BB` collapses
too and `LB` returns to the CRB (Section 4).

## 4. Real detector: `LB → CRB` once SER drops below about 1 %

100 trials, per-trial medians, delay block:

| SNR | OFDM SER | OFDM LB/CRB | BB/MCRB0 | OTFS SER | OTFS LB/CRB | BB/MCRB0 |
|---|---|---|---|---|---|---|
| 0 | 4.0e-1 | 8.71 | 3.09 | 3.8e-1 | 8.85 | 1.46 |
| 4 | 2.6e-1 | 4.71 | 2.36 | 2.5e-1 | 9.55 | 3.63 |
| 8 | 1.2e-1 | 1.74 | 0.65 | 9.4e-2 | 3.22 | 1.49 |
| 12 | 6.2e-2 | 1.40 | 0.35 | 3.2e-2 | 1.13 | 0.11 |
| 16 | 2.7e-2 | 1.105 | 0.10 | 5.2e-3 | 1.001 | 0.000 |
| 20 | 9.4e-3 | 1.028 | 0.025 | 7.3e-5 | 1.000 | 0 |
| 24 | 4.9e-3 | 1.015 | 0.015 | 8.4e-6 | 1.000 | 0 |

`LB/CRB < 1.10` is first reached at 18 dB for OFDM (SER 1.6e-2) and 14 dB for OTFS (SER 5.8e-3);
`LB/CRB < 1.02` at 22 dB (SER 5.8e-3) and 14 dB respectively. The rule of thumb is
**SER ≈ 1e-2 for a 10 % penalty and ≈ 1e-3 for 2 %**, on both waveforms and in both blocks. The
decomposition of `LB` follows the same threshold: the CRB share `fCRB` goes 0.115 → 0.71 → 0.97
for OFDM at 0/12/20 dB, and the misspecification excess of the noise term is never more than 27 %
of `LB` (OTFS at 0 dB) and below 3 % from 8 dB on.

## 5. Which factor dominates: `BB`, except at low SNR; and `MCRB@η̄` is not a substitute

`BB` is 76 % of `LB` at 0 dB, 30 % at 10 dB, 2.4 % at 20 dB (OFDM, real detector), and ≥ 93 % at
every `p` in the BSC sweep at 20 dB. So:

- **At high SNR with a real detector** everything collapses onto the CRB and the bound carries no
  information beyond it.
- **At fixed error rate** (a BSC, or a receiver whose SER has floored) `BB` is the whole story.
- **At low SNR** the shrink term `N0·G` still carries 40–90 % of `LB`, so both terms are needed.

Evaluating the bound at the true parameter instead of at `η₀` is not a shortcut. `MCRB@η̄ / LB`
(10 dB, CP-OFDM, delay) is 0.71 and 0.82 at `p = 0.05, 0.1`, then 1.49, 3.98, 7.06, 6.11, 4.21,
5.12 at `p = 0.15…0.45`: it **understates** by 20–30 % at small `p` and **overstates** by 4–10×
beyond `p = 0.15`, and 70–99 % of it is the `F` term that vanishes at `η₀`. Report it only as the
legacy reference.

## 6. OFDM vs OTFS: equal in the CRB, equal in the noise term, different only through `BB`

The genie CRBs coincide (ratio 1.000 ± 0.007 on shared draws, consistent with `liu_comparison.md`:
no genie gap on frame-sized grids). At **equal error rate** the noise terms also coincide: the
paired OTFS/OFDM `G` ratio in the BSC sweep is 1.00–1.05 for `p ≤ 0.45`, and the `BB` ratio
scatters 0.52–1.6 with no trend. So any waveform difference in `LB` comes from the *error pattern*
its receiver produces, not from the waveform's own information geometry.

With the real receivers (paired, 100 trials):

| SNR | delay LB, OTFS/OFDM | trials with OTFS worse | Doppler LB, OTFS/OFDM | trials with OTFS worse |
|---|---|---|---|---|
| 0 | 0.83 | 43 % | 17.7 | 94 % |
| 4 | 1.95 | 68 % | 33.7 | 98 % |
| 8 | 1.32 | 63 % | 15.6 | 95 % |
| 12 | 0.90 | 33 % | 1.18 | 65 % |
| 16 | 0.93 | 21 % | 0.94 | 22 % |
| 24 | 0.99 | 21 % | 0.99 | 17 % |

- **Delay**: OTFS is 1.3–2× *worse* at 4–8 dB despite an equal or lower BER, because an MP failure
  corrupts a whole DD-spread symbol while a one-tap MMSE error is local; it is 5–10 % better in
  variance (2–5 % in RMSE) at 12–20 dB, where OFDM still has a 3–10 % bias share to lose. Neither
  ordering is resolved per trial (fractions 21–68 %).
- **Doppler**: OTFS is 15–35× worse below 10 dB, and this is entirely `BB` (94–98 % of OTFS's
  Doppler `LB` there, against `MCRB0` ratios of only 1.1–1.5). MP failures bias `ν` systematically.
  Crossover at 12 dB; above it the two are equal to 1 %.

**So the BER advantage of OTFS buys almost nothing in the bound.** At 16–20 dB, where OTFS's BER
is 5–140× lower, OFDM's own bias share is only 2–10 % of its `LB` — that is the entire prize. And
in the window where OTFS's BER advantage is largest in relative terms (4–10 dB), its errors are
the more damaging kind, so it loses on Doppler and on delay. This is the sharpest cross-result of
the project: **BER and sensing accuracy are not the same currency.**

## 7. Error placement matters as much as error rate

The BSC-matched control arm re-scatters the *same number* of symbol errors i.i.d. over the data
cells. Ratio of the real to the BSC-matched `LB` (delay trace, per-trial medians, 100 trials):

| SNR | 0 | 4 | 8 | 12 | 16 | 20 | 24 |
|---|---|---|---|---|---|---|---|
| OFDM | 1.01 | 0.96 | 0.73 | 0.71 | 0.58 | 0.61 | 0.48 |
| OTFS | 1.24 | 2.28 | 1.41 | 1.00 | 1.00 | 1.00 | 1.00 |

OFDM's real errors are **half as damaging** as i.i.d. errors of the same count above 8 dB: the
one-tap MMSE errs where its post-equalisation SNR is lowest, i.e. in deep fades, and those cells
carry little channel energy, so they pull `η₀` less. OTFS is the opposite below 10 dB (up to 2.3×
in delay and 12–26× in Doppler): MP failures are structured and Doppler-biasing. Conclusion: a BSC
at the measured BER predicts the **noise term** accurately (Section 2) but mis-predicts `BB` by up
to an order of magnitude in either direction.

## 8. Data-aided vs pilot-only sensing

A 20 % pilot comb has exactly `1/0.2 = 5×` the full-grid CRB. The data-aided `LB` falls below it
(delay) from 4 dB for OFDM and 8 dB for OTFS, and saturates at the pilot fraction 0.203–0.209 once
mismatch is gone. On MT40 at 20 dB the crossover is at an error rate of 0.5–1 % (16-PSK) and
1–2.5 % (16-QAM); at 0 dB the data-aided bound stays below the pilot-only CRB for every `p ≤ 0.1`.
So spending the whole frame on data and tolerating detection errors beats reserving 20 % for
pilots, provided the SER is under roughly 1 %.

## 9. MT40: Liu's OFDM advantage survives in `MCRB0`, but not visibly in `LB`

Because the shrink `ρ` is waveform-independent and the flips are paired across waveforms, `MCRB0`
inherits the `p = 0` CRB ordering exactly: the paired OFDM/OTFS `MCRB0` ratio is 0.979, 0.980,
0.966, 0.925, 0.940 across `p = 0…0.1` (16-PSK, 2-D), i.e. the same ~2 % (0.08 dB) advantage
measured at `p = 0`, held at every `p`. The `LB` ratio, however, scatters 0.75–1.13 over the same
grid, because `BB` is a per-draw quantity with heavy tails and 20 draws cannot resolve a few
percent. At equal `p`, 16-PSK gives the *larger* bound (its Gray-bit BSC injects ~4× the symbol
errors); at equal symbol-error rate the ordering reverses.

## 10. Practical limits of these numbers

- The mean-field law is reliable to `p = 0.4`, usable in the median at 0.45, and meaningless at
  0.5: the conditioning `κ` is flat at 90–160 for `p ≤ 0.4`, 466 at 0.45 and 2.7e4 at 0.5, where
  `A(η₀)` collapses onto its own zero-mean fluctuation (the Ω_N wall of `kappa_validity.tex`).
- With Rayleigh gains everything is heavy-tailed. Quote medians of per-trial ratios; the trial
  mean of `J⁻¹` has no finite expectation (`tr J⁻¹ ∝ 1/|α|²`).
- What 100 trials support: the shrink law, the `F = 0` structure, `BB` dominance at high SNR, the
  SER thresholds, the OTFS Doppler penalty (94–98 % of trials), the real-vs-BSC direction, the
  pilot crossovers. What they do not: the sign of the delay `LB` ordering between waveforms at any
  SNR, the `BB` exponent beyond ±0.4, MT40 orderings per draw, anything at `p = 0.5`.

## What would strengthen this

1. A finer BSC grid `p = 0:0.002:0.1` to see the small-`p` regime and the 20 dB `BB = N0·G`
   crossing directly instead of by extrapolation.
2. Dev tier (200–500 trials) on the full 9-path EVA, which is the gated bound channel and has more
   Doppler diversity than the reduced one used here.
3. The `nearest` mode of `bsc_match` alongside `uniform`, to separate error *placement* from error
   *energy* in Section 7.
4. More MT40 geometries (one placement so far) and ≥ 200 draws, to make the `LB` ordering a
   per-draw statement.

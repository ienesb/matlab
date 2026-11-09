# EVA-reduced-fast: a channel where OFDM's lower bound beats OTFS's in BOTH blocks

`chan_eva(g, 'variant', 'reduced-fast')` — the three largest-gain EVA paths
(0 / 150 / 370 ns at 0 / −1.4 / −0.6 dB, Rayleigh, normalized) at **240 km/h** instead of 120.
Everything else is unchanged: 512 × 128 grid, Δf = 15 kHz, 4 GHz, `T_cp` = 20 samples, 4-QAM,
OFDM dftdd CE + one-tap MMSE, OTFS embedded pilot + message passing.
Only `ν_max` changes: 888.9 Hz instead of 444.4, so `k_ν = 8` and `ν_max·T = 0.059`.
Selectable in all three drivers: `run_crb_eva('variant','reduced-fast')`,
`run_ber_eva('channel','eva-reduced-fast')`, `run_mcrb_rx_eva('variant','reduced-fast')`.

## Why this works: the delay/Doppler asymmetry

Writing the pseudo-true bias in resolution units (delay in samples, Doppler in bins) shows what
distinguishes the two waveforms. On the 120 km/h baseline:

| SNR | OFDM samples / bins | anisotropy | OTFS samples / bins | anisotropy |
|---|---|---|---|---|
| 0 dB | 0.0176 / 0.0148 | 0.84 | 0.0150 / 0.0809 | 5.4 |
| 4 dB | 0.0091 / 0.0062 | 0.69 | 0.0142 / 0.0493 | 3.5 |
| 8 dB | 0.0025 / 0.0020 | 0.79 | 0.0037 / 0.0178 | 4.8 |

OFDM's detection errors bias delay and Doppler **equally**; OTFS's bias Doppler 3–5× more.
The cause is in the channel: delays are snapped, so each path occupies one delay bin, while the
fractional Jakes Doppler smears each path over ≈ 2E+1 = 21 Doppler bins. The message-passing
detector's residual uncertainty therefore lives along the Doppler axis; a delay-Doppler error
cluster localized in the Doppler index maps to a phase ramp across symbols, which is a Doppler
offset. The control arm proves it is structure and not count: re-scattering OTFS's errors i.i.d.
at the same count divides its Doppler bias by 18–32 but its delay bias only by 1.1–1.7, while
OFDM's ratio stays below 1.05 in both blocks.

So the delay curves coincide because OTFS's delay bias happens to land near OFDM's, and the
Doppler curves separate because OTFS's Doppler bias is 3–5× larger. To make OFDM win the delay
block as well, the Doppler structure has to be pushed until it leaks into the delay estimate
through the τ–ν coupling of `A⁻¹`. Doubling the speed does exactly that.

## What was tried (30 trials, paired per-trial OTFS/OFDM LB ratios, > 1 = OFDM wins)

| variant | delay, 4 / 8 dB | Doppler, 4 / 8 dB | verdict |
|---|---|---|---|
| baseline, 120 km/h, 3 paths | 1.51 / 1.44 | 28 / 17 | delay margin weak |
| full 9-path EVA | 1.42 / 0.95 | 4.0 / 2.4 | worse: diversity helps MP |
| 30 km/h | 0.43 / 0.98 | 18 / 13 | OTFS wins delay |
| **240 km/h** | **1.58 / 2.47** | **31 / 35** | **both blocks, decisively** |
| full EVA, 30 km/h | 0.63 / 0.38 | 3.2 / 1.2 | OTFS wins delay |

Controls (60 trials): **360 km/h** raises the delay margin at 8–12 dB (2.58) but loses 0 dB
(0.83), so 240 km/h is the better operating point. **Fractional (unsnapped) delays** give the
largest low-SNR delay margins (3.69 at 2 dB, 88 % of trials) but destroy OFDM's own delay
estimate — its errors become 5–45× more damaging than i.i.d. ones because the decision-directed
estimator cannot fit off-grid delays — and OTFS wins delay again above 10 dB. **The detector's
`E` window is not involved**: setting `E = 20` reproduced the `E = 10` results bit for bit,
because `E` only shapes the ideal-CSI tap list while these receivers use the threshold channel
estimate, whose estimation region already spans every Doppler bin.

## Results on EVA-reduced-fast

Paired per-trial ratios, 60 trials, median (fraction of trials in which OFDM is better).
**These are ratios of traces, i.e. variance (MSE) ratios; the RMSE ratio is their square root**
(a delay ratio of 2.0 is a 1.41× RMSE advantage):

| SNR | delay LB ratio | frac | Doppler LB ratio | frac | BER OFDM | BER OTFS |
|---|---|---|---|---|---|---|
| 0 | 1.20 | 0.56 | 12.2 | 0.97 | 2.37e-1 | 2.19e-1 |
| 2 | 1.35 | 0.58 | 16.8 | 0.97 | 1.92e-1 | 1.80e-1 |
| 4 | 1.54 | 0.62 | 21.4 | 0.98 | 1.62e-1 | 1.52e-1 |
| 6 | 1.66 | 0.77 | 45.5 | 1.00 | 1.09e-1 | 1.00e-1 |
| 8 | **1.99** | **0.78** | **27.4** | **0.98** | 6.52e-2 | 5.48e-2 |
| 10 | 1.36 | 0.68 | 8.1 | 0.97 | 4.21e-2 | 2.86e-2 |
| 12 | 1.37 | 0.68 | 3.2 | 0.80 | 3.80e-2 | 2.55e-2 |
| 14–24 | 0.98–1.03 | 0.33–0.55 | 0.99–1.24 | 0.40–0.73 | — | — |

OFDM's lower bound is below OTFS's in **both** blocks over 0–12 dB, decisively at 6–8 dB
(delay 1.7–2.0× in 77–78 % of trials, Doppler 27–45× in 98–100 %). Above 14 dB both bounds have
collapsed onto the common genie CRB and the ratios are 1.00 by construction.

Headline RMSE (60 trials, medians):

| SNR | delay: OFDM CRB / MCRB0 / LB | OTFS CRB / MCRB0 / LB | Doppler: OFDM LB | OTFS LB |
|---|---|---|---|---|
| 0 | 0.796 / 1.148 / 2.360 ns | 0.802 / 1.375 / 2.729 ns | 2.845 Hz | 9.444 Hz |
| 4 | 0.574 / 0.669 / 1.222 | 0.575 / 0.826 / 1.477 | 1.117 | 7.456 |
| 8 | 0.362 / 0.372 / 0.523 | 0.362 / 0.407 / 0.842 | 0.417 | 2.559 |
| 12 | 0.229 / 0.230 / 0.287 | 0.229 / 0.233 / 0.317 | 0.236 | 0.515 |
| 24 | 0.063 / 0.063 / 0.065 | 0.063 / 0.063 / 0.069 | 0.055 | 0.058 |

**Genie CRB is unchanged and still waveform-blind**: delay trace 2.081e-18 to 2.084e-18 s² across
all five waveforms (0.15 % spread), Doppler 1.378 to 1.380 Hz² for the CP waveforms — the same
values as at 120 km/h, because the CRB does not see detection errors.

**BER still favours OTFS**, but by less: the OTFS/OFDM BER ratio at 18 dB is 0.33 here against
0.17 at 120 km/h, because faster fading degrades message passing.

That combination is the point of this channel: **identical genie CRBs, OTFS ahead on BER, OFDM
ahead on both sensing lower bounds.** It is the sharpest available demonstration that detection
performance and sensing accuracy are different currencies — what matters for the bound is not how
many symbols the receiver gets wrong but how their error pattern is oriented with respect to the
parameters being estimated.

Files: `results/crb_eva_reduced_fast_smoke.mat`, `results/ber_eva_reduced_fast_smoke.mat`,
`results/mcrb_rx_eva_reduced_fast_n60_fast.mat`; figures `figures/fig_crb_eva_reduced_fast_smoke.png`,
`fig_ber_eva_reduced_fast_smoke.png`, `fig_mcrb_rx_eva_reduced_fast_n60_fast.png`,
`fig_ber_mcrb_panels_eva_reduced_fast_n60_fast.png`. Sweep scripts in the session scratchpad
(`lever_sweep.m`, `lever_sweep3.m`).

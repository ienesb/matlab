# EVA-reduced: BER, CRB and MCRB/LB

The channel is the existing `EVA-reduced` variant, `chan_eva(g, 'reduced', true)`: the **three
largest-gain paths of the EVA profile**, normalized to unit total power, with the same Rayleigh
gains and Jakes Doppler as full EVA.

| path | delay | rel. power | normalized `P_l` | tap `l_p` |
|---|---|---|---|---|
| 1 | 0 ns | 0 dB | 0.4053 | 0 |
| 2 | 150 ns | −1.4 dB | 0.2936 | 1 |
| 3 | 370 ns | −0.6 dB | 0.3011 | 3 |

(The full EVA powers in descending order are 0, −0.6, −1.4, −1.5, −3.6, −7.0, −9.1, −12.0,
−16.9 dB, so these are the top three; they are also well separated in delay, unlike the 0/30 ns
pair of full EVA that collides on tap 0.) Grid 512 × 128, Δf = 15 kHz, `T_cp` = 20 samples,
4 GHz, 120 km/h, Gray 4-QAM.

Runs: `run_crb_eva('tier','smoke','reduced',true,'snr_db',0:2:24)`,
`run_ber_eva('tier','smoke','channel','eva-reduced','snr_db',0:2:24)` and
`run_mcrb_rx_eva('tier',100,'snr_db',0:2:24)`. The BER and MCRB runs share seed 9100 and the
substream key `1000 + 10·SNR_dB`, so they see **identical channel draws at each SNR** (verified,
5 of 5 trials checked).

## 1. Genie CRB — all five waveforms coincide

Delay and Doppler traces at `σ² = 1`, mean over 20 trials:

| waveform | delay trace [s²] | Doppler trace [Hz²] |
|---|---|---|
| CP-OFDM | 2.0572e-18 | 1.3802 |
| CP-OTFS | 2.0575e-18 | 1.3815 |
| RCP-OFDM | 2.0617e-18 | 1.4938 |
| RCP-OTFS (sampled) | 2.0577e-18 | 1.4939 |
| RCP-OTFS (continuous) | 2.0623e-18 | 1.4953 |

Delay spread across waveforms 0.25 %, ICI vs no-ICI 0.01 %. The 8 % higher Doppler bound of the
reduced-CP waveforms is the shorter symbol period (`T` vs `T + T_cp`, ratio 1.079), not a waveform
effect. CP-OFDM with 4-QAM sits exactly on its `(E J)⁻¹` closed form (0.00 % Jensen gap, as
expected for constant-modulus data). **No OFDM/OTFS sensing advantage exists here**, consistent
with `liu_comparison.md`: 65 536 samples for 12 unknowns is far from the identifiability edge
where Liu's effect lives.

## 2. BER

20 trials, `N0 = 1/SNR`, matched data-referenced pilot budgets:

| SNR | OFDM dftdd | OFDM dft | OFDM interp | OFDM ideal 1-tap | OFDM ideal banded | OTFS CE | OTFS ideal |
|---|---|---|---|---|---|---|---|
| 0 | 2.48e-1 | 2.52e-1 | 2.89e-1 | 2.18e-1 | 2.18e-1 | 2.28e-1 | 2.14e-1 |
| 4 | 1.64e-1 | 1.71e-1 | 2.03e-1 | 1.47e-1 | 1.46e-1 | 1.48e-1 | 1.37e-1 |
| 8 | 6.15e-2 | 7.04e-2 | 9.50e-2 | 5.89e-2 | 5.83e-2 | 4.45e-2 | 3.64e-2 |
| 12 | 2.79e-2 | 3.28e-2 | 4.70e-2 | 2.70e-2 | 2.65e-2 | 1.28e-2 | 1.02e-2 |
| 16 | 1.31e-2 | 1.53e-2 | 2.68e-2 | 1.28e-2 | 1.23e-2 | 2.01e-3 | 1.38e-3 |
| 20 | 8.30e-3 | 9.84e-3 | 1.79e-2 | 8.14e-3 | 7.30e-3 | 4.53e-5 | 1.78e-5 |
| 24 | 2.90e-3 | 3.56e-3 | 1.07e-2 | 2.87e-3 | 2.39e-3 | 3.87e-7 | 7.74e-7 |

OTFS is 1 to 4 decades below OFDM above 12 dB; the estimated-channel receivers track their
ideal-CSI references within 10 % (dftdd) and 15 % (OTFS CE). OFDM's curve is diversity-1 and
flattens: with only 3 paths and 4-QAM its one-tap MMSE cannot exploit the delay diversity that MP
extracts in the delay-Doppler domain. The non-monotonic OTFS points at 14–22 dB (2.0e-3 at 16 dB
against 4.5e-5 at 20 dB) are the tap-truncation floor frames documented in
`raviteja_comparison.md`: about 2 % of draws floor at ~1e-2, so a 20-trial point at that level is
a lottery. Take those four points as ≲ 2e-3, not as a curve.

## 3. MCRB / LB from the real detector

100 trials, medians of per-trial values, RMSE in ns and Hz. `CRB` here is the genie bound on the
same draws (median; its mean is 1.465 ns at 0 dB, matching the CRB sweep's 1.434 ns — the gap
between the two statistics is the `1/|α|²` heavy tail, not a modelling difference).

**Delay [ns]**

| SNR | OFDM CRB | OFDM MCRB0 | OFDM LB | OTFS CRB | OTFS MCRB0 | OTFS LB |
|---|---|---|---|---|---|---|
| 0 | 0.809 | 1.141 | 2.773 | 0.809 | 1.391 | 2.566 |
| 4 | 0.587 | 0.679 | 1.418 | 0.589 | 0.810 | 2.056 |
| 8 | 0.379 | 0.385 | 0.547 | 0.379 | 0.415 | 0.718 |
| 12 | 0.221 | 0.222 | 0.289 | 0.221 | 0.224 | 0.252 |
| 16 | 0.151 | 0.151 | 0.168 | 0.151 | 0.151 | 0.154 |
| 20 | 0.100 | 0.100 | 0.105 | 0.100 | 0.100 | 0.100 |
| 24 | 0.060 | 0.060 | 0.061 | 0.060 | 0.060 | 0.060 |

**Doppler [Hz]**

| SNR | OFDM CRB | OFDM MCRB0 | OFDM LB | OTFS CRB | OTFS MCRB0 | OTFS LB |
|---|---|---|---|---|---|---|
| 0 | 0.698 | 0.930 | 2.064 | 0.697 | 1.181 | 9.426 |
| 4 | 0.500 | 0.570 | 0.944 | 0.499 | 0.742 | 5.841 |
| 8 | 0.313 | 0.329 | 0.422 | 0.314 | 0.359 | 2.058 |
| 12 | 0.178 | 0.177 | 0.210 | 0.177 | 0.178 | 0.243 |
| 16 | 0.129 | 0.128 | 0.135 | 0.129 | 0.129 | 0.130 |
| 20 | 0.086 | 0.086 | 0.089 | 0.086 | 0.086 | 0.086 |
| 24 | 0.051 | 0.051 | 0.052 | 0.051 | 0.051 | 0.051 |

## 4. Reading the three together

- **The genie CRB is waveform-blind, the BER is not, and the bound sits in between.** OTFS's
  1-to-4-decade BER advantage does not translate: at 16 dB its delay LB is 0.154 ns against OFDM's
  0.168 ns, an 8 % gain, and by 20 dB both are on the common CRB.
- **Mismatch is a low-SNR phenomenon on this channel.** `LB/CRB` is 3.4 (OFDM) and 3.2 (OTFS) at
  0 dB, 1.31 and 1.14 at 12 dB, and 1.05 and 1.00 at 20 dB. The crossing of the 10 % mark happens
  at SER ≈ 1e-2, i.e. 18 dB for OFDM and 14 dB for OTFS, exactly the threshold established in
  `mcrb_conclusions.md`.
- **OTFS is the worse Doppler sensor below 12 dB**, by 4.6× in RMSE at 0 dB and 4.9× at 8 dB,
  entirely through the pseudo-true bias: its MP failures displace `ν` systematically while its
  `MCRB0` is only 1.1–1.3× OFDM's. The two are equal from 12 dB on.
- **OFDM's BER floor costs it nothing in the bound above 20 dB.** Its BER stalls near 3e-3 at
  24 dB, but that SER is already low enough for `LB/CRB = 1.015`.
- 3 paths instead of 9 changes the absolute levels (delay CRB 0.81 ns vs 9.99 ns median at 0 dB on
  full EVA) and reduces the diversity that OTFS exploits, so the OTFS BER advantage here is
  smaller than on full EVA. The conclusions of `mcrb_conclusions.md` are unaffected.

Files: `results/crb_eva_reduced_smoke.mat`, `results/ber_eva_reduced_smoke.mat`,
`results/mcrb_rx_eva_reduced_n100_n100.mat`, `results/mcrb_bsc_eva_reduced_smoke.mat`;
figures `figures/fig_crb_eva_reduced_smoke.png`, `fig_ber_eva_reduced_smoke.png`,
`fig_mcrb_rx_eva_reduced_n100_n100.png`, `fig_ber_mcrb_panels_eva_reduced_n100_n100.png`.
The earlier no-Doppler BER run under the same name is kept as
`results/ber_eva_reduced_smoke_nodoppler.mat`.

# OTFS Embedded-Pilot Channel Estimation — Fig. 9 reproduction

Self-contained MATLAB reproduction of **Figure 9** of

> P. Raviteja, K. T. Phan, Y. Hong, *"Embedded Pilot-Aided Channel
> Estimation for OTFS in Delay–Doppler Channels,"* IEEE Trans. Veh.
> Technol., vol. 68, no. 5, May 2019. (the **Embedded** paper)

The OTFS input/output relation and the message-passing (MP) detector are
taken from the companion paper bundled here:

> P. Raviteja, K. T. Phan, Q. Jin, Y. Hong, E. Viterbo, *"Low-Complexity
> Iterative Detection for OTFS,"* IEEE WCNC 2018. (the **MP** paper)

## Figure 9

BER vs `SNR_d` for the **fractional-Doppler, full-guard** embedded-pilot
scheme. Five curves:

| Curve | Meaning |
|-------|---------|
| `SNR_p = 40/45/50 dB` | OTFS with the proposed threshold-based channel estimation + MP detection |
| OTFS, Ideal | OTFS with perfectly known channel + MP detection |
| OFDM, Ideal | OFDM with perfectly known (ICI) channel + MMSE — reference |

System: `N = 128` (Doppler), `M = 512` (delay), `Δf = 15 kHz`, `f_c = 4 GHz`,
4-QAM, EVA channel, UE 120 km/h (`k_ν = 4`), `l_τ = 20`, threshold `T = 3σ`.

## How it works

1. **`gen_eva_channel`** — draws one EVA realization: 9 taps with the 3GPP
   power-delay profile, integer delay taps, Rayleigh gains, and a Jakes'
   Doppler `ν_i = ν_max cos θ_i` per tap → (fractional) Doppler tap.
2. **`build_effective_dd`** — builds the exact 2D effective delay-Doppler
   response `G[k',l']` (MP eq. 13 / Embedded eq. 3). With a finite `E` it
   truncates the inter-Doppler-interference spread for the detector.
3. **`apply_dd_channel`** — applies `G` to the DD symbols as one
   Doppler-domain circular convolution per occupied delay tap, adds AWGN.
4. **`estimate_channel`** — threshold-based embedded-pilot estimation in the
   pilot delay band (Embedded eq. 9): `Ghat = y/x_p` where `|y| ≥ 3σ`.
5. **`mp_detect_otfs`** — vectorised Gaussian-approximation MP detector
   (MP eqs. 18–21), with pilot/guard symbols fixed as known nodes.
6. **`ofdm_ideal`** — OFDM reference: exact per-symbol frequency channel
   `H_n = W H_t W^H` (banded, with ICI), known-channel MMSE equalisation.

## Run

```matlab
cd impl/otfs_ber
run_sim                       % full run -> results.mat, fig_9.png
run_sim('quick', true)        % fast smoke run
run_sim('n_trials', 30)       % override any params_fig9 field
run_sim('nworkers', 8)        % parallel pool size (parfor over trials)
```

Headless: `matlab -batch "cd('impl/otfs_ber'); run_sim"`.

All parameters live in **`params_fig9.m`**. Trials run in parallel via
`parfor`, each with its own reproducible RNG substream (`p.seed`).

## Files

```
params_fig9.m        parameters (system, EVA, sweep, detector)
gen_eva_channel.m    EVA delay-Doppler channel realization
build_effective_dd.m exact / truncated effective DD response  (MP eq. 13)
apply_dd_channel.m   channel application + AWGN
estimate_channel.m   threshold-based embedded-pilot CE         (Embedded eq. 9)
dd_to_taps.m         DD response -> MP tap list
mp_detect_otfs.m     vectorised MP detector                    (MP eqs. 18-21)
ofdm_ideal.m         OFDM ideal-channel reference (MMSE)
qam4_mod.m/qam4_demod.m  4-QAM Gray (de)mapping
run_sim.m            driver (parfor over trials) -> Fig. 9
plot_fig9.m          rendering
```

## Status / validation

`fig_9.png` (48 trials) reproduces Fig. 9 well — see the BER table below
against the paper's published curves:

| SNR_d (dB) | 10 | 12 | 14 | 16 | 18 |
|---|---|---|---|---|---|
| OTFS, Ideal | 1.8e-2 | 3.2e-3 | 9.7e-4 | 2.9e-4 | 1.2e-4 |
| CE SNR_p=50 | 1.9e-2 | 3.3e-3 | 1.9e-3 | 4.4e-4 | 9.0e-5 |
| CE SNR_p=45 | 2.5e-2 | 1.1e-2 | 2.0e-3 | 9.5e-4 | 3.4e-4 |
| CE SNR_p=40 | 3.2e-2 | 1.6e-2 | 1.1e-2 | 4.8e-3 | 3.4e-3 |
| OFDM, Ideal | 4.5e-2 | 3.0e-2 | 2.0e-2 | 1.3e-2 | 8.1e-3 |

Matches the paper on all the salient features:

- correct curve ordering (`SNR_p = 40` worst, flattening toward a floor →
  `45` → `50 ≈ ideal`), with all curves clustering near 2–3e-2 at 10 dB;
- steep OTFS waterfall vs. the shallow, diversity-1 OFDM-Ideal reference
  (~4.5e-2 → ~8e-3 over 10–18 dB, close to the paper's ~3e-2 → ~6e-3);
- channel estimation at `SNR_p = 50 dB` is essentially on the ideal curve,
  the paper's central claim.

The OTFS-Ideal tail is a little steeper than the published curve at high
SNR_d (near-ML MP over the full DD diversity); exact path-gain
normalization, EVA delay-tap rounding, and the Doppler-spread truncation
`E` shift the absolute level slightly. Increase `p.n_trials` for smoother
low-BER tails.

### Notes on modelling choices

- **Doppler:** `θ_i ~ U[-π,π]` (Embedded paper, Sec. IV). The MP paper writes
  `U(0,π)`; both give the same `cos θ` distribution.
- **Threshold:** Fig. 9 uses `T = 3σ_p` on `ŷ = y/x_p`, identical to `T = 3σ`
  on `|y|` with `σ = 1`; we threshold `|y|`.
- **OFDM reference** uses MMSE with the exact known `H_n` rather than MP; at
  `ν_max/Δf ≈ 0.03` the ICI is small and both give the same diversity-1
  reference curve.

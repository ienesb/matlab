# Channel catalog — `eva/` scripts

Every channel model defined across the `eva/` scripts, verified against the code.
A "channel" here = the set of paths `(alpha_p, tau_p, nu_p)` plus how they are
drawn per trial. Error models (BSC vs realistic detector) and receiver choices
are **not** channels — see the note at the bottom. Last refreshed 2026-08-17.

## 1. EVA family (realistic multipath — the receiver/BER/MCRB workhorse)

| # | Channel | Paths / delays | Powers | Doppler | Defined in |
|---|---|---|---|---|---|
| 1 | **Reduced EVA** (primary) | 3 taps, [0, 150, 370] ns → integer taps [0,1,3] at Δf=15 kHz (or [0,2,6] at 30 kHz) | [0, −1.4, −0.6] dB, normalized Σ\|a\|²=1 | per-path ν = ν_max·cos(U), ν_max ≈ 444 Hz (v=120 km/h, f_c=4 GHz); random phase per trial | `mcrb_ce_ici.m`, `ber_ce_channels.m` ('eva'), `ber_ce_ici.m`, `compute_mcrb_mc.m` ('eva'/'evace'), `run_fim_avg.m` ('EVA'), `mcrb_avg.m`, `compute_mmle_mc.m`, `pseudotrue_diff*.m`, `ber_ofdm_ce_vs_perfect.m`, `mcrb_sweep_grid.m`, `ber_confirm_config.m`, `diag_*.m` |
| 2 | **Full EVA** | 9 taps, [0, 30, 150, 310, 370, 710, 1090, 1730, 2510] ns | [0, −1.5, −1.4, −3.6, −0.6, −9.1, −7.0, −12.0, −16.9] dB | same as reduced | `mcrb_ce_ici.m` / `ber_ce_ici.m` (`use_full_eva=true`; rarely used — worse-conditioned A) |

**Reduced-EVA variants (toggles, not separate channels):**
- **Doppler scale:** nominal (~444 Hz) / `doppler_scale=0` (ν=0 → ICI vanishes automatically) / arbitrary scale (`ber_ce_channels`, `mcrb_ce_ici`, `ber_rayleigh_eva(…,doppler_on)`).
- **ICI model:** `use_ici=false` keeps Doppler but drops intra-symbol leakage (`mcrb_ce_ici`) — measured to change bounds/BER by <0.3 % at 444 Hz.
- **Tap fading model** (`ber_rayleigh_eva(fading,…)`):
  - `'constmag'` (default everywhere else): fixed \|amp\|, random phase — per-subcarrier gain still Rayleigh-like, diversity order ≈ 0.99 (`plot_hist_gain_constmag.m`).
  - `'rayleigh'`: α_p ~ CN(0, power_p) — genuine Rayleigh taps.
- **Grid variants** (`mcrb_sweep_grid.m`): same delays on other (N, M, Δf) grids, e.g. Δf=30 kHz → B=15.4 MHz, taps [0,2,6] (the Part-15 sweet spot); the OTFS delay guard `G=l_tau+1+guard_extra` is a receiver knob, not a channel change.

## 2. Paper family (Liu et al. 2026 reproduction — CRB/MCRB)

| # | Channel | Paths | Doppler | Grid | Defined in |
|---|---|---|---|---|---|
| 3 | **paper** | 40 delay-only targets, uniform in [0, (N−1)/B] with ≥1-sample separation, unit amplitude, seed 20260514 | none (ν≡0, known) | N=128, M=1, B=160 MHz | `run_fim_avg.m`, `compute_mcrb_mc.m` / `compute_mcrb_mc_true.m` ('paper') |
| 4 | **paperv2** | same 40 delays (identical rng draw) | per-path ν ~ U[±ν_max], ν_max = 0.5/(M·T) ≈ 39 kHz (seed 20260707) | N=128, M=16, B=160 MHz, ICI on | `run_fim_avg.m`, `compute_mcrb_mc(_true).m` ('paperv2') |
| 5 | **paperdopp** | Doppler-only **dual** of 'paper': 40 targets, **τ≡0**, distinct Dopplers filling [0, (M−1)/(M·T)] with ≥1-bin separation | yes (the estimated parameter) | **N=1**, M=128, Δf=15 kHz | `run_fim_avg.m` ('paperdopp') |
| 6 | **paperv2comms** | 4 taps [0, 2, 5, 8] samples | [0, −1.4, −3.6, −6.0] dB; ν ~ U[±0.5/(M·T)] | N=128, M=16, B=160 MHz | `ber_ce_channels.m` ('paperv2comms') |

## 3. Idealized / control channels

| # | Channel | Definition | Purpose | Defined in |
|---|---|---|---|---|
| 7 | **AWGN** | single unit-gain tap, τ=0, ν=0, \|H\|≡1 (random phase) | validation: OFDM = OTFS = Q(√SNR) | `ber_awgn.m` |
| 8 | **Single-path flat Rayleigh** | 1 tap, α ~ CN(0,1) per frame, τ=0, ν=0 | control: no diversity to harvest → OFDM = OTFS on the order-1 Rayleigh curve | `ber_rayleigh_flat.m` |
| 9 | **Rayleigh-ized reduced EVA** | reduced-EVA delays/powers, but α_p ~ CN(0, power_p) | true-Rayleigh multipath: OFDM order-1, OTFS harvests delay diversity | `ber_rayleigh_eva.m` ('rayleigh') |

## 4. FIM Δf-sweep channel

| # | Channel | Definition | Defined in |
|---|---|---|---|
| 10 | **5-path random channel** (Δf-independent) | 5 paths, delays uniform in [0, 300 ns] with ≥40 ns separation (seed 20260712), unit amplitudes, ν ~ U[±ν_max]; drawn once in physical units and reused for every Δf so the sweep isolates the grid | `run_fim_sweep_deltaf_fixedB.m`, `…_fixedB_fixedNM.m`, `…_fixedN.m` |

## Not channels (easy to confuse)

- **`'evace'`** (`compute_mcrb_mc.m`): the *reduced-EVA channel* (#1) with the abstract **BSC error model** replacing the real detector — an error-model scenario, not a new channel.
- **BSC vs realistic detector**, CE method (`interp`/`dft`/`dftdd`), OTFS guard size, `Emax` DD truncation, pilot budget: receiver/error-model knobs applied *on* a channel.
- **Integer-tap rounding**: all EVA-family delays are rounded to the nearest tap of the active grid (the Raviteja assumption), so the *same* physical delays land on different taps at different Δf.

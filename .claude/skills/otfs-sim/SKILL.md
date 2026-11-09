---
name: otfs-sim
description: Use when running, debugging, or extending the OTFS/OFDM BER and CRB/MCRB simulations in crb/EVA. Covers EVA channel parameter setup, the message-passing (MP) detector, DD/TF domain conversions, and result plotting. Trigger on tasks about eva_ber_mp, mp_detector, MCRB, sfft/isfft, or reproducing Raviteja 2019 figures.
allowed-tools: Bash(matlab:*), Read, Edit, Glob, Grep
---

# OTFS/OFDM Simulation Helper

Domain knowledge for the EVA-channel OTFS vs OFDM study in `crb/EVA/`.
Reproduces results from Raviteja et al., IEEE TVT 2019.

## Core parameters (initialize_parameters.m)
- N = 512 subcarriers, M = 128 symbols, 4-QAM
- Carrier 4 GHz, subcarrier spacing 15 kHz, velocity 120 km/h (nu_max ≈ 540 Hz)
- EVA: 9 paths, delays 0–2510 ns (max delay tap 20), max Doppler tap 4
- SNR 10–18 dB in 2 dB steps, 50 Monte Carlo trials
- MP detector: 30 iterations, damping factor **0.7**

## Key files
| File | Role |
|------|------|
| `eva_ber_mp2.m` | Main BER script (OFDM CSI, OTFS CSI+MP, OTFS CE+MP) |
| `mp_detector2.m` | Message-passing iterative detector |
| `channel.m` | Multipath fading channel in TF domain |
| `sfft.m` / `isfft.m` | TF <-> DD domain transforms |
| `getMCRB_OFDM.m` / `getMCRB_TF.m` | Mismatched CRB calculations |

## Conventions that bite
- Channel runs in **time-frequency (TF)** domain; OTFS converts to/from
  **delay-Doppler (DD)** via `sfft`/`isfft`. Mixing domains is the usual bug.
- **Integer delay model**: delays are rounded to the nearest tap. Skipping
  this rounding produces wrong CE BER curves.
- Damping factor must stay at 0.7 — lowering it (e.g. 0.5) degrades the CE
  BER convergence. (See project memory on the CE BER fix.)
- Pilot is a single DD-grid impulse at `(l_p, k_p)` with a guard region.

## Common tasks
- **Run BER curves**: `matlab -batch "cd('crb/EVA'); eva_ber_mp2"`,
  results in `eva_ber_mp_results.mat`. Uses `parfor`.
- **Plot MLE vs MCRB**: run `plot_mle_vs_mcrb.m`.
- **Debug an error floor**: check integer delay rounding and damping first,
  then verify SFFT round-trip on the pilot.

When editing simulation scripts, match the existing MATLAB style (1-based
indexing, vectorized ops, no unasked-for refactors of the detector loop).

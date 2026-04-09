# CLAUDE.md — OTFS/OFDM Research Project

## Project Overview
MATLAB research project comparing **OFDM vs OTFS** modulation in realistic multipath/mobility scenarios. Main focus: Cramer-Rao Bound (CRB) analysis, channel estimation, and BER performance under the EVA (Extended Vehicular A) channel model. Reproduces results from Raviteja et al. IEEE TVT 2019.

## Directory Structure
- `crb/EVA/` — **Active development area.** BER simulations, CRB/MCRB calculations, MP detector, channel estimation for EVA channel.
- `crb/ici/` — Inter-Carrier Interference specific analysis.
- `beam-space/` — Beam-space OTFS/OFDM signal processing.
- `bistatic/` — Bistatic radar/communication simulations.
- `ofdm/`, `otfs/` — Modulation-specific implementations.
- `channel/` — Channel modeling utilities.
- `ltfat/` — Large Time-Frequency Analysis Toolbox (external library).

## Key Files in `crb/EVA/`
| File | Purpose |
|------|---------|
| `initialize_parameters.m` | EVA channel parameters (N=512, M=128, fc=4GHz, deltaf=15kHz, 9-path) |
| `eva_ber_mp2.m` | Main BER script: OFDM(CSI), OTFS(CSI+MP), OTFS(CE+MP) — reproduces Fig. 9 |
| `mp_detector2.m` | Message-passing iterative detector |
| `channel.m` | Applies multipath fading channel in TF domain |
| `sfft.m` / `isfft.m` | Symplectic FFT/IFFT (TF <-> DD domain) |
| `getFIM_OFDM_ICI.m` | Fisher Information Matrix for OFDM with ICI |
| `getMCRB_OFDM.m` | Mismatched CRB for OFDM |

## Simulation Parameters
- **N** = 512 subcarriers, **M** = 128 symbols
- **Modulation**: 4-QAM
- **Carrier**: 4 GHz, subcarrier spacing 15 kHz
- **Velocity**: 120 km/h (nu_max ~ 540 Hz)
- **EVA channel**: 9 paths, delays 0–2510 ns, max delay tap 20, max Doppler tap 4
- **SNR range**: 10–18 dB (2 dB steps), 50 Monte Carlo trials
- **MP detector**: 30 iterations, damping factor 0.7

## Conventions
- Channel operates in **time-frequency (TF)** domain; OTFS uses `sfft`/`isfft` to convert to/from **delay-Doppler (DD)** domain.
- Pilot placement: single impulse at `(l_p, k_p)` in DD grid with guard region.
- Integer delay model: delays are rounded to nearest tap for simulation (matching paper's assumption).

## Research Tasks (from `crb/todo.txt`)
1. CRB functions for OFDM/OTFS (no ICI/ISI) — done
2. MCRB functions for OFDM/OTFS — done
3. Compare monostatic/genie CRB for OTFS vs OFDM — done
4. CRB functions for OFDM with ICI channel — in progress

## Running Simulations
- Open MATLAB, `cd` to `crb/EVA/`
- Run `eva_ber_mp2` for BER curves (uses `parfor`, benefits from Parallel Computing Toolbox)
- Results plotted automatically; saved to `eva_ber_mp_results.mat`

## Reference Papers
- Raviteja et al., "Embedded Pilot-Aided Channel Estimation for OTFS in Delay-Doppler Channels," IEEE TVT 2019
- Raviteja et al., "Low-complexity iterative detection for OTFS," IEEE WCNC 2018

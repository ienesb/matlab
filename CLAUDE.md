# CLAUDE.md — OTFS vs OFDM research project (MATLAB)

## What this is
MATLAB study of **OFDM vs OTFS** for joint communication and sensing (ISAC): genie CRB,
BER with realistic receivers, and misspecified CRB / total lower bound (MCRB, LB) when the
sensing estimator uses detected data symbols. Derivations live in `../latex/` (`io/` signal
models → `fim/` FIM and CRB → `mcrb/` MCRB, pseudo-true parameter, LB). Two channels, defined in
`crb/channels.md`: **EVA** (9-path 3GPP Extended Vehicular A, delay + Doppler, Raviteja et al.
IEEE TVT 2019) and **MT40** (40 static delay-only targets, Liu et al. arXiv:2605.14451).

## Where the code is
| Folder | Status |
|---|---|
| `crb/` | **Active.** Clean re-implementation of every derivation. Start here: `crb/README.md` (conventions, layout, old→new name map, status), `crb/channels.md` (channels, numerology, translation table). |
| `old/eva/` | Previous validated implementation (moved 2026-09-06). `old/eva/utils/` (63 functions) is the **regression oracle** used by `crb/test/` through `oracle_paths`; production code in `crb/` never calls it. `old/eva/EXPERIMENTS.md` is the old experiment log. |
| `old/crb/EVA/` | The original scripts (`eva_ber_mp2`, `initialize_parameters`, …). Historical only. |
| `impl/otfs_ber/` | Raviteja-faithful per-edge MP detector + Fig. 9 reproduction (48-trial reference table in its README). |
| `impl/2605.14451v1/` | Liu et al. Fig. 1 reproduction (`run_sim.m`, `verify_my_derivation.m`), source of the MT40 geometry and precoders. |
| `fim/` | Stale older copy: its ICI Doppler derivative lacks the `1/(j2π)` normalization. **Never port from it.** |
| `beam-space/`, `bistatic/`, `ofdm/`, `otfs/`, `channel/`, `isi-ici/`, `ltfat/` | Unrelated / earlier explorations and external libs. |

## Running things (`crb/`)
```
cd matlab/crb; setup_crb              % adds cfg/ io/ fim/ rx/ mcrb/ run/ plot/ test/; asserts no name collisions
run_all_tests                         % identities + regressions (~1 min)
run_all_tests('full')                 % + Monte-Carlo closed-form checks (~5 min here, 260 checks)
```
Headless (Linux, R2026b): `cd matlab/crb && /usr/local/bin/matlab -batch "setup_crb; run_all_tests"`.
All results in `crb/results/` were produced on Windows with R2025b; the workspace moved to Linux on 2026-11-09.
**Toolboxes (2026-11-09):** Optimization and Parallel Computing Toolboxes were installed on Linux with
`mpm install --release=R2026b --destination=/usr/local/MATLAB/R2026b --products Optimization_Toolbox Parallel_Computing_Toolbox`;
`run_all_tests('full')` passed 249/249 (280 s), `run_manuscript_checks` 115/115. Since 2026-11-18 the full suite has
256 checks (7 new ones for the Gauss-Newton estimator `mmle_cp(..., 'method','lsq')` and the `pseudo_true` `target`
option): 256/256; 260/260 since 2026-11-20 (4 pinned-lattice checks in `test_rx`). This PC has 4 cores, so the default `Processes` pool is 4 workers (6 or 8 workers give the same
throughput, measured 2026-11-18); `run/run_manuscript_confirmation.m` now opens the profile's default pool.
The 8-frame September benchmark reproduces digit for digit here but takes 170 s instead of 64 s (2.7x slower).
**Phase 2 runs (2026-11-18):** `run/run_phase2_eva.m` is `run_mcrb_rx_eva`'s frame loop (bit-identical legacy
outputs for the same seed, SNR and trial) plus both control arms, the estimator arms, a speed option, a genie-CSI
detector option, checkpoints after every chunk and resume; `run/run_phase2_campaign.m` is the queue (250 of the
500 planned frames per point since the trim of 2026-11-18, option `frames`). Results are
`results/CKPT_p2_*.mat` (`p2_*.mat` only once a run has all 500 frames; for 60 km/h use the file with the suffix
`_lat37x2`, made on 2026-11-20 with the OFDM pilot lattice pinned by `'ofdm_lattice', [37 2]`, not the first one; the genie-CSI detector runs of 2026-11-19 have the suffix `_geniecsi`); see `../latex/submission_plan.md` (Phase 2) for the design, the commands and the status.
All Phase 2 runs are finished (2026-11-20); `../latex/tools/export_phase2.py` turns the nine valid files into the paper's figures, tables and number macros.
**Estimator supplement (2026-11-26):** `run/run_phase2_mmle.m` reruns only the estimator arms on the frames of a finished Phase 2 run (each regenerated
frame is checked against the stored true parameter and bit-error count) and writes `results/p2mmle_*.mat`; the two files for 120 and 240 km/h hold the
estimator at 2, 6, 10 and 24 dB (250 frames, 4.7 h for both), and the exporter reads them together with the campaign files. Log: `results/phase2_mmle.log`.
Drivers in `crb/run/` take `'tier', 'smoke'|'dev'|'paper'` (or a trial count) and write
`crb/results/*.mat` + `crb/figures/*.png`; every driver prints its acceptance table (PASS/WARN).
Parallel Computing Toolbox: drivers use `parfor` (8 workers on the old PC, 4 here). Optimization Toolbox: `lsqnonlin`
in `pseudo_true` (built-in LM fallback exists).

## Conventions that bite (asserted by `crb/test/`)
- **Grid letters**: `N` = subcarriers = delay axis = rows, `M` = OFDM symbols = Doppler axis =
  columns; TF `X` and DD `x` are both `N×M`. **Raviteja's paper swaps the letters** (his
  `M = 512` subcarriers, `N = 128` symbols). EVA: `N = 512`, `M = 128`.
- **Three timing constants**: `T = 1/Δf`, sample `Ts = 1/(NΔf)`, symbol `TS = T + Tcp` (CP) or
  `T` (reduced CP). Never alias `Ts` and `TS`. EVA `Tcp = 20 Ts = 2.604 µs` → Doppler bin 112.8 Hz.
- EVA `ν_max = 444.44 Hz` (120 km/h, 4 GHz), `l_τ = 20`, `k_ν = 4`; the 0 ns and 30 ns paths
  collide on tap 0. Gains Rayleigh `CN(0, P_l)` per trial (powers normalized to sum 1).
- Parameter vector `η_p = [α_R; α_I; τ; ν]` stacked `4P`; inversions in **bin units** (`eta_scale`);
  the `ν` block is masked explicitly when `ν ≡ 0` (MT40).
- ICI kernel derivatives are **reduced** (`I' = (1/j2π) ∂I/∂ν`), limits `1, T/2, T²/3`.
- `use_ici = false` (drop Doppler leakage) ≠ `ν = 0` (derivative kernels survive) ≠
  `doppler_scale = 0` (channel toggle).
- CP-OTFS anything = CP-OFDM at `X = ISFFT(x)` (unitary `isfft_dd`/`sfft_dd`), exact per realization.
- MCRB / pseudo-true / MMLE exist for CP-OFDM and CP-OTFS only; RCP waveforms are FIM/CRB only.
- Bound names: `CRB` genie; `MCRB0 = A⁻¹BA⁻¹` at the solved pseudo-true `η₀`; `BB = (η̄−η₀)(η̄−η₀)ᵀ`;
  `LB = MCRB0 + BB`. Affine in noise: build at `σ² = 1`, then `LB(N0) = N0·G + F + BB`.
- MP detector: per-edge damping **0.7**, 20 iterations (Raviteja WCNC'18 eq. 21); the old
  posterior-damping mode (0.3 / 30) exists only as a regression option.
- Common random numbers: `crn_stream(seed, sweep_idx, trial)`; never keyed on the waveform.

## Helpers
`.claude/skills/implement-paper` (one self-contained `impl/<paper>/` folder per reproduced paper) is current.
The stale `run-ber` command and `otfs-sim` skill (pointing to the retired `crb/EVA/`, swapped N/M,
wrong `ν_max ≈ 540 Hz`) were removed on 2026-11-09; use `crb/run/run_ber_eva` and `crb/run/run_ber_fig9`.
`.claude/settings.local.json` holds only Windows-path permissions from the old PC.

## Reference papers (`crb/papers/`)
- Raviteja, Phan, Hong, "Embedded Pilot-Aided Channel Estimation for OTFS in Delay–Doppler Channels," IEEE TVT 2019 (EVA, Fig. 9).
- Raviteja et al., "Low-complexity iterative detection for OTFS," IEEE WCNC 2018 (MP detector).
- Liu et al., "CP-OFDM Achieves Lower Ranging CRB Than Frequency-Spread Waveforms in the Large-Sample Regime," arXiv:2605.14451 (MT40, Fig. 1).

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
run_all_tests('full')                 % + Monte-Carlo closed-form checks (~3.5 min, 249 checks)
```
Headless: `"C:/Program Files/MATLAB/R2025b/bin/matlab.exe" -batch "cd('.../matlab/crb'); setup_crb; run_all_tests"`.
Drivers in `crb/run/` take `'tier', 'smoke'|'dev'|'paper'` (or a trial count) and write
`crb/results/*.mat` + `crb/figures/*.png`; every driver prints its acceptance table (PASS/WARN).
Parallel Computing Toolbox: drivers use `parfor` (8 workers). Optimization Toolbox: `lsqnonlin`
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

## Stale helpers (do not follow)
`matlab/.claude/commands/run-ber.md` and `matlab/.claude/skills/otfs-sim/SKILL.md` still point
to `crb/EVA/` (now `old/crb/EVA/`), swap N/M and quote `ν_max ≈ 540 Hz`. Use `crb/run/run_ber_eva`
and `crb/run/run_ber_fig9` instead.

## Reference papers (`crb/papers/`)
- Raviteja, Phan, Hong, "Embedded Pilot-Aided Channel Estimation for OTFS in Delay–Doppler Channels," IEEE TVT 2019 (EVA, Fig. 9).
- Raviteja et al., "Low-complexity iterative detection for OTFS," IEEE WCNC 2018 (MP detector).
- Liu et al., "CP-OFDM Achieves Lower Ranging CRB Than Frequency-Spread Waveforms in the Large-Sample Regime," arXiv:2605.14451 (MT40, Fig. 1).

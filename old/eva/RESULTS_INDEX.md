# Results index — `eva/results/*.mat`

Metadata read directly from each file's stored variables (via the `nTrials`,
`scenario`/`chan_desc`/`channel_set`, `configs`, `N`/`M`, `mod_type`,
`SNR_dbs`/`p_list` fields). Last refreshed 2026-08-17.
Note: MATLAB `.fig` files live in `eva/figures/`; `.png`/`.mat` in `eva/results/`.

## Current / relevant `.mat` files

| File | Script | Trials | Grid / mod | Channel(s) | Waveforms | Sweep | Notes |
|---|---|---|---|---|---|---|---|
| `mcrb_ce_ici_results.mat` | `mcrb_ce_ici.m` | 50 | N=512, M=128 (4-QAM) | reduced EVA (3 taps, ν_max=444 Hz) | OFDM-ICI / CP-OTFS-ICI | SNR 0–25 dB (26 pts) | CE=`dftdd`, **+BSC-matched arm**; realistic detector; MCRB + true BER/SER |
| `mcrb_ce_ici_results_2000.mat` | `mcrb_ce_ici.m` | **2000** | N=512, M=128, 4-QAM | reduced EVA (3 taps, ν_max=444 Hz, ICI) | OFDM-ICI / CP-OTFS-ICI | SNR 0–25 dB (26 pts) | high-trial confirm of the delay/Doppler placement asymmetry (Part 7) |
| `mcrb_ce_ici_nodopp_2000_results.mat` | `mcrb_ce_ici.m` (`doppler_scale=0`) | **2000** | N=512, M=128, 4-QAM | reduced EVA, **NO Doppler (ν=0, no ICI)** | OFDM-ICI / CP-OTFS-ICI | OFDM 0–25, **OTFS 0–15** dB | stores MCRB & BB separately (`rmse_tau_{mcrb,bb,lb,crb}`); Doppler-as-sensing-diversity (Part 10) |
| `mcrb_ce_ici_noici_1000t_results.mat` | `mcrb_ce_ici.m` (`use_ici=false`, `doppler_scale=1`) | 1000 | N=512, M=128, 4-QAM | reduced EVA, **Doppler ~444 Hz but NO ICI** | OFDM / CP-OTFS | OFDM 0–25, OTFS 0–15 dB | ICI negligible (=full-ICI within 0.3%, Part 11); delay/Doppler waveform split (Part 12); MCRB & BB separate |
| `compute_mcrb_mc_evace_results.mat` | `compute_mcrb_mc.m` (`scenario='evace'`) | 50 | N=512, M=128, 4-QAM | reduced EVA (3 taps, ν_max=444 Hz, random Doppler+phase, ICI) | OFDM-ICI / CP-OTFS-ICI | p 0–0.5 (26 pts) | **abstract-BSC twin of `mcrb_ce_ici`** on the SAME channel; SER-matched by `plot_bsc_vs_ce.m` |
| `compute_mcrb_mc_true_paper_results.mat` | `compute_mcrb_mc_true.m` | 10000 | N=128, M=1, 16-PSK | paper, 10 delay-only targets | OFDM / SC / OTFS | p 0–0.1 (51 pts) | true-param (BB=0); paper Fig. 1(b) |
| `compute_mcrb_mc_true_paperv2_results.mat` | `compute_mcrb_mc_true.m` | 5000 | N=128, M=16, 16-PSK | paperv2 (delay+Doppler, ν_max=39.1 kHz, ICI) | OFDM / SC / OTFS | p 0–0.1 (51 pts) | true-param; CRN-paired; delay **&** Doppler |
| `compute_mcrb_mc_paper_results.mat` | `compute_mcrb_mc.m` | 2000 | N=128, M=1, 16-PSK | paper, delay-only | OFDM / SC / OTFS | p 0–0.1 (51 pts) | pseudo-true (+BB) |
| `compute_mcrb_mc_paperv2_results.mat` | `compute_mcrb_mc.m` | 2000 | N=128, M=16, 16-PSK | paperv2 (delay+Doppler, ICI) | OFDM / SC / OTFS | p 0–0.1 (51 pts) | pseudo-true; delay & Doppler |
| `ber_ce_channels_results.mat` | `ber_ce_channels.m` | 20 | N=512×128 (eva) & 128×16 (pv2c), 4-QAM | eva, paperv2comms | OFDM {LS+interp / DFT / DFT+DD} / CP-OTFS-MP | SNR −10–35 dB (46 pts) | BER only (no MCRB) |
| `ber_ce_channels_nodopp_results.mat` | `ber_ce_channels.m` (`doppler_scale=0`) | 20 | same, **ν=0** | eva, paperv2comms | same as above | SNR −10–35 dB | no-Doppler BER: gap is diversity, no OFDM floor |
| `ber_ofdm_ce_vs_perfect_results.mat` | `ber_ofdm_ce_vs_perfect.m` | 500 | N=512, M=128, 4-QAM | reduced EVA (ν_max=444 Hz) | OFDM {perfect CSI / interp / DFT / DFT+DD} (+OTFS overlay) | SNR −10–50 dB (61 pts) | CE vs genie; DFT+DD reaches perfect CSI; ~4e-4 ICI floor (Part 8) |
| `ber_ofdm_ce_vs_perfect_nodopp_results.mat` | `ber_ofdm_ce_vs_perfect.m` (`nu_scale=0`) | 500 | same, **ν=0** | reduced EVA, no ICI | same | SNR −10–50 dB | no floor (→1e-6); confirms floor = ICI |
| `ber_awgn_results.mat` | `ber_awgn.m` | 500 | N=512, M=128, 4-QAM | **AWGN** (flat unit tap) | OFDM {perfect/interp/DFT/DFT+DD} / CP-OTFS-MP + theory | SNR −6–12 dB | validation: all collapse onto Q(√SNR) |
| `ber_rayleigh_flat_results.mat` | `ber_rayleigh_flat.m` | 1000 | N=512, M=128, 4-QAM | **single-path flat Rayleigh** (no Dopp/ICI) | OFDM (perfect CSI) / CP-OTFS-MP + theory | SNR 0–34 dB | Part 14a: OFDM = OTFS (no diversity) |
| `ber_{rayleigh,constmag}_eva_results.mat` | `ber_rayleigh_eva.m` (`fading`) | 1000 | N=512, M=128, 4-QAM | reduced-EVA 3-tap, Rayleigh **or** const-\|amp\| taps, no Dopp/ICI | OFDM (perfect CSI) / CP-OTFS-MP | SNR 0–34 dB | Part 14b–c: OFDM order-1; OTFS harvests delay diversity; `_dopp` variant adds Doppler+ICI |
| `mcrb_sweep_{grid,bw,guard,confirm}_results.mat` | `mcrb_sweep_grid.m` | 200–800 | N,M,Δf grid, 4-QAM | reduced EVA + Doppler+ICI | OFDM / CP-OTFS (LB ratio) | SNR 0–6 / 0–4 dB | Part 15: low-SNR LB ratio vs (N,M,Δf) / bandwidth / OTFS guard G / G13-vs-G9 confirm |
| `mcrb_confirm_g13_{ext,full}_results.mat` | `mcrb_sweep_grid.m` (G=13, Emax=8) | 500 | N=512, M=128, Δf=30 kHz (B=15.4 MHz) | reduced EVA + Doppler+ICI | OFDM / CP-OTFS | SNR 0–10 / **0–25** dB | Part 15c: full LB arc — OFDM Doppler edge 0–7 dB, OTFS wins elsewhere |
| `ber_confirm_config_g{13,9}_results.mat` | `ber_confirm_config.m` | 100–200 | N=512, M=128, Δf=30 kHz | reduced EVA + Doppler+ICI, OTFS guard G=13/9 | OFDM (dftdd) / CP-OTFS (MP) | SNR 0–25 dB | BER companion: G13 OTFS cliffs; G9 (shrunk guard) OTFS floors |
| `ber_vs_doppler_results.mat` | `ber_vs_doppler.m` | 20 | N=512, M=128, 4-QAM | reduced EVA, Doppler swept | OFDM (DFT CE) / CP-OTFS (MP) | v={0,3,10,30,60,120,240} km/h @ SNR={20,30} dB | BER vs mobility |

### Notes
- **Two channel families, not directly comparable.** The `paper`/`paperv2`
  MCRB files use the **abstract 16-PSK BSC** on the synthetic paper grid (`p` =
  symbol-error rate). The `mcrb_ce_ici` / `ber_*` files use the **realistic EVA
  channel** with a full CE + demod receiver (4-QAM, `SNR` sweep). Different
  grids, modulation, and error model (the Part 6 comparison in `EXPERIMENTS.md`).
- `true` files (10000/5000 trials) have more trials than the `pseudo-true` ones
  (2000): the pseudo-true `lsqnonlin` solve is slower.
- `ber_*` files store no `N`/`M` (grid is set per-channel inside the loop); those
  entries are filled from the scripts. `channel: reduced EVA [0 1 3] taps` for
  `mcrb_ce_ici` and `ber_vs_doppler` (they don't save a channel-name string).
- **BSC-vs-CE pairing.** `compute_mcrb_mc_evace_results.mat` (abstract i.i.d. BSC)
  and `mcrb_ce_ici_results.mat` (realistic detector) are the SAME reduced-EVA
  channel; `plot_bsc_vs_ce.m` matches the BSC `p` to the CE detector's SER at each
  SNR and overlays the two MCRBs (→ `mcrb_ce_ici_bsc_vs_ce_{delay,doppler,ser}.png`).
  Finding: at low SNR the real detector's errors bias the **delay ~3.7×** but the
  **Doppler only ~1.4×** above a same-count BSC — the errors sit in the delay-set
  frequency fades (delay reads frequency, so it feels the correlated placement;
  Doppler reads time, which the CE tracks out → BSC-like). Verified by
  `diag_error_geometry.m` (error clustering, freq > time) and the in-script
  `do_bsc` arm; above ~12 dB both ratios drop below 1 (deep-fade nulling). See the
  Part 6 discussion in `EXPERIMENTS.md`.

## Sidecar / backup files
| File | Note |
|---|---|
| `ber_ce_channels_results_otfsref.mat` | OTFS-curve backup, reused by `ber_ce_channels.m` when only re-running OFDM |
| `eva_fim_avg_results.mat` | `run_fim_avg.m` → `plot_crlb.m` (FIM/CRLB). **Currently the `paperdopp` (Doppler-only) run** → `plot_crb_doppler.m`/`crb_doppler_vs_snr.png`; regenerable for `paper`/`eva`/`paperv2` by switching the selectors (Part 13) |

## Stale — old `store_all=true` dumps (safe to delete, ~55 GB)
Full per-trial `nE²·nTrials·np·nC` arrays from **before** the paper scenarios
switched to `store_all=false`. Regenerable and unused by current work.

| File | Size | Date |
|---|---|---|
| `compute_mcrb_mc_true_results.mat` | 47 GB | Jul 6 |
| `compute_mcrb_mc_results_mc1000eva.mat` | 5.8 GB | Jul 2 |
| `compute_mcrb_mc_results_mc8000.mat` | 807 MB | Jun 30 |
| `compute_mcrb_mc_results_mc100eva.mat` | 584 MB | Jul 1 |
| `compute_mcrb_mc_results.mat` | 471 MB | Jul 5 |
| `compute_mcrb_mc_results_mc1000p.mat` | 447 MB | Jun 30 |

## Other / older experiments (small)
`mcrb_ce_ici_results_{0to20,100a,100b,ext,hi,lo}.mat`, `ber_ce_ici_results.mat`,
`compute_mmle_mc_results{,1000,1000v2,500}.mat`, `mcrb_avg_results.mat`,
`pseudotrue_diff{,_ofdm}_results.mat`, `fim_sweep_deltaf_*` (×3).

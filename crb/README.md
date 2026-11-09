# `matlab/crb/` — CRB, BER and MCRB/LB for OFDM vs OTFS

Self-contained MATLAB implementation of the derivations in `../../latex/` (`io/` → `fim/` →
`mcrb/`) on the two channels described in [`channels.md`](channels.md): **EVA** (9-path, delay +
Doppler) and **MT40** (40 static delay-only targets). Plan and status:
`~/.claude/plans/in-latex-folder-there-mossy-puppy.md`.

```
setup_crb                 % adds cfg/ io/ fim/ rx/ mcrb/ util/ run/ plot/ test/; asserts no name collisions
run_all_tests             % Phase 0 + Phase 1 identities (< 1 min)
run_all_tests('full')     % + Monte-Carlo closed-form checks (~1 min)
```
Headless (Linux, R2026b): `cd matlab/crb && /usr/local/bin/matlab -batch "setup_crb; run_all_tests"`.

## Conventions (asserted by `test/`)

| Item | Convention |
|---|---|
| Grid | `N` subcarriers = delay = rows, `M` symbols = Doppler = columns; TF `X` and DD `x` are `N×M`. (Raviteja's paper swaps the letters.) |
| Timing | `g = make_grid(N, M, deltaf, Tcp, 'cp'|'rcp')`: `T = 1/Δf`, sample `Ts = 1/(NΔf)`, symbol `TS = T + Tcp` (cp) or `T` (rcp). Never alias `Ts` and `TS`. |
| Parameters | `η_p = [α_R; α_I; τ; ν]`, stacked `4P` (`eta_pack/eta_unpack/eta_index`). |
| Kernel | `ici_kernel(q, nu, T)` returns the peak-1 kernel and the **reduced** derivatives `I' = (1/j2π)∂I/∂ν`, `I''`; limits `1, T/2, T²/3`. `isi_kernel` gives the RCP split. Both are `exp_moments` of `∫₀¹uᵏe^{Du}du`. |
| Transforms | `isfft_dd` / `sfft_dd` unitary; every CP-OTFS quantity = the CP-OFDM one at `X = ISFFT(x)`; `rcp-otfs-cont` = RCP-OFDM at `ISFFT(x)` (`rcp_ofdm_fim.tex:319`). |
| FIM | `fim(S, wf, g, cr, use_ici, sigma2) = (2/σ²)Re{GᴴG}`, `G = grad_stack(...)`. `crb_from_fim` inverts in bin units (`eta_scale`), masks the ν block when `ν ≡ 0`, and returns `kappa`. |
| Regimes | `use_ici = false` drops only the Doppler-induced leakage (derivative kernels set to zero); `ν = 0` is a different regime (the derivative kernels survive). |
| RNG | `crn_stream(seed, sweep_idx, trial)` (Threefry substreams keyed on trial, never on the waveform). |
| Tiers | drivers take `'tier'` = `'smoke'` (20) / `'dev'` (200–500) / `'paper'` (2000+). |

Waveform names: `cp-ofdm`, `cp-otfs`, `rcp-ofdm`, `rcp-otfs` (sampled branch, integer delays only),
`rcp-otfs-cont` (continuous integrate-and-dump branch). MCRB/pseudo-true exist for the CP
waveforms only.

## Layout

| Folder | Contents |
|---|---|
| `cfg/` | `make_grid`, `chan_eva`, `chan_mt40`, `chan_controls`, `chan_realize`, `eta_*`, `gen_symbols`, `liu_precoders`, `crn_stream` |
| `io/` | kernels (`exp_moments`, `ici_kernel`, `isi_kernel`), transforms, field builders (`fields_cp_ofdm`, `fields_rcp_ofdm`, `fields_rcp_otfs`), means (`io_*`), `rcp_otfs_time_sim` oracle |
| `fim/` | `grad_stack`, `fim`, `fim_expected`, `crb_from_fim`, `crb_pilot_only`, `ct_crlb_floor`, `ici_moments`, `power_sums`, `dirichlet_moment`, `kappa_filter` |
| `rx/` | `modem` (Gray 4/16-QAM, 16-PSK; no toolbox); OTFS `otfs_pilot_frame`, `otfs_threshold_ce`, `mp_detect_dd` (per-edge MP, Raviteja WCNC'18 damping); OFDM `ofdm_pilot_lattice`, `pilot_spacing_match`, `ofdm_ce` (interp / dft / dftdd), `ofdm_eq_1tap`, `ofdm_eq_banded`; BSC arms `bsc_symbols`, `bsc_match`; one-frame driver `ber_frame` |
| `mcrb/` | `mcrb_AB_cp_ofdm`, `mcrb_AB_cp_otfs` (A, B, J̃, b, A_Δ at any evaluation point), `pseudo_true` (VARPRO / full-4P / closed form), `stab_quad_bin`, `mcrb_assemble` + `mcrb_eval` (affine-in-noise pieces G, F, BB, Gc), `mcrb_AB_expected`, `covb_cp_ofdm`, `covb_cp_otfs` (16×8 oracles), `mmle_cp` (mismatched MLE, overlay only) |
| `run/` | CRB: `diag_cond_eva`, `run_crb_mt40_block` (Liu Fig. 1), `run_crb_eva`, `run_crb_mt40_frame`; BER: `run_ber_fig9` (Raviteja Fig. 9), `run_ber_eva`; MCRB: `run_covb_experiment`, `run_mcrb_bsc_eva`, `run_mcrb_rx_eva`, `run_mcrb_bsc_mt40` |
| `plot/` | `plot_liu_fig1`, `plot_crb_vs_snr`, `plot_fig9`, `plot_ber_vs_snr`, `plot_mcrb_vs_p`, `plot_mcrb_vs_snr`, `plot_ber_mcrb_panels`, `plot_mcrb_mt40` |
| `test/` | `run_all_tests`, `test_*`, `oracles/` (independent direct-DD builders), `oracle_paths` (adds `old/eva/utils` inside tests only) |
| `results/`, `figures/` | outputs (`.mat`, `.png`) |

## Old → new name map

| old/eva/utils | here |
|---|---|
| `get_ici_coeffs`, `_derivative`, `_derivative2` (full derivatives) | `ici_kernel` (reduced derivatives) |
| `get_isi_coeffs`, `_derivative` | `isi_kernel` |
| `isfft_u`, `sfft_u` | `isfft_dd`, `sfft_dd` |
| `ofdm_mean` (+3 private copies) | `io_cp_ofdm` via `fields_cp_ofdm` |
| `getFIM_OFDM_multi_{noICI,ICI}`, `getFIM_CP_OTFS_*`, `getFIM_RCP_OFDM_multi`, `getFIM_RCP_OTFS_multi` | `grad_stack` + `fim` |
| `getFIM_*_expected` | `fim_expected` |
| `crb_multi`, `pinv_equilibrated`, `stab_quad` (5 copies) | `crb_from_fim`, `stab_quad_bin` |
| `packEta` (6 closures), `unpackEta` | `eta_pack`, `eta_unpack`, `eta_index`, `eta_scale` |
| `match_pilot_spacing`, pilot code in `mcrb_ce_ici.m` | `pilot_spacing_match`, `ofdm_pilot_lattice`, `otfs_pilot_frame` |
| `build_effective_dd_eva` (WCNC eq. 13 with twist) | `dd_response_io` (taps from the io layer: no twist, `TS`) |
| `estimate_channel` / threshold CE in `mcrb_ce_ici.m` | `otfs_threshold_ce` |
| `mp_detector` (posterior damping 0.3, 30 it) | `mp_detect_dd` (`'edge'` = `impl/otfs_ber/mp_detect_otfs` transposed, 0.7, 20 it; `'posterior'` = old, bit-exact) |
| `ce_ofdm` (3 local copies) | `ofdm_ce` |
| `ofdm_banded_mmse` | `ofdm_eq_banded` |
| `bsc_qpsk`, `bsc_16psk`, `bsc_match` locals | `bsc_symbols`, `bsc_match` |
| `getAB_OFDM_core` (+8 wrappers), `getAB_CP_OTFS_*` | `mcrb_AB_cp_ofdm`, `mcrb_AB_cp_otfs` |
| `pseudoTrue_OFDM`, `pseudoTrue_CP_OTFS` | `pseudo_true` |
| `getAB_*_expected`, `getMCRB_covb_*` | `mcrb_AB_expected`, `covb_cp_ofdm`, `covb_cp_otfs` |
| `mmle_OFDM`, `mmle_CP_OTFS` | `mmle_cp` |
| LB assembly in `mcrb_ce_ici.m` (lines 240-315) | `mcrb_assemble`, `mcrb_eval` |

## Receiver knobs (not channels)

| Knob | Default | Where |
|---|---|---|
| OTFS pilot | `(l_τ, M/2)`, full-Doppler guard rows `0..2l_τ` (5248 cells at full EVA), `guard_extra` shrinks it | `otfs_pilot_frame` |
| pilot budget | `'matched'`: `E_pil = n_nondata · 10^(boost/10)` (data-referenced), shared by both waveforms; OFDM lattice matched in cell count within 5 % and in energy exactly | `ber_frame`, `ofdm_pilot_lattice` |
| OTFS CE | threshold `3·sqrt(N0)` on `|y|` in the pilot band | `otfs_threshold_ce` |
| DD spread | `E = 10`, `tap_tol = 0.02` (receiver approximation only) | `dd_response_io` |
| MP | per-edge damping 0.7, 20 iterations, single precision, early stop `1e-4` | `mp_detect_dd` |
| OFDM CE | `dftdd` (DFT-denoised + 2 decision-directed passes); `dft`, `interp` alternatives; genie modes `ideal1tap` (the exact diagonal `Σ_p α_p I_p(0) e^{…}` of the ICI operator, ICI treated as noise) and `idealbanded` (`W_ici = 8`) | `ofdm_ce`, `ofdm_eq_*`, `ber_frame` |
| BSC | Gray-bit flips for 4-QAM / 16-PSK, uniform symbol errors for 16-QAM; `bsc_match` re-scatters the real detector's error count i.i.d. (`'uniform'`: random different point, the legacy control; `'nearest'`: one Gray bit flipped, placement-only control) | `bsc_symbols`, `bsc_match` |
| pseudo-true | VARPRO over `[τ; ν]` in bin units, Kaufman Jacobian, `lsqnonlin` (`ftol 1e-9`, `steptol 1e-8`, `opttol 1e-7`), `ν` masked when `ν ≡ 0`; `info.exitflag` is stored by every driver and a non-converged solve invalidates the trial | `pseudo_true` |
| trial statistics | every MCRB driver keeps the per-trial traces and reports **both** the mean and the median of the per-trial composite bounds (`mcrb_stats`: `N0·G_t + F_t (+ BB_t)` first, then the statistic, never medians of the pieces). Headline = median for Rayleigh gains (the trial mean of `J⁻¹` has no finite expectation: `tr G ∝ 1/|α|²`), mean for constant-magnitude gains (the legacy references). `valid` = `pc.valid` ∧ converged pseudo-true ∧ `kappa ≤ kappa_thresh`; the drop fraction is stored and printed | `mcrb_stats`, `run_mcrb_*` |
| common random numbers | `crn_stream(seed, key, trial)`; `run_ber_eva` and `run_mcrb_rx_eva` share seed 9100 and key `1000 + 10·SNR_dB`, so they see identical channel / noise / data draws at equal SNR | `crn_stream` |
| options | every driver rejects unknown option names (`parse_opts`); `'tier'` may be a name or a trial count; `'tag'` is appended to the output file names | `parse_opts` |

## Cov(b) experiment (plan F3.0, `run_covb_experiment`, 40 000 trials, 16×8 legacy preset, p = 0.1, no ICI)

- **Single target:** the Kronecker-trace closed form of `cp_otfs_mcrb.tex` §4 evaluates to the
  collapse `4p(1−p)|α|²J/σ²` to machine precision (1e-16), as does the OFDM form, and the MC
  `Cov(b)` agrees within its standard error (3e-3 vs SE 6e-3). So for `P = 1` there is **no
  contradiction** between `cp_otfs_mcrb.tex` §3-4 and `mcrb_expectation_comparison.tex` §7.3.
- **Three targets:** both closed forms match MC (7e-3 vs SE 1e-2); the OFDM and OTFS closed forms
  differ by 3.2 % on the kinematic block. The paired (common-random-number) MC difference
  OTFS − OFDM is non-zero: kinematic-block trace `t = 3.7` (99 % critical value 3.25 with 9 batch
  dof), per-entry `max |z| = 9.1` over the 21 unique entries (`χ² = 383`), and its trace
  `6.7e12 ± 1.8e12` agrees with the closed-form difference `9.4e12`. Conclusion: the collapse
  `Cov(b) = 4p(1−p)|α|²J/σ²` is exact only for `P = 1`; for `P > 1` the OTFS `Cov(b)` of
  `cp_otfs_mcrb.tex` §3-4 is the right one and differs from the OFDM one by a few percent
  (the §7.3 extension to CP-OTFS holds only in the single-target no-ICI case). Runtime code never
  uses either closed form (all drivers use the realization `B`), so nothing downstream depends on it.

Deviations from the old code that are deliberate: the sampled RCP-OTFS builder has **no** OFDM-side
ICI kernel (the old `use_ici=true` path double-counted intra-symbol Doppler; the old `use_ici=false`
path is the exact sampled model and is what `rcp-otfs` with `use_ici=true` reproduces); the DD
receiver taps will come from the io layer (no twist phase, `TS` in the Doppler kernel).

## Results and comparisons

- [`eva_reduced_results.md`](eva_reduced_results.md): BER, genie CRB and MCRB/LB on the
  **EVA-reduced** channel (the three largest-gain EVA paths), all on the same 0–24 dB grid and
  shared channel draws.
- [`grid_sweep_results.md`](grid_sweep_results.md): whether the grid `(N, M)` can widen OFDM's
  delay margin over OTFS. It cannot: every bias falls as `(NM)^-1/2`, so each margin is a ratio of
  prefactors the grid does not set. The anisotropy itself does follow `N^0.34 M^0.23`, but in the
  direction that *closes* the delay margin. Includes the theory-versus-measurement scorecard.
- [`eva_reduced_fast_results.md`](eva_reduced_fast_results.md): **EVA-reduced-fast** (the same
  channel at 240 km/h), the variant where OFDM's lower bound beats OTFS's in *both* the delay and
  the Doppler block while the genie CRBs stay equal and OTFS keeps the lower BER; includes the
  delay/Doppler anisotropy measurement and the channel-parameter sweep behind it.

- [`mcrb_conclusions.md`](mcrb_conclusions.md): what the MCRB/LB runs say, cross-read against the
  CRB and BER results (F = 0 at the pseudo-true point so LB = N0·G + BB; the 1/(1−2p)² shrink law;
  BB dominance and the SER ≈ 1e-2 threshold; OFDM vs OTFS; error placement; data-aided vs
  pilot-only), from the 100-trial real-detector run and the BSC/MT40 sweeps.

## Comparisons with the papers

- [`liu_comparison.md`](liu_comparison.md): our MT40 CRB results vs Liu et al. Fig. 1 (exact
  vector-extracted paper data; OFDM gains agree to 0.01–0.08 dB, identities exact, a common 12 %
  offset of the absolute levels from the unstated placement rule), plus the EVA-channel answer
  (no OFDM/OTFS CRB gap on frame-sized grids; the gap scales as unknowns/samples).
- [`raviteja_comparison.md`](raviteja_comparison.md): our EVA BER results vs Raviteja et al.
  Fig. 9 (exact vector-extracted paper data; all curves within the frame-level confidence
  intervals at 10–15 dB with a 200-trial check at the paper's 15 dB vertex, OFDM at all SNRs;
  the 16–18 dB OTFS tail is governed by ~2 %-rate tap-truncation floor frames and needs ~1000
  frames per point).

## Status (2026-09-07)

- All three phases are implemented; `run_all_tests('full')` = 249/249 checks (identities, finite
  differences of every gradient and Hessian column, 1e-12 regressions against `old/eva/utils`,
  independent direct-DD oracles, Monte-Carlo closed forms, Liu identities, MP bit-exactness);
  `run_all_tests('slow')` adds the legacy end-to-end MCRB regression (~20 min).
- The Phase-3 drivers were reviewed (four lenses, adversarially verified) and corrected on
  2026-09-07: per-trial composite statistics (`mcrb_stats`), median headline for Rayleigh gains,
  validity/drop fraction with the pseudo-true exit flag, entrywise mean-field check, MT40-1D
  framing with Liu's precoders and the F1.1/F1.3 cross-checks, pilot-only identifiability guard,
  the exact one-tap genie diagonal, the paired flat-Rayleigh check, the calibrated Cov(b) test,
  unknown-option rejection, and shared draws between the BER and MCRB drivers.
- Conditioning gate (`results/diag_cond_eva.mat`): full 9-path EVA with Rayleigh gains passes
  (0 % of draws with `kappa > 1e10`, 0.5 % above 1e8), so full EVA is the bound channel;
  EVA-reduced is the fast default of the MCRB drivers (`'reduced', false` switches).
- Figures produced so far (smoke tier unless noted): `fig_liu_fig1.png` (F1.1, Liu Fig. 1 levels
  and identities reproduced), `fig_crb_eva_smoke.png` (F1.2), `fig_crb_mt40_*.png` (F1.3),
  `fig_ber_fig9_*.png` (F2.1), `fig_ber_<channel>_smoke.png` (F2.2), `fig_covb_experiment.png`
  (F3.0), `fig_mcrb_vs_p_*.png` (F3.1), `fig_mcrb_rx_*.png` + `fig_ber_mcrb_panels_*.png` (F3.2),
  `fig_mcrb_mt40_*.png` (F3.3). Each driver prints its acceptance table; see the run logs.
- **Fig. 9 OTFS-ideal tail (F2.1) — settled.** Our 50-trial run (`results/ber_fig9_dev.mat`)
  gives 1.4e-2 / 4.2e-3 / 3.7e-4 / 2.7e-5 / 2.1e-6 at 10–18 dB against the `impl/otfs_ber` README
  table 1.8e-2 / 3.2e-3 / 9.7e-4 / 2.9e-4 / 1.2e-4 (OFDM-ideal agrees within 7 % everywhere, CE
  curves and ordering agree); an independent-seed re-run of the impl code (seed 2020, 48 trials)
  gives 1.1e-2 / 2.8e-3 / 5e-4 / 2e-4 / 1e-4. A paired per-frame comparison on those 48 draws
  (`scratchpad/fig9_paired.m`, identical channel / data / noise per trial) shows: (i) the MP
  detectors are bit-identical; (ii) the io-layer chain is slightly *worse* than impl's DD-twist
  chain on every frame (mean 4.1e-4 vs 1.4e-4 at 18 dB), so the channel model is not what lowers
  our tail; (iii) at 16–18 dB the entire impl tail is ONE frame (trial 16: BER 6.6e-3 at 14, 16
  and 18 dB, i.e. an SNR-independent floor; 2.0e-2 in the io chain) — a strong channel (total
  power 2.3) whose second-strongest path sits at a half-integer fractional Doppler (`k_real`
  −3.50, 3.49), where the Dirichlet inter-Doppler spread is widest and the `E = 10` tap
  truncation (Raviteja's, shared by both codes) leaves a signal-proportional model mismatch.
  Without that frame both chains are at ~1e-5. Such frames occur in ~2 % of EVA draws, so
  50-trial tails differ by 10× by chance (our dev run happened to draw none at 16–18 dB; the
  seed-2019 impl run also carried a deep-fade draw reused for all `SNR_d`). The expected
  high-SNR floor is therefore ≈ 0.02 × 7e-3 ≈ 1e-4, and a paper-tier run (≥ 200 trials) or a
  larger `E` is needed for a tail below that. Compare only points with ≥ 20× the measurement
  floor and expect the 16–18 dB points to be dominated by rare truncation-floor frames.
- **Legacy end-to-end regression (`test_regression_old`, `run_mcrb_rx_eva('legacy', true)`, 50
  trials, EVA-reduced, constmag, `T_cp = 0`, dftdd, posterior MP 0.3 × 30) against
  `old/eva/results/mcrb_ce_ici_results_2000.mat`:** OFDM BER within 1 % at every SNR
  (2.29e-1 … 1.32e-3), genie delay CRB within 3.3 %, delay LB within 25 % except two points (OTFS
  0 dB 27 %, OFDM 12 dB 25 %; the trial mean of LB is heavy-tailed through BB and the old 50-trial
  file itself differs from the 2000-trial one by up to 30 % there), OTFS LB within 3 % above
  12 dB. OTFS BER in this configuration is 1.3–1.5× the old value at 10–12 dB (5.1e-3 vs 3.7e-3,
  4.9e-4 vs 2.8e-4): the old chain applied the DD-circular twist model, this one the io layer.
- **MCRB results at smoke tier (EVA-reduced, 20 trials, median over trials).** F3.1: `p = 0`
  identities exact; MCRB0 monotone in `p`; `BB` dominates MCRB0 at 20 dB for every `p ≥ 0.05`
  (×10–×100); OFDM and OTFS MCRB0 within 1.5 % for `p ≤ 0.3`; the entrywise mean-field check
  `E{A(η̄)} = −(1−2p)E{J}` holds within 2 % on every within-target entry for `0 < p ≤ 0.2`
  (no-ICI, and the four asserted ICI entries; the reported gain-Doppler and (ν,ν) ICI deviations
  are 0.2–0.5 %). F3.2 (real detector, dftdd / MP): the two waveforms' genie CRBs agree to 0.2 %,
  MCRB0 is within 0.2 % of the CRB above 18 dB, LB/CRB = 1.01–1.05 at 25 dB, and the OFDM
  real/BSC-matched LB ratio goes from 1.5 at 0 dB to 0.73 at 20–25 dB (the old "real ≫ BSC at low
  SNR, real < BSC at high SNR" pattern, weaker than the old 2.5–3.7× because Rayleigh gains and the
  median statistic are used here; the legacy configuration is reproduced by `'legacy', true`).
  OTFS's Doppler LB is 4–12× above OFDM's below 10 dB (MP failures bias ν, as in the old log) and
  equal above 15 dB.
- **MT40 MCRB (F3.3, smoke tier, 20 data draws, median).** `p = 0` at `L = 40` reproduces the
  stored F1.1 genie CRBs (1-D, Liu precoders) within 0.5–4.5 % and the F1.3 ones (2-D frame)
  within 1 %, the Jensen bound exactly; the ordering OFDM < OTFS < SC (1-D) and OFDM < OTFS (2-D)
  holds at `p = 0` for both constellations (OFDM/OTFS 0.67 / 0.80 in 1-D). At `L = 10`, `p > 0`
  (20 draws) the MCRB0 ordering is kept for 16-PSK in 2-D but not strictly in 1-D or for 16-QAM
  (ratios within 0.8–1.2, i.e. inside the draw-to-draw scatter). The data-aided LB crosses the
  20 %-comb pilot-only CRB at 20 dB at `p ≈ 0.5–1 %` (16-PSK) and `1–2.5 %` (16-QAM), the old
  "≈ 1 % at 20 dB"; at 0 dB it stays below it for all `p ≤ 0.1`. `LB / MCRB@η̄` = 1.1–3.7
  (old 1.2–3.1). With `L = 40` the 20 % comb (26 pilots, 52 real equations) cannot identify the
  120 unknowns, so the pilot-only line is omitted there.
- **RCP-OTFS sampled branch on MT40-int (F1.3).** With 40 integer-spaced targets the sampled
  RCP-OTFS delay-CRB trace (0.018 samples²) is ~3.5× below the per-symbol-CP waveforms and the
  continuous RCP-OTFS branch (0.068), although the single-target FIMs of all branches agree to 1.6 %.
  The integer geometry itself is harder for the CP waveforms (0.031 with continuous delays vs 0.068
  with rounded ones, `kappa` 3e4 vs 6e3): the sampled model's Doppler-index-dependent wrap phase
  decorrelates the 40 delay-derivative columns. Reported, not asserted (plan R13).

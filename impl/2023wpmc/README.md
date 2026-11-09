# CRLB Analysis for OTFS and OFDM — WPMC 2023

Self-contained MATLAB reproduction of

> B. Wang, J. Zhu, X. She, P. Chen,
> **"Cramér-Rao Lower Bound Analysis for OTFS and OFDM Modulation Systems,"**
> 2023 26th Int. Symp. on Wireless Personal Multimedia Communications (WPMC).
> (`2023wpmc.pdf` in this folder.)

The paper derives the Cramér-Rao lower bound (CRLB) of **delay (τ)** and
**Doppler (ν)** estimation for four point-to-point signal models and shows that
OTFS outperforms OFDM for sensing.

## What is implemented

Four closed-form channel-coefficient matrices `Ψ^p` (the unified model
`y = Σ_p Ψ^p x + n`, eq. (27)), each transcribed from the paper:

| scheme | file | paper eq. | global const |
|---|---|---|---|
| CP-Free OFDM   | `psi_cpfree_ofdm.m` | (6),(7)   | `1/√M` |
| CP-OFDM        | `psi_cp_ofdm.m`     | (10),(11) | `1/√M` |
| Zak-OTFS       | `psi_zak.m`         | (38),(40) | `1/√(MN)` |
| Two-Steps OTFS | `psi_twosteps.m`    | (25),(26) | `1/√(MN)` *(see note)* |

Each builder reproduces the paper's Dirichlet kernels, phase terms and the
I₁/I₂ inter-symbol-interference split. The derivations were re-derived from
Sec. II and match the paper's printed `Ψ` structure.

**FIM / CRLB** (`compute_crlb_scheme.m`): per path the unknowns are
`θ = [|h|, ∠h, τ, ν]` (Sec. III). With AWGN `n ~ CN(0, σ²I)` in the estimation
domain (eq. (46)) and `σ² = 10^(−SNR/10)`, the Slepian–Bangs FIM is
`J_ij = (2/σ²)·Re{(∂μ/∂θ_i)^H (∂μ/∂θ_j)}`, `μ = Ψx` (eq. (48)).
`∂/∂|h| = Ψ/|h|` and `∂/∂∠h = jΨ` analytically; `∂/∂τ`, `∂/∂ν` by central
finite differences of the exact closed-form `Ψ` (equivalent to eqs.
(30)–(36),(43)–(44),(51)–(52), but far less transcription-prone — and safe
because τ stays within one I₁/I₂ piece for the step used).

Because the four parameters differ by ~20 orders of magnitude (τ~1e-6, ν~500,
|h|~1), the FIM is severely scaled. We apply **Jacobi (diagonal)
preconditioning** before inversion (exact, not an approximation):
`[J⁻¹]_ii = d_i² [J̃⁻¹]_ii` with `d = 1/√diag(J)`. This drops `cond` from
~1e17 to ~30.

The data `x` is a **fixed random unit-power 4-QAM block** (`NMC = 1`, seeded).
Set `NMC > 1` in `run_sim.m` to Monte-Carlo average the FIM (→ the `E{xxᴴ}=I`
form). CRLB curves are smooth in SNR regardless, since for fixed data the FIM
kernel is SNR-independent and only `σ²` scales.

## Figures

Run, then inspect `fig_1.png`, `fig_2.png`, `fig_3.png`:

```
matlab -batch "cd('impl/2023wpmc'); run_sim"
```

- **fig_1** — CRLB of τ and ν vs SNR, single path, `N=M=12, τ=3.33 µs,
  ν=500 Hz, SCS=15 kHz` (paper Fig. 1).
- **fig_2** — CRLB vs SNR for sizes `M=N ∈ {8,12,16}` (paper Fig. 2: effect of
  PRB size; larger MN ⇒ lower CRLB).
- **fig_3** — CRLB vs SNR, two-path channel `h₁=0.7e^{jπ/3}, h₂=0.3e^{j3π/4},
  τ₁=3.33 µs/ν₁=500 Hz, τ₂=5 µs/ν₂=2.5 kHz` (paper Fig. 3).

## Validation against the paper

Reproduced qualitative claims (CRLB @ 0 dB, Fig. 1 / Fig. 3):

- **OTFS ≫ OFDM for sensing** — both OTFS schemes are ~12× lower in delay and
  ~10× lower in Doppler than OFDM. ✓ (paper's central result)
- **Zak-OTFS best for delay** — τ-CRLB 2.232e-16 < Two-Steps 2.279e-16. ✓
- **Two-Steps OTFS best for Doppler** — in the two-path case (Fig. 3),
  ν-CRLB 28.6 < Zak 29.5. ✓ (paper: "Zak / Two-Steps best for delay / Doppler
  respectively")
- **CP-Free ≈ CP-OFDM** — curves within ~1 dB. ✓ (paper: "comparable")
- **Stronger path ⇒ lower CRLB** — Fig. 3 path 1 (|h|=0.7) below path 2
  (|h|=0.3). ✓
- **Slope** — CRLB ∝ 1/SNR (−1 decade per 10 dB), straight lines on the
  semilog axes. ✓

The CRLB **magnitudes** depend on the paper's unstated SNR convention
(`σ² = 10^(−SNR/10)`, transmit-referenced), so absolute levels may be offset
from the paper by a constant; trends, crossovers and ordering match.

## Modeling notes / discrepancies

1. **Two-Steps constant.** Eq. (26) prints `1/MN`; we use `1/√(MN)`. Zak-OTFS
   and Two-Steps share the *same transmit signal* and differ only by an
   *invertible* receiver, so their CRLBs must be ~equal. `1/(MN)` makes
   Two-Steps exactly `MN = 144×` worse than Zak — physically impossible and
   opposite to the paper's own conclusion. `1/√(MN)` (matching eq. (40))
   restores Zak ≈ Two-Steps, with the small residual difference coming from the
   I₁/I₂ split and the `e^{j2πντ}` coupling. We interpret the printed `1/MN`
   as a typo for `1/√(MN)`.

2. **Receiver invertibility.** Under a physically consistent (time-domain)
   noise model, an invertible receiver cannot change the CRLB; the scheme
   differences the paper reports follow from adding `σ²I` noise in each
   scheme's own (non-unitary) estimation domain, i.e. they are driven by the
   per-scheme `Ψ` normalization plus the TX waveform (OFDM identity vs OTFS
   ISFFT). This implementation follows the paper's output-domain `σ²I` model.

3. **`f_c` does not enter.** The carrier-delay phase `e^{−j2πf_cτ}` is a
   constant per-path phase absorbed into `∠h_p` (an estimated parameter); it
   does not appear in any `Ψ` (consistent with eqs. (7),(11),(26),(40)). The
   `f_c = 3 GHz` in the captions is therefore not used in the CRLB.

## Files

```
2023wpmc.pdf            source paper
run_sim.m               driver: produces all three figures + results.mat
sys_params.m            system parameter struct (M,N,Δf,T,CP,…)
psi_build.m             dispatcher -> per-scheme Ψ^p
psi_cpfree_ofdm.m       Ψ^p, CP-Free OFDM   (eq. 7)
psi_cp_ofdm.m           Ψ^p, CP-OFDM        (eq. 11)
psi_zak.m               Ψ^p, Zak-OTFS       (eq. 40)
psi_twosteps.m          Ψ^p, Two-Steps OTFS (eq. 26)
compute_crlb_scheme.m   FIM + CRLB for one scheme (all SNR points)
fig_1.png fig_2.png fig_3.png   reproduced figures
results.mat             saved CRLB data
```

## Status

Complete. All three figures reproduce the paper's trends and scheme ordering.
The only deviation from the printed equations is the Two-Steps constant
(`1/MN` → `1/√(MN)`), justified above.

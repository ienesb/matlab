---
name: implement-paper
description: Use when implementing a research paper in MATLAB to reproduce its figures. Each paper lives in its own impl/<subdir>/ folder containing the source PDF. Trigger on requests like "implement the paper in impl/<name>", "reproduce Figure N from this paper", or working with any PDF under impl/. Enforces a self-contained, single-directory implementation.
allowed-tools: Bash(matlab:*), Read, Write, Edit, Glob, Grep
---

# Implement a Paper in MATLAB

Reproduce a research paper's results/figures in MATLAB, one paper per
`impl/<subdir>/` folder.

## Hard rules (non-negotiable)
1. **Self-contained.** ALL code for a paper lives inside its own
   `impl/<subdir>/` folder. Create every `.m` file there.
2. **No external scripts.** Do NOT call, `addpath`, copy, or depend on any
   script outside that subdirectory — not from `crb/`, `otfs/`, `ltfat/`,
   or another paper's folder. If a utility is needed, reimplement it locally.
   The only allowed externals are core MATLAB / standard toolbox functions.
3. **Match the figures.** The goal is reproducing the paper's figures as
   closely as the paper allows. Treat each target figure as the acceptance test.

## Workflow
1. **Locate & read the PDF.** Find the `*.pdf` in the target subdirectory and
   read it (use the `pages` arg for long PDFs). Extract:
   - System model, key equations, and the algorithm(s) to implement
   - All numerical parameters (sizes, SNR ranges, iteration counts, seeds)
   - Exactly which figures are being reproduced and what each axis shows
2. **Plan.** Briefly summarize the algorithm and list the figures to produce.
   Confirm the parameter table with the user if anything is ambiguous or
   missing from the paper.
3. **Implement (inside the subdir only):**
   - `main_<paper>.m` or `run_sim.m` — top-level driver
   - One function per logical block (model, estimator/detector, metrics)
   - Hard-code paper parameters with comments citing the paper section/eq #
4. **Reproduce figures.** Generate each target figure with `plot`/`semilogy`
   matching the paper's axes, legends, and scale. Save to the subdir
   (e.g. `fig_<N>.png` + the data `.mat`).
5. **Validate.** Compare against the paper: do trends, crossover points, and
   magnitudes match? Note any discrepancy and the likely cause
   (missing parameter, integer rounding, normalization, etc.).

## Per-paper folder layout (target)
```
impl/<subdir>/
  <subdir>.pdf          # source paper (provided by user)
  run_sim.m             # driver: runs everything, produces all figures
  <local functions>.m   # self-contained helpers
  results.mat           # saved simulation data
  fig_*.png             # reproduced figures
  README.md             # what's implemented, which figures, how to run, status
```

## MATLAB conventions
- 1-based indexing, vectorized ops, `parfor` for Monte Carlo loops.
- Set and record RNG seed for reproducibility.
- Run headless to check: `matlab -batch "cd('impl/<subdir>'); run_sim"`.
- Keep style consistent within the paper's folder.

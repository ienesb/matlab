---
description: Run the EVA BER simulation and summarize the resulting curves
argument-hint: [snr-range, e.g. 10:2:18]
allowed-tools: Bash(matlab:*), Read, Glob
---

Run the EVA BER message-passing simulation.

Target SNR range (optional): $ARGUMENTS

## Steps
1. Confirm the channel/sim parameters in @crb/EVA/initialize_parameters.m
   (expect N=512, M=128, 4-QAM, EVA 9-path). Note anything that differs.
2. Execute the main script from `crb/EVA/`:
   `matlab -batch "cd('crb/EVA'); eva_ber_mp2"`
   If an SNR range was provided above, mention that it should override the
   default `SNRdB` vector before running.
3. After it finishes, read `crb/EVA/eva_ber_mp_results.mat` (or list the
   newest `*_results.mat`) and report the BER values for each curve:
   OFDM(CSI), OTFS(CSI+MP), OTFS(CE+MP).
4. Flag any curve that looks like an error floor or diverges from the
   expected Raviteja 2019 Fig. 9 trend.

Keep the summary to a short table plus one line of interpretation.

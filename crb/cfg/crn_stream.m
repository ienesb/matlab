function rs = crn_stream(seed, sweep_idx, trial)
% CRN_STREAM  Counter-based random stream for common-random-number pairing.
%   rs = crn_stream(seed, sweep_idx, trial)
% Substream = 10000*sweep_idx + trial is keyed on the sweep point and the
% trial only -- NEVER on the waveform/config -- so every waveform sees the
% same data, gains, Doppler and noise draws (old/eva convention).  Threefry
% streams are reproducible inside parfor.  sweep_idx >= 0, trial >= 1.
    if nargin < 2, sweep_idx = 0; end
    if nargin < 3, trial = 1; end
    rs = RandStream('Threefry', 'Seed', seed);
    rs.Substream = 10000*sweep_idx + trial;
end

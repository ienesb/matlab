function run_phase2_campaign(varargin)
% RUN_PHASE2_CAMPAIGN  The simulation queue of submission-plan Phase 2
% (latex/submission_plan.md), most useful points first.  Every entry is one
% resumable run_phase2_eva call restricted to a few SNR points ('pass'); a
% finished entry returns immediately, so the same command restarts the queue
% after a stop, a crash or a reboot:
%
%   cd matlab/crb
%   setsid nohup nice -n 10 /usr/local/bin/matlab -batch "setup_crb; run_phase2_campaign" \
%       >> results/phase2_campaign.log 2>&1 &
%   touch results/STOP_PHASE2      % stop cleanly after the current chunk (delete it to resume)
%
% Runs (QPSK, 512 x 128, standard receivers, controls 'uniform' and 'nearest';
% every run is laid out for 500 frames per point and the campaign fills the
% first 'frames' = 250 of them, see below):
%   main   120 and 240 km/h, EVA-reduced, SNR 0:2:24 dB, seed 9100, estimator
%          arms at 0:4:20 dB.  Frames 1-100 (120 km/h) and 1-60 (240 km/h)
%          are the frames of the September paper runs.
%   sweep  60 and 360 km/h, SNR 0:4:16 dB, seed 9100, estimator arms at 4, 8 dB.
%   conf   all four speeds, SNR 4 and 8 dB, fresh seed 9101.
%   eva9   full 9-path EVA, 120 km/h, SNR 0:4:12 dB, seed 9100, estimator arms at 4, 8 dB
%          started at the true parameter only (the periodogram/CLEAN start does not
%          resolve the nine paths: 8-frame pilot of 2026-11-18).
%   genie  120 and 240 km/h, SNR 4 and 8 dB, seed 9100, detectors with the true
%          channel ('csi','genie'; the receiver-tier benchmark of latex/todo.md).
% Stages, in this order: A = 4 and 8 dB of main and sweep (the speed trend, the
% mechanism control and the estimator at the comparison points); B = 0-12 dB of
% main; C = conf at 120 and 240 km/h; D = rest of sweep; E = eva9 at 4 and 8 dB;
% F = 16, 20, 24 dB of main.  Not in the default (trimmed on 2026-11-18 to end
% the campaign in about 30 h): G = genie; X = conf at 60 and 360 km/h, eva9 at
% 0 and 12 dB, main at 14, 18, 22 dB.  L (added 2026-11-20) = the 60 km/h sweep
% again with the OFDM pilot lattice pinned to (37, 2), the lattice of the other
% three speeds: the search in pilot_spacing_match returns (4, 19) at 60 km/h,
% which makes that run incomparable (file suffix _lat37x2; the first 60 km/h
% file is kept).
%
%   run_phase2_campaign                                  % stages 'ABCDEF', 250 frames per point
%   run_phase2_campaign('stages', 'GX')                  % the trimmed parts, 250 frames
%   run_phase2_campaign('stages', 'L')                   % 60 km/h with the pinned OFDM lattice
%   run_phase2_campaign('stages', 'A', 'frames', 500)    % top up the 4 and 8 dB points to 500 frames
% 'frames' is the number of frames per point to reach in this call; the files
% stay CKPT_p2_*_n500_* until every point of a run has 500, and frames 1-250
% are the same frames whenever they are computed.  'tier' other than 500 is
% for rehearsals only (it changes the file names).
    o = struct('stages', 'ABCDEF', 'workers', [], 'tier', 500, 'eva9_tier', [], 'frames', 250);
    o = parse_opts(o, varargin, 'run_phase2_campaign');
    if isempty(o.eva9_tier), o.eva9_tier = o.tier; end
    grid = 0:2:24;
    main  = @(v) {'tier', o.tier, 'v_kmh', v, 'snr_db', grid, 'seed', 9100, 'mmle_snr', 0:4:20};
    sweep = @(v) {'tier', o.tier, 'v_kmh', v, 'snr_db', 0:4:16, 'seed', 9100, 'mmle_snr', [4 8]};
    conf  = @(v) {'tier', o.tier, 'v_kmh', v, 'snr_db', [4 8], 'seed', 9101};
    sweep60p = [sweep(60), {'ofdm_lattice', [37 2]}];
    genie = @(v) {'tier', o.tier, 'v_kmh', v, 'snr_db', [4 8], 'seed', 9100, 'csi', 'genie'};
    eva9  = {'tier', o.eva9_tier, 'v_kmh', 120, 'reduced', false, 'snr_db', 0:4:12, 'seed', 9100, 'mmle_snr', [4 8], ...
             'mmle_arms', {'e2e_truth', 'fresh_truth'}};
    Q = { ...
        'A', main(120),  [8 4];   'A', main(240),  [8 4];   'A', sweep(60), [8 4];   'A', sweep(360), [8 4]; ...
        'B', main(120),  [0 12];  'B', main(240),  [0 12];  'B', main(120), [2 6 10]; 'B', main(240), [2 6 10]; ...
        'G', genie(120), [8 4];   'G', genie(240), [8 4]; ...
        'C', conf(120),  [8 4];   'C', conf(240),  [8 4]; ...
        'D', sweep(60),  [0 12 16]; 'D', sweep(360), [0 12 16]; ...
        'E', eva9,       [8 4]; ...
        'F', main(120),  [16 20 24];  'F', main(240), [16 20 24]; ...
        'X', conf(60),   [8 4];   'X', conf(360),  [8 4];   'X', eva9, [0 12]; ...
        'X', main(120),  [14 18 22];  'X', main(240), [14 18 22]; ...
        'L', sweep60p,   [8 4];   'L', sweep60p,   [0 12 16]};
    if isempty(gcp('nocreate'))
        c = parcluster('Processes');
        if ~isempty(o.workers), c.NumWorkers = o.workers; end
        parpool(c, c.NumWorkers);
    end
    t0 = tic;
    for k = 1:size(Q, 1)
        if ~contains(o.stages, Q{k, 1}), continue; end
        fprintf('\n=== %s  campaign stage %s, entry %d of %d, pass %s, %d frames (elapsed %.2f h)\n', string(datetime('now')), Q{k, 1}, k, size(Q, 1), mat2str(Q{k, 3}), o.frames, toc(t0)/3600);
        run_phase2_eva(Q{k, 2}{:}, 'pass', Q{k, 3}, 'frames', o.frames);
    end
    fprintf('\nPHASE2_CAMPAIGN_FINISHED stages %s, %d frames, at %s (%.2f h in this call)\n', o.stages, o.frames, string(datetime('now')), toc(t0)/3600);
end

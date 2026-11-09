function plot_ber_ce_ici(matfile)
% PLOT_BER_CE_ICI  Load the saved BER results from ber_ce_ici.m and (re)draw the
% BER-vs-SNR_d figure.  Kept separate so the (long) simulation and the plotting
% are decoupled -- restyle here without re-running the sweep.
%
% Usage:
%   plot_ber_ce_ici                       % loads eva/ber_ce_ici_results.mat
%   plot_ber_ce_ici('other_results.mat')  % load a specific results file

    thisdir = fileparts(mfilename('fullpath'));
    if nargin < 1 || isempty(matfile)
        matfile = fullfile(thisdir, 'results', 'ber_ce_ici_results.mat');
    end
    assert(isfile(matfile), ...
        'Results file not found: %s\nRun ber_ce_ici first.', matfile);

    S = load(matfile);

    % ---- required fields ----
    SNR_dbs        = S.SNR_dbs;
    SNRp_dB        = S.SNRp_dB;
    BER_ofdm_ideal = S.BER_ofdm_ideal;
    BER_ofdm_ce    = S.BER_ofdm_ce;
    BER_otfs_ideal = S.BER_otfs_ideal;
    BER_otfs_ce    = S.BER_otfs_ce;         % nSNR x nP

    % ---- optional fields (older result files may lack these) ----
    evaname      = getfield_default(S, 'evaname',      'EVA');
    SNRp_ofdm_dB = getfield_default(S, 'SNRp_ofdm_dB', NaN);
    nMonteCarlo  = getfield_default(S, 'nMonteCarlo',  NaN);
    Ngrid        = getfield_default(S, 'N',            NaN);
    Mgrid        = getfield_default(S, 'M',            NaN);

    nP = numel(SNRp_dB);

    figure; hold on; grid on;
    semilogy(SNR_dbs, BER_ofdm_ideal, '-d', 'LineWidth', 2, 'MarkerSize', 8);
    semilogy(SNR_dbs, BER_ofdm_ce,    '-x', 'LineWidth', 2, 'MarkerSize', 8);
    semilogy(SNR_dbs, BER_otfs_ideal, '-o', 'LineWidth', 2, 'MarkerSize', 8);
    mk = {'--s','--^','--v','--*','--p'};
    for pj = 1:nP
        semilogy(SNR_dbs, BER_otfs_ce(:,pj), mk{mod(pj-1,numel(mk))+1}, ...
                 'LineWidth', 2, 'MarkerSize', 8);
    end
    set(gca, 'YScale', 'log');
    xlabel('SNR_d [dB]');  ylabel('BER');

    if isnan(SNRp_ofdm_dB)
        ofdm_ce_lbl = 'OFDM (CE)';
    else
        ofdm_ce_lbl = sprintf('OFDM (CE, SNR_p=%d dB)', SNRp_ofdm_dB);
    end
    ce_lbl = arrayfun(@(x) sprintf('OTFS (CE, SNR_p=%d dB)', x), SNRp_dB, ...
                      'UniformOutput', false);
    legend([{'OFDM (ideal, banded MMSE)', ofdm_ce_lbl, 'OTFS (ideal, MP)'}, ce_lbl], ...
           'Location', 'southwest');

    ttl = sprintf('%s ICI channel — CE + demodulation BER', evaname);
    if ~isnan(Ngrid)
        ttl = sprintf('%s  (N=%d, M=%d, %d trials)', ttl, Ngrid, Mgrid, nMonteCarlo);
    end
    title(ttl);

    % ---- console recap ----
    fprintf('Loaded %s\n', matfile);
    fprintf('%6s  %10s %10s %10s', 'SNR_d', 'OFDMid', 'OFDMce', 'OTFSid');
    for pj = 1:nP, fprintf('  OTFSce%d', SNRp_dB(pj)); end
    fprintf('\n');
    for si = 1:numel(SNR_dbs)
        fprintf('%5d   %10.3e %10.3e %10.3e', SNR_dbs(si), ...
            BER_ofdm_ideal(si), BER_ofdm_ce(si), BER_otfs_ideal(si));
        for pj = 1:nP, fprintf(' %10.3e', BER_otfs_ce(si,pj)); end
        fprintf('\n');
    end
end

function v = getfield_default(S, name, dflt)
    if isfield(S, name), v = S.(name); else, v = dflt; end
end

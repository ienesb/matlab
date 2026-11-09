function hf = plot_ber_ce_channels(matfile)
% PLOT_BER_CE_CHANNELS  BER vs SNR for OFDM-ICI and CP-OTFS-ICI across channels,
% from ber_ce_channels.m output.  Replots WITHOUT re-running the simulation.
%
%   plot_ber_ce_channels                    % default results file
%   plot_ber_ce_channels(matfile)          % a specific results file
%   hf = plot_ber_ce_channels(...)          % handle to the figure
%
% Zero-BER points (no bit errors observed in the whole sweep) cannot be shown on
% a log axis; they are masked out and reported on the console, with the
% measurement floor 1/(total bits) drawn as a reference when available.
    thisdir = fileparts(mfilename('fullpath'));
    if nargin < 1 || isempty(matfile)
        matfile = fullfile(thisdir, 'results', 'ber_ce_channels_results.mat');
    end
    assert(isfile(matfile), ...
        'Results file not found: %s\nRun ber_ce_channels first.', matfile);
    S = load(matfile);

    SNR = S.SNR_dbs;  nSNR = numel(SNR);
    nC  = size(S.configs,1);  nCh = numel(S.channel_set);
    % Default palette: OFDM variants in warm tones, OTFS in blue.  Falls back to
    % `lines` if a run defines more configs than we have hand-picked colours.
    cols = [0.95 0.55 0.35;    % OFDM (LS+interp)  -- pale: the weak CE
            0.85 0.33 0.10;    % OFDM (DFT CE)     -- strong: the improved CE
            0.00 0.45 0.74];   % CP-OTFS (MP)      -- blue
    if nC > size(cols,1), cols = lines(nC); end
    cols = cols(1:nC,:);
    styl = {'-','--',':','-.'};
    mrk  = {'o','^','s','d','v','p'};

    % ---- console summary + zero-BER report ---------------------------------
    fprintf('\nBER vs SNR  (%d-QAM, SNR_p=%d dB, %d trials/point)\n', ...
            S.M_order, S.SNRp_dB, S.nTrials);
    for chi = 1:nCh
        fprintf(' [%s]\n', S.channel_set{chi});
        for ci = 1:nC
            z = SNR(S.ber(:,ci,chi) == 0);
            fprintf('   %-12s BER: %s\n', S.configs{ci,1}, ...
                    strjoin(arrayfun(@(x) sprintf('%.1e',x), S.ber(:,ci,chi).', ...
                            'UniformOutput',false), ' '));
            if ~isempty(z)
                fprintf('     (no bit errors observed at SNR = %s dB -> masked on the log axis)\n', ...
                        num2str(z));
            end
        end
    end

    % ---- figure -------------------------------------------------------------
    hf = figure('Name','BER vs SNR (realistic CE + demod)', ...
                'Position',[200 200 840 560],'Color','w');
    lg = cell(1, nC*nCh);  k = 0;
    for chi = 1:nCh
        for ci = 1:nC
            y = S.ber(:,ci,chi);
            y(y == 0) = NaN;                 % zeros cannot be drawn on a log axis
            semilogy(SNR, y, [styl{chi} mrk{ci}], 'Color', cols(ci,:), ...
                     'LineWidth',1.9, 'MarkerSize',6); hold on;
            k = k + 1;
            lg{k} = sprintf('%s \\times %s', S.channel_set{chi}, S.configs{ci,1});
        end
    end
    % Measurement floor 1/(total bits), when the cell counts were saved.
    if isfield(S,'n_data_bits')
        yfloor = 1 ./ (S.nTrials * S.n_data_bits);
        yline(min(yfloor), ':', 'measurement floor', 'Color',[0.4 0.4 0.4], ...
              'LineWidth',1.2, 'LabelHorizontalAlignment','left');
    end
    hold off; grid on; box on;
    xlabel('SNR  [dB]');  ylabel('BER');
    title(sprintf(['BER vs SNR: realistic CE + demod  (%d-QAM, SNR_p=%d dB, %d trials)\n' ...
                   'OFDM: pilot LS+interp -> 1-tap MMSE   |   OTFS: embedded pilot -> MP'], ...
                  S.M_order, S.SNRp_dB, S.nTrials));
    legend(lg, 'Location','southwest');  set(gca,'FontSize',11);
    if exist('theme','file'), theme(gcf,"light"); end

    % PNG name derived from the results file so distinct runs (e.g. _nodopp)
    % don't overwrite each other's figure.
    [resdir, base] = fileparts(matfile);
    pngname = [regexprep(base, '_results$', '') '.png'];
    exportgraphics(hf, fullfile(resdir, pngname), ...
                   'Resolution',150, 'BackgroundColor','white');
    fprintf('\nSaved %s in %s\n', pngname, resdir);
end

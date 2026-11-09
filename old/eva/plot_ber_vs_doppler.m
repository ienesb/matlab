function hf = plot_ber_vs_doppler(matfile)
% PLOT_BER_VS_DOPPLER  BER vs mobility (velocity -> nu_max) at fixed SNR, from
% ber_vs_doppler.m output.  Replots WITHOUT re-running the sweep.
%
%   plot_ber_vs_doppler                 % default results file
%   plot_ber_vs_doppler(matfile)        % a specific file
%   hf = plot_ber_vs_doppler(...)       % figure handle
%
% v = 0 is ICI-free by construction, so any OFDM-vs-OTFS gap there is pure
% DIVERSITY (uncoded OFDM = order 1), not Doppler.  Zero-BER points are masked
% on the log axis and reported on the console.
    thisdir = fileparts(mfilename('fullpath'));
    if nargin < 1 || isempty(matfile)
        matfile = fullfile(thisdir, 'results', 'ber_vs_doppler_results.mat');
    end
    assert(isfile(matfile), 'Results file not found: %s\nRun ber_vs_doppler first.', matfile);
    S = load(matfile);

    v = S.v_kmh(:);  SNR = S.SNR_dbs;  nS = numel(SNR);
    nC = size(S.configs,1);  names = S.configs(:,1);
    cols  = [0.85 0.33 0.10;    % OFDM (red-orange)
             0.00 0.45 0.74];   % OTFS (blue)
    styl  = {'-','--',':','-.'};        % one line style per SNR

    fprintf('\nBER vs mobility  (%d-QAM, %d trials/point)\n', S.M_order, S.nTrials);
    fprintf('%-14s', 'v [km/h] ->'); fprintf('%8.0f', v); fprintf('\n');
    for si = 1:nS
        for ci = 1:nC
            fprintf('%s @%2ddB : ', pad(names{ci},12), SNR(si));
            fprintf('%8.1e', S.ber(:,ci,si)); fprintf('\n');
        end
    end

    hf = figure('Name','BER vs Doppler','Position',[200 200 820 560],'Color','w');
    lg = {};
    for si = 1:nS
        for ci = 1:nC
            y = S.ber(:,ci,si);  y(y==0) = NaN;         % zeros off the log axis
            semilogy(v, y, [styl{si} 'o'], 'Color', cols(ci,:), ...
                     'LineWidth',1.9, 'MarkerSize',6, 'MarkerFaceColor', cols(ci,:)); hold on;
            lg{end+1} = sprintf('%s, SNR=%d dB', names{ci}, SNR(si)); %#ok<AGROW>
        end
    end
    hold off; grid on; box on;
    xlabel('velocity  [km/h]');  ylabel('BER');
    title(sprintf(['BER vs mobility, EVA grid, pilot lattice frozen (%d-QAM, %d trials)\n' ...
                   'v=0 is ICI-free: any gap there is pure DIVERSITY, not Doppler'], ...
                   S.M_order, S.nTrials));
    legend(lg,'Location','southeast','FontSize',9); set(gca,'FontSize',11);
    if exist('theme','file'), theme(hf,"light"); end

    png = fullfile(fileparts(matfile),'ber_vs_doppler.png');
    exportgraphics(hf, png, 'Resolution',150, 'BackgroundColor','white');
    fprintf('\nSaved %s\n', png);
end

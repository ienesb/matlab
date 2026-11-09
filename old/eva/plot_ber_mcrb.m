function hf = plot_ber_mcrb(matfile)
% PLOT_BER_MCRB  Detection BER and the sensing bounds, side by side vs SNR.
%
% Three stacked panels sharing the SNR axis, from mcrb_ce_ici.m output:
%   (1) detector BER (and SER, dotted) -- what the demodulator achieves
%   (2) delay   RMSE bound: MCRB (solid) vs genie CRB (dashed)
%   (3) Doppler RMSE bound: MCRB (solid) vs genie CRB (dashed)
%
% The point of stacking them: BER and MCRB are NOT the same currency.  The BER
% panel can show orders of magnitude between the waveforms while the bound panels
% show them on top of each other -- sensing integrates the whole frame and is
% noise-limited (N0*G) once the symbol-error rate is small enough that the
% pseudo-true bias BB stops mattering.
%
%   plot_ber_mcrb                 % default results file
%   plot_ber_mcrb(matfile)        % a specific results file
    thisdir = fileparts(mfilename('fullpath'));
    if nargin < 1 || isempty(matfile)
        matfile = fullfile(thisdir, 'results', 'mcrb_ce_ici_results.mat');
    end
    assert(isfile(matfile), 'Results file not found: %s\nRun mcrb_ce_ici first.', matfile);
    S = load(matfile);

    SNR = S.SNR_dbs;  nSNR = numel(SNR);  nC = size(S.configs,1);
    it = S.idx_tau;   in = S.idx_nu;
    names = S.configs(:,1);
    col = lines(nC);

    % ---- bounds ------------------------------------------------------------
    [rt_m, rt_c, rn_m, rn_c] = deal(zeros(nSNR,nC));
    for ci = 1:nC
        for si = 1:nSNR
            N0  = 1/db2pow(SNR(si));
            LB  = N0*S.G_bar(:,:,si,ci) + S.F_bar(:,:,si,ci) + S.BB_bar(:,:,si,ci);
            CRB = N0*S.Gc_bar(:,:,si,ci);
            rt_m(si,ci)=sqrt(sum(real(diag(LB (it,it)))));
            rt_c(si,ci)=sqrt(sum(real(diag(CRB(it,it)))));
            rn_m(si,ci)=sqrt(sum(real(diag(LB (in,in)))));
            rn_c(si,ci)=sqrt(sum(real(diag(CRB(in,in)))));
        end
    end
    has_ser = isfield(S,'ser_data');

    hf = figure('Position',[100 60 820 900],'Color','w');

    % ---- (1) BER -----------------------------------------------------------
    subplot(3,1,1); hold on; grid on; box on;
    lg1 = {};
    for ci = 1:nC
        plot(SNR, nz(S.ber_data(:,ci)), '-o', 'Color',col(ci,:),'LineWidth',2,'MarkerSize',5);
        lg1{end+1} = sprintf('%s: BER', names{ci}); %#ok<AGROW>
    end
    if has_ser
        for ci = 1:nC
            plot(SNR, nz(S.ser_data(:,ci)), ':', 'Color',col(ci,:),'LineWidth',1.4);
            lg1{end+1} = sprintf('%s: SER', names{ci}); %#ok<AGROW>
        end
    end
    set(gca,'YScale','log'); ylabel('detection BER / SER');
    title(sprintf('%s ICI — detection vs sensing  (P=%d, %d trials, OFDM CE = %s)', ...
                  S.evaname, S.P, S.nTrials, getfield_default(S,'ce_mode','?')));
    legend(lg1,'Location','southwest','FontSize',8); set(gca,'FontSize',10);

    % ---- (2) delay ---------------------------------------------------------
    subplot(3,1,2); hold on; grid on; box on;
    lg2 = {};
    for ci = 1:nC
        plot(SNR, rt_m(:,ci)*1e9, '-o', 'Color',col(ci,:),'LineWidth',2,'MarkerSize',5);
        plot(SNR, rt_c(:,ci)*1e9, '--', 'Color',col(ci,:),'LineWidth',1.5);
        lg2 = [lg2, {sprintf('%s: MCRB',names{ci}), sprintf('%s: CRB',names{ci})}]; %#ok<AGROW>
    end
    set(gca,'YScale','log'); ylabel('delay RMSE bound [ns]');
    legend(lg2,'Location','southwest','FontSize',8); set(gca,'FontSize',10);

    % ---- (3) Doppler -------------------------------------------------------
    subplot(3,1,3); hold on; grid on; box on;
    for ci = 1:nC
        plot(SNR, rn_m(:,ci), '-o', 'Color',col(ci,:),'LineWidth',2,'MarkerSize',5);
        plot(SNR, rn_c(:,ci), '--', 'Color',col(ci,:),'LineWidth',1.5);
    end
    set(gca,'YScale','log'); xlabel('SNR [dB]'); ylabel('Doppler RMSE bound [Hz]');
    legend(lg2,'Location','southwest','FontSize',8); set(gca,'FontSize',10);

    if exist('theme','file'), theme(hf,"light"); end
    resdir = fileparts(matfile);
    png = fullfile(resdir,'mcrb_ce_ici_ber_mcrb.png');
    exportgraphics(hf, png, 'Resolution',150, 'BackgroundColor','white');

    % ---- recap -------------------------------------------------------------
    fprintf('\n%5s |', 'SNR');
    for ci=1:nC, fprintf(' %-11s BER    RMSEtau[ns] |', names{ci}); end
    fprintf('\n');
    for si = 1:nSNR
        fprintf('%5d |', SNR(si));
        for ci=1:nC, fprintf(' %19.3e %10.3e |', S.ber_data(si,ci), rt_m(si,ci)*1e9); end
        fprintf('\n');
    end
    fprintf('\nSaved %s\n', png);
end

function v = getfield_default(S, f, d)
    if isfield(S,f), v = S.(f); else, v = d; end
end

function y = nz(y)
    y(y == 0) = NaN;   % zeros cannot be drawn on a log axis
end

% plot_crlb.m
% --------------------------------------------------------------------------
% Loads the multi-target CRBs saved by run_fim_avg.m and plots the delay
% CRB (trace over paths and path-1, in s^2 or range m^2) versus SNR for the
% waveform family {CP-OFDM, RCP-OFDM, CP-OTFS, RCP-OTFS} x {ICI, no-ICI},
% together with the continuous-time floor and the uniform-spectrum
% reference.  When the saved run jointly estimated Doppler (e.g. channel
% 'paperv2'), the analogous Doppler CRB figures (Hz^2) are added.
%
% The FIM is linear in 1/sigma^2 so the CRB is linear in sigma^2.  We sweep
% SNR analytically from the saved sigma^2 = 1 results:
%   SNR := |alpha_1|^2 / sigma^2,  CRB(SNR) = CRB(sigma^2=1) * sigma^2(SNR).
% --------------------------------------------------------------------------

clear; clc; close all;
thisdir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisdir, 'utils'));
resdir = fullfile(thisdir, 'results');

% ===== Units switch =====================================================
units = 'range';      % 'delay' (s^2)  or  'range' (m^2)

S = load(fullfile(resdir, 'eva_fim_avg_results.mat'));
P = S.P;

c_light    = 3e8;
alpha_ref2 = abs(S.alpha_vec(1))^2;

snr_db   = -10:2:30;
snr_lin  = 10.^(snr_db/10);
sigma2_s = alpha_ref2 ./ snr_lin;

% ===== Curve family ======================================================
% {field suffix in the results struct, legend label, line style}
% Older results files may lack RCP-OFDM; those curves are skipped.
curves = { ...
    'OF_I',    'CP-OFDM (ICI)',      'b-o' ; ...
    'RCPOF_I', 'RCP-OFDM (ICI)',     'c-x' ; ...
    'CPOT_I',  'CP-OTFS (ICI)',      'r-s' ; ...
    'RCPOT_I', 'RCP-OTFS (ICI)',     'm-d' ; ...
    'OF_N',    'CP-OFDM (no-ICI)',   'b--^' ; ...
    'RCPOF_N', 'RCP-OFDM (no-ICI)',  'c--*' ; ...
    'CPOT_N',  'CP-OTFS (no-ICI)',   'r--v' ; ...
    'RCPOT_N', 'RCP-OTFS (no-ICI)',  'm--p' };
keep   = cellfun(@(sfx) isfield(S, ['CRBt_' sfx]), curves(:,1));
curves = curves(keep, :);
nC     = size(curves, 1);

% Per-curve statistics at sigma^2 = 1: trace over paths and path-1 value.
trace_tau = zeros(1,nC);  crb1_tau = zeros(1,nC);
trace_nu  = zeros(1,nC);  crb1_nu  = zeros(1,nC);
for k = 1:nC
    Ct = S.(['CRBt_' curves{k,1}]);
    Cn = S.(['CRBn_' curves{k,1}]);
    trace_tau(k) = sum(mean(Ct, 1));   crb1_tau(k) = mean(Ct(:,1));
    trace_nu(k)  = sum(mean(Cn, 1));   crb1_nu(k)  = mean(Cn(:,1));
end

% ===== Continuous-time CRLB benchmark (ct_crlb_benchmark.tex) ============
% Modulation- and waveform-independent lower bound for a single target of
% gain |alpha_1|^2, at the frame's occupied bandwidth W, observation time T,
% and full-frame energy SNR.  Linear in sigma^2, evaluated across the sweep.
%
%   Fundamental floor (band-edge loading) : CRLB(tau) = 1 / (2 pi^2 W^2 SNR)
%   Uniform-spectrum reference (flat load): CRLB(tau) = 3 / (2 pi^2 W^2 SNR)
W_bench   = S.N_o * S.deltaf_o;        % bandwidth  [Hz]
T_bench   = S.M_o * S.Ts_o;            % duration   [s]
NM_bench  = S.N_o * S.M_o;             % frame size (degrees of freedom)
snr_frame = NM_bench * snr_lin;        % |alpha_1|^2 NM / sigma^2

crb_tau_floor = 1 ./ (2*pi^2 * W_bench^2 .* snr_frame);   % absolute floor [s^2]
crb_tau_flat  = 3 * crb_tau_floor;                        % uniform ref    [s^2]
crb_nu_floor  = 1 ./ (2*pi^2 * T_bench^2 .* snr_frame);   % Doppler floor  [Hz^2]
crb_nu_flat   = 3 * crb_nu_floor;                         % Doppler unif.  [Hz^2]

% Trace benchmark: per-path floor summed over all P paths.
sum_inv_pow     = sum(1 ./ abs(S.alpha_vec).^2);
trace_tau_floor = sum_inv_pow / (2*pi^2 * W_bench^2 * NM_bench) .* sigma2_s;
trace_tau_flat  = 3 * trace_tau_floor;
trace_nu_floor  = sum_inv_pow / (2*pi^2 * T_bench^2 * NM_bench) .* sigma2_s;
trace_nu_flat   = 3 * trace_nu_floor;

% Unit conversion: tau [s] -> R [m] via R = c/2 * tau, CRB_R = (c/2)^2 CRB_tau.
switch lower(units)
    case 'delay'
        u_scale = 1;
        yl_tr   = 'tr\{CRB_{\tau}\}  =  \Sigma_p CRB(\tau_p)  [s^2]';
        yl_1    = 'CRB(\tau_1)  [s^2]';
        ttl_tr  = 'Trace of delay CRB';    ttl_1 = 'Delay CRB';
        png_tr  = 'crlb_trace_tau.png';    png_1 = 'crlb_tau1.png';
    case 'range'
        u_scale = (c_light/2)^2;
        yl_tr   = 'tr\{CRB_R\}  =  (c/2)^2 \Sigma_p CRB(\tau_p)  [m^2]';
        yl_1    = '(c/2)^2 CRB(\tau_1)  [m^2]';
        ttl_tr  = 'Trace of range CRB';    ttl_1 = 'Range CRB';
        png_tr  = 'crlb_trace_range.png';  png_1 = 'crlb_range1.png';
    otherwise
        error('Unknown units "%s" (use ''delay'' or ''range'')', units);
end

bench_leg = {'CT floor (\Sigma_p, edge)', 'Uniform-spectrum ref (3\times floor)'};
cfg = sprintf('(params: %s, channel: %s)', S.param_mode, S.channel_type);

% ---- Figure helper: values (1 x nC) at sigma2=1, benchmarks, labels -----
make_fig = @(vals, bfloor, bflat, ylab, ttl, logscale, png) ...
    plot_curve_family(snr_db, sigma2_s, vals, curves, bfloor, bflat, ...
                      ylab, ttl, bench_leg, logscale, ...
                      ternary(isempty(png), '', fullfile(resdir, png)));

% ===== Delay CRB figures ==================================================
make_fig(u_scale*trace_tau, u_scale*trace_tau_floor, u_scale*trace_tau_flat, ...
         yl_tr, sprintf('%s across P=%d paths  %s', ttl_tr, P, cfg), true, png_tr);
print_ranking('Delay CRB ordering (trace over paths, best -> worst)', ...
              curves, trace_tau, 'tr CRB(tau)', 's^2');
make_fig(u_scale*trace_tau, u_scale*trace_tau_floor, u_scale*trace_tau_flat, ...
         yl_tr, sprintf('%s across P=%d paths  %s', ttl_tr, P, cfg), false, '');
make_fig(u_scale*crb1_tau, u_scale*crb_tau_floor, u_scale*crb_tau_flat, ...
         yl_1, sprintf('%s for path 1 (\\tau_1 = %.1f ns)  %s', ...
                       ttl_1, S.tau_vec(1)*1e9, cfg), true, png_1);
fprintf('Saved %s and %s in %s\n', png_tr, png_1, resdir);

% ===== Doppler CRB figures (only when nu was jointly estimated) ==========
if isfield(S, 'estimate_doppler') && S.estimate_doppler
    make_fig(trace_nu, trace_nu_floor, trace_nu_flat, ...
             'tr\{CRB_{\nu}\}  =  \Sigma_p CRB(\nu_p)  [Hz^2]', ...
             sprintf('Trace of Doppler CRB across P=%d paths  %s', P, cfg), ...
             true, 'crlb_trace_nu.png');
    print_ranking('Doppler CRB ordering (trace over paths, best -> worst)', ...
                  curves, trace_nu, 'tr CRB(nu)', 'Hz^2');
    make_fig(crb1_nu, crb_nu_floor, crb_nu_flat, ...
             'CRB(\nu_1)  [Hz^2]', ...
             sprintf('Doppler CRB for path 1 (\\nu_1 = %+.1f Hz)  %s', ...
                     S.nu_vec(1), cfg), true, 'crlb_nu1.png');
    fprintf('Saved crlb_trace_nu.png and crlb_nu1.png in %s\n', resdir);
end

% ===== Paper precoded-waveform figure (Fig. 1 reproduction) =============
% Single-block ranging CRB for OFDM / SC / OTFS under the paper's decoupled
% model (delay steering over the full N; waveform enters only through the
% frequency power |Q^H s|^2).  OFDM sits on the Jensen bound for constant-
% modulus (PSK) constellations; QAM lifts it slightly above.
if isfield(S,'do_precoded') && S.do_precoded
    snr_p  = -10:1:6;                          % paper Fig. 1 SNR axis (= 1/sigma^2)
    s2_p   = alpha_ref2 ./ 10.^(snr_p/10);
    trc    = @(C) sum(mean(C,1));              % trace over paths at sigma^2 = 1
    y_ofdm = u_scale * trc(S.CRBt_POF) * s2_p;
    y_sc   = u_scale * trc(S.CRBt_PSC) * s2_p;
    y_otfs = u_scale * trc(S.CRBt_POT) * s2_p;
    y_jen  = u_scale * sum(S.CRBt_PJEN) * s2_p;
    figure('Name','Paper Fig.1 reproduction','Position',[220 200 720 520]);
    plot(snr_p, y_ofdm, 'o-', 'Color',[0.85 0.33 0.10],'LineWidth',1.8,'MarkerSize',6); hold on;
    plot(snr_p, y_sc,   's-', 'Color',[0.00 0.45 0.74],'LineWidth',1.8,'MarkerSize',6);
    plot(snr_p, y_otfs, 'd--','Color',[0.93 0.69 0.13],'LineWidth',1.8,'MarkerSize',6);
    plot(snr_p, y_jen,  'kx', 'LineWidth',1.5,'MarkerSize',9);
    hold off; grid on; box on;
    xlabel('SNR = 1/\sigma^2  [dB]'); ylabel('Ranging CRB  [m^2]');
    title(sprintf('Paper Fig. 1 reproduction  (%s, N=%d, L=%d)', S.modulation, S.N_o, P));
    legend({'OFDM','SC','OTFS','Jensen bound'}, 'Location','northeast');
    theme(gcf,"light");
    set(gca, 'FontSize', 11);
    exportgraphics(gcf, fullfile(resdir,'crlb_paper_fig1.png'), 'Resolution',150, 'BackgroundColor','white');
    fprintf('Saved crlb_paper_fig1.png in %s\n', resdir);
    fprintf(['Precoded trace CRB (sigma^2=1):  OFDM=%.4e  SC=%.4e  OTFS=%.4e  Jensen=%.4e s^2\n' ...
             '  OTFS vs OFDM: +%.1f%% (%.2f dB) ;  SC vs OFDM: +%.1f%% (%.2f dB)\n'], ...
            trc(S.CRBt_POF), trc(S.CRBt_PSC), trc(S.CRBt_POT), sum(S.CRBt_PJEN), ...
            100*(trc(S.CRBt_POT)/trc(S.CRBt_POF)-1), 10*log10(trc(S.CRBt_POT)/trc(S.CRBt_POF)), ...
            100*(trc(S.CRBt_PSC)/trc(S.CRBt_POF)-1), 10*log10(trc(S.CRBt_PSC)/trc(S.CRBt_POF)));
end

% ------------------------------------------------------------------------
function plot_curve_family(snr_db, sigma2_s, vals, curves, bfloor, bflat, ...
                           ylab, ttl, bench_leg, logscale, png_path)
% One CRB-vs-SNR figure: nC waveform curves (vals scaled by sigma2 sweep)
% plus the two continuous-time benchmark lines.  Curves are drawn (and
% therefore listed in the legend) sorted best (lowest CRB) to worst.
    [vals, ord] = sort(vals, 'ascend');
    curves      = curves(ord, :);
    nC = size(curves, 1);
    figure('Name', ttl, 'Position', [200, 240, 760, 500]);
    for k = 1:nC
        y = vals(k) * sigma2_s;
        if logscale, semilogy(snr_db, y, curves{k,3}, 'LineWidth', 1.6, 'MarkerSize', 6);
        else,        plot    (snr_db, y, curves{k,3}, 'LineWidth', 1.6, 'MarkerSize', 6);
        end
        hold on;
    end
    hold off; grid on;
    xlabel('SNR [dB]'); ylabel(ylab); title(ttl);
    legend(curves(:,2).', 'Location', 'northeast');
    theme(gcf,"light");
    set(gca, 'FontSize', 11);
    if ~isempty(png_path)
        exportgraphics(gcf, png_path, 'Resolution', 150, 'BackgroundColor', 'white');
    end
end

function print_ranking(header, curves, vals, qty, unit)
% Console ranking of the curves by their sigma^2 = 1 value (lowest = best),
% with the percentage penalty relative to the best curve.
    [vs, ord] = sort(vals, 'ascend');
    fprintf('\n%s:\n', header);
    for r = 1:numel(vs)
        fprintf('  %d. %-20s %s = %.4e %s   (+%.2f%% vs best)\n', ...
                r, curves{ord(r), 2}, qty, vs(r), unit, 100*(vs(r)/vs(1) - 1));
    end
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end

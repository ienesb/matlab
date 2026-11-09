% ber_vs_doppler.m
% --------------------------------------------------------------------------
% What happens to the OFDM BER floor if we reduce the maximum Doppler shift?
%
% Sweeps the mobility (velocity -> nu_max) at FIXED SNR on the EVA grid, for
%   OFDM (DFT CE -> 1-tap MMSE)   and   CP-OTFS (embedded pilot -> MP).
%
% The pilot lattice is FROZEN at the nominal v=120 km/h design (k_nu=4) so the
% pilot overhead is identical at every velocity -- the ONLY thing changing is
% the Doppler.  Otherwise match_pilot_spacing would re-tune dM_p as k_nu shrinks
% and confound the comparison.
%
% Expectation being tested:
%   * OFDM's floor is ICI-driven, ICI ~ (pi*nu_max*T)^2/3, so the floor should
%     fall as nu_max^2 and vanish at low mobility (v -> 0 is a pure, ICI-free
%     frequency-selective Rayleigh channel).
%   * OTFS should be almost flat in velocity: its MP detector models the Doppler
%     taps explicitly, so Doppler is not an impairment for it.
%   * v = 0 isolates DIVERSITY from Doppler: with no ICI at all, any residual
%     OFDM-vs-OTFS gap is pure diversity (uncoded OFDM subcarriers see a flat
%     fade -> diversity order 1 -> BER ~ 1/(4*SNR); OTFS spreads over the whole
%     DD grid -> full diversity).
% --------------------------------------------------------------------------
clear; clc; close all;
thisdir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisdir, 'utils'));
resdir = fullfile(thisdir, 'results');
warning('off','MATLAB:nearlySingularMatrix');

%% ---- toggles ---------------------------------------------------------------
v_kmh    = [0 3 10 30 60 120 240];      % mobility sweep
SNR_dbs  = [20 30];                     % fixed SNRs (floor region)
SNRp_dB  = 50;   nTrials = 20;   seed_base = 9200;
n_iter_mp = 30;  damping = 0.3;  M_order = 4;
results_name = 'ber_vs_doppler_results.mat';

nbits = log2(M_order);
constellation = qammod((0:M_order-1).', M_order, UnitAveragePower=true);
configs = {'OFDM (DFT CE)','ofdm'; 'CP-OTFS (MP)','otfs'};
nC = size(configs,1);  nV = numel(v_kmh);  nS = numel(SNR_dbs);

%% ---- EVA grid --------------------------------------------------------------
N = 512;  M = 128;  deltaf = 15e3;  Ts = 1/deltaf;  T = 1/deltaf;
fc = 4e9;  c = 3e8;
delay_taps = round([0 150 370]*1e-9 * N * deltaf);
amp = (10.^([0 -1.4 -0.6]/20)).';  amp = amp/norm(amp);
tau = (delay_taps/(N*deltaf)).';  P = numel(tau);  l_tau = max(delay_taps);

% Pilot lattice FROZEN at the nominal v=120 km/h design.
nu_max_nom = (120/3.6)*fc/c;
k_nu_nom   = ceil(nu_max_nom * M / deltaf);
l_p = l_tau;  k_p = M/2;
guard_rows = 1:(2*l_tau+1);  ce_rows = (l_p+1):(2*l_tau+1);
data_mask_otfs = true(N,M);  data_mask_otfs(guard_rows,:) = false;
otfs_nondata = numel(guard_rows)*M;
[dN_p, dM_p, n_pilot_cells] = match_pilot_spacing(N, M, otfs_nondata, l_tau, k_nu_nom);
n_pilots = 1:dN_p:N;  m_pilots = 1:dM_p:M;
pilot_mask_ofdm = false(N,M);  pilot_mask_ofdm(n_pilots,m_pilots) = true;
data_mask_ofdm  = ~pilot_mask_ofdm;

fprintf('EVA grid, pilot lattice frozen at v=120 km/h (k_nu=%d, spacing %d x %d)\n', ...
        k_nu_nom, dN_p, dM_p);
fprintf('v [km/h] -> nu_max [Hz] -> nu_max*T -> predicted ICI [dB]\n');
for vi = 1:nV
    nm = (v_kmh(vi)/3.6)*fc/c;
    fprintf('  %5.0f -> %7.1f -> %8.2e -> %6.1f\n', v_kmh(vi), nm, nm*T, ...
            10*log10(max((pi*nm*T)^2/3, eps)));
end

ber = nan(nV, nC, nS);  tAll = tic;

for si = 1:nS
    N0  = 1/db2pow(SNR_dbs(si));
    x_p = sqrt(db2pow(SNRp_dB) * N0);
    a_p = sqrt(db2pow(SNRp_dB) * N0 / n_pilot_cells);
    for vi = 1:nV
        nu_max = (v_kmh(vi)/3.6)*fc/c;
        for ci = 1:nC
            modn = configs{ci,2};
            be_t = zeros(1,nTrials);
            parfor t = 1:nTrials
                warning('off','MATLAB:nearlySingularMatrix');
                rs = RandStream('Threefry','Seed',seed_base);
                rs.Substream = 1000*vi + t;      % paired across configs & SNR
                RandStream.setGlobalStream(rs);

                nu    = (nu_max * cos(2*pi*rand(1,P) - pi)).';
                alpha = amp .* exp(1j*2*pi*rand(P,1));
                sym = randi(M_order, N, M) - 1;
                Xd  = qammod(sym, M_order, UnitAveragePower=true);
                Z   = sqrt(N0/2)*(randn(N,M) + 1j*randn(N,M));

                if strcmp(modn,'ofdm')
                    d = Xd;  d(pilot_mask_ofdm) = a_p;
                    Y = ofdm_mean(d, deltaf, Ts, alpha, tau, nu, true) + Z;
                    % DFT CE: LS fit of the (l_tau+1) delay taps across the comb
                    H_ls  = Y(n_pilots,m_pilots) ./ d(n_pilots,m_pilots);
                    np0   = (n_pilots(:)-1);
                    Abas  = exp(-1j*2*pi * np0       * (0:l_tau) / N);
                    Afull = exp(-1j*2*pi * (0:N-1).' * (0:l_tau) / N);
                    Hp    = Afull * (Abas \ H_ls);
                    F_H   = griddedInterpolant({(1:N).', m_pilots}, Hp, 'linear','nearest');
                    [Nq,Mq] = ndgrid(1:N,1:M);
                    H_hat = F_H(Nq,Mq);
                    Xsoft = conj(H_hat).*Y ./ (abs(H_hat).^2 + N0);
                    sh = qamdemod(Xsoft(data_mask_ofdm), M_order, UnitAveragePower=true);
                    be_t(t) = mean(int2bit(sym(data_mask_ofdm),nbits) ~= int2bit(sh,nbits),'all');
                else
                    d = zeros(N,M);  d(data_mask_otfs) = Xd(data_mask_otfs);
                    d(l_p+1,k_p+1) = x_p;
                    Gtrue = build_effective_dd_eva(alpha, delay_taps, nu, N, M, deltaf, Inf);
                    Y = apply_dd_channel_eva(d, Gtrue) + Z;
                    h_ce = zeros(N,M);
                    h_ce(ce_rows,:) = Y(ce_rows,:) / x_p;
                    h_ce(abs(h_ce) < 3*sqrt(N0)/x_p) = 0;
                    [rc,cc,gc] = find(h_ce);
                    if isempty(gc)
                        be_t(t) = 0.5;
                    else
                        dl_c = mod(rc-(l_p+1),N);  dk_c = mod(cc-(k_p+1),M);
                        known = zeros(N,M);  known(l_p+1,k_p+1) = x_p;
                        xhat = mp_detector(Y, dl_c, dk_c, gc, N0, constellation, ...
                                           n_iter_mp, data_mask_otfs, known, damping);
                        sh = qamdemod(xhat(data_mask_otfs), M_order, UnitAveragePower=true);
                        be_t(t) = mean(int2bit(sym(data_mask_otfs),nbits) ~= int2bit(sh,nbits),'all');
                    end
                end
            end
            ber(vi,ci,si) = mean(be_t);
        end
        fprintf('SNR=%2d dB | v=%3.0f km/h (nu_max=%6.1f Hz) | OFDM=%.3e  OTFS=%.3e   (%.0fs)\n', ...
                SNR_dbs(si), v_kmh(vi), nu_max, ber(vi,1,si), ber(vi,2,si), toc(tAll));
    end
end

save(fullfile(resdir,results_name),'ber','v_kmh','SNR_dbs','configs','nTrials','M_order','SNRp_dB');

%% ---- figure ----------------------------------------------------------------
figure('Name','BER vs Doppler','Position',[200 200 820 560],'Color','w');
cols = [0.85 0.33 0.10; 0.00 0.45 0.74];  styl = {'-','--'};
lg = {};
for si = 1:nS
    for ci = 1:nC
        y = ber(:,ci,si);  y(y==0) = NaN;
        semilogy(v_kmh, y, [styl{si} 'o'], 'Color', cols(ci,:), ...
                 'LineWidth',1.9,'MarkerSize',6,'MarkerFaceColor',cols(ci,:)); hold on;
        lg{end+1} = sprintf('%s, SNR=%d dB', configs{ci,1}, SNR_dbs(si)); %#ok<SAGROW>
    end
end
hold off; grid on; box on;
xlabel('velocity  [km/h]   (\nu_{max} = 3.7 Hz per km/h)');
ylabel('BER');
title(sprintf(['BER vs mobility, EVA grid, pilot lattice frozen (%d-QAM, %d trials)\n' ...
               'v=0 is ICI-free: any gap there is pure DIVERSITY, not Doppler'], ...
               M_order, nTrials));
legend(lg,'Location','southeast');  set(gca,'FontSize',11);
if exist('theme','file'), theme(gcf,"light"); end
exportgraphics(gcf, fullfile(resdir,'ber_vs_doppler.png'),'Resolution',150,'BackgroundColor','white');
fprintf('\nDone in %.0fs. Saved ber_vs_doppler.png\n', toc(tAll));
warning('on','MATLAB:nearlySingularMatrix');

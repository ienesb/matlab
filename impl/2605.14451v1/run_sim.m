% run_sim.m
%
% Reproduces Fig. 1 of:
%   F. Liu, Y. Xiong, Y.-F. Liu, J. Yang, C. Masouros, S. Jin,
%   "CP-OFDM Achieves Lower Ranging CRB Than Frequency-Spread
%   Waveforms in the Large-Sample Regime," arXiv:2605.14451 (2026).
%
% Plotted curves: SC, AFDM, OTFS, OFDM, Jensen Bound, for 16-QAM and
% 16-PSK constellations. All formulas follow Sections II--III and the
% parameter setup in Section V.C of the paper.

clear; clc; close all;

%% ===== System parameters (paper, Sec. V.C) ======================
N  = 128;        % number of subcarriers / transmitted symbols (eq. (1))
L  = 40;         % number of targets
N1 = 32;         % OTFS Kronecker factor 1 (eq. (8c)), N1*N2 = N
N2 = 4;          % OTFS Kronecker factor 2
c1 = 1/16;       % AFDM chirp rate 1 (eq. (8d)), 2*N*c1 in Z
c2 = 1/8;        % AFDM chirp rate 2
B  = 160e6;      % effective bandwidth (Hz)
c0 = 3e8;        % speed of light
range_factor = (c0/(2*B))^2;     % c^2/(4B^2): normalized delta-tau^2 -> m^2

SNR_dB  = -10:1:6;               % SNR = 1/sigma^2 (paper, Fig. 1 axis)
SNR_lin = 10.^(SNR_dB/10);
n_snr   = length(SNR_dB);

n_mc = 500;      % Monte Carlo trials over random data symbols

%% ===== Random target geometry (distinct delays, Lemma 1) ========
% L random distinct delays in [0, 127] (paper: 40 targets randomly placed
% in the normalized delay interval). A minimum separation of 1 sample is
% imposed to enforce Lemma 1 distinctness (tau_k ~= tau_l) numerically.
rng(20260514);
min_sep   = 1.0;                              % minimum delay separation
usable    = 127 - (L-1)*min_sep;              % "free" length to distribute
u         = sort(usable*rand(L,1));           % unconstrained order stats
tau       = u + min_sep*(0:L-1).';            % impose min separation
beta      = ones(L,1);                        % unit complex amplitudes

%% ===== Build Jacobian H = [H_tau, H_a, H_b] (eqs. (29)-(31)) ====
n_idx     = (0:N-1).';
% a(tau_l)_n = exp(-j*2*pi*(n-1)*tau_l/N)              (eq. (29))
a_mat     = exp(-1j*2*pi*(n_idx)*(tau.'/N));    % N x L
% a-dot(tau_l) = Diag{-j*2*pi*(n-1)/N} * a(tau_l)      (eq. (30))
adot_mat  = (-1j*2*pi*n_idx/N) .* a_mat;         % N x L
H_tau     = adot_mat .* beta.';                  % N x L  (eq. (31))
H_a       = a_mat;                               % N x L
H_b       = 1j*a_mat;                            % N x L
H         = [H_tau, H_a, H_b];                   % N x 3L

% Selection matrix T = T_tau^T * T_tau, T_tau = [I_L, 0_{L x 2L}] (eq. (37))
T = blkdiag(eye(L), zeros(2*L));                 % 3L x 3L

%% ===== Modulation bases U (eqs. (8a)-(8d)) and Q^H = F_N U =======
F_N  = fft(eye(N))  / sqrt(N);                   % F_N (eq. (9))
F_N1 = fft(eye(N1)) / sqrt(N1);                  % F_{N_1}
Lam_c1 = diag(exp(-1j*2*pi*c1*(n_idx.^2)));      % Lambda_{c1} (eq. (10))
Lam_c2 = diag(exp(-1j*2*pi*c2*(n_idx.^2)));      % Lambda_{c2}

U_sc   = eye(N);                                 % SC   : U = I_N        (8a)
U_ofdm = F_N';                                   % OFDM : U = F_N^H      (8b)
U_otfs = kron(F_N1', eye(N2));                   % OTFS : F_{N1}^H (x) I_{N2} (8c)
U_afdm = Lam_c1' * F_N' * Lam_c2';               % AFDM : Lc1^H F_N^H Lc2^H (8d)

% Q := U^H F_N^H  =>  Q^H = F_N U  (definition above eq. (29))
waveforms = {'SC','AFDM','OTFS','OFDM'};
QH = { F_N*U_sc, F_N*U_afdm, F_N*U_otfs, F_N*U_ofdm };
n_wf = numel(waveforms);

%% ===== Jensen lower bound (Proposition 1, eqs. (40)-(43)) =======
% E[|Q^H s|^2] = 1_N  (Q unitary, E[|s_n|^2]=1) =>
%   bar(J) = (2/sigma^2) Re{H^H H}
% CRB_Jensen(sigma^2) = (sigma^2/2) * Tr(T * inv(Re{H^H H}))
ReHH = real(H'*H);
crb_jensen_per_invSNR = trace(T/ReHH) / 2;       % times sigma^2 = 1/SNR
crb_jensen = (crb_jensen_per_invSNR ./ SNR_lin) * range_factor;

%% ===== Monte Carlo CRB per waveform, 16-QAM & 16-PSK ============
% f(U,s) = Tr{T J^{-1}}, J = (2/sigma^2) Re{H^H Diag(|Q^H s|^2) H} (eqs. 32,39)
constellations = {'16-QAM', '16-PSK'};
crb = zeros(n_snr, n_wf, 2);     % SNR x waveform x constellation

for ci = 1:2
    cname = constellations{ci};
    fprintf('Constellation %s:\n', cname);
    acc = zeros(1, n_wf);   cnt = zeros(1, n_wf);

    for tr = 1:n_mc
        s = gen_symbols(N, cname);
        for w = 1:n_wf
            p   = abs(QH{w}*s).^2;               % |Q^H s|^2          (eq. (32))
            J0  = real(H' * (p .* H));           % Re{H^H Diag(p) H}  (no 2/sig^2)
            val = trace(T / J0);                 % Tr{T J0^{-1}}
            if isfinite(val) && val > 0          % restrict to Omega_N (invertible FIM)
                acc(w) = acc(w) + val;  cnt(w) = cnt(w) + 1;
            end
        end
    end
    for w = 1:n_wf
        fprintf('   %-5s valid trials = %d / %d\n', waveforms{w}, cnt(w), n_mc);
        % CRB(m^2) = (sigma^2/2) * <Tr{T J0^{-1}}> * c^2/(4B^2)
        crb(:,w,ci) = (acc(w)/cnt(w)) ./ SNR_lin / 2 * range_factor;
    end
end

% Convenience aliases (per-constellation columns) for plotting / saving
crb_sc   = squeeze(crb(:,1,:));
crb_afdm = squeeze(crb(:,2,:));
crb_otfs = squeeze(crb(:,3,:));
crb_ofdm = squeeze(crb(:,4,:));

save('crb_new_results.mat', 'SNR_dB', 'crb', 'crb_jensen', ...
     'crb_sc', 'crb_afdm', 'crb_otfs', 'crb_ofdm', 'waveforms');

%% ===== Figures (a) 16-QAM and (b) 16-PSK =======================
make_fig(SNR_dB, crb_sc(:,1), crb_afdm(:,1), crb_otfs(:,1), crb_ofdm(:,1), ...
         crb_jensen, '-',  '(a) Ranging CRB under 16-QAM constellation', ...
         'fig_crb_16qam.png');
make_fig(SNR_dB, crb_sc(:,2), crb_afdm(:,2), crb_otfs(:,2), crb_ofdm(:,2), ...
         crb_jensen, 'kx', '(b) Ranging CRB under 16-PSK constellation', ...
         'fig_crb_16psk.png');


%% ===== Local functions ==========================================
function make_fig(snr, sc, afdm, otfs, ofdm, jensen, jstyle, ttl, fname)
    fig = figure('Color','w','Position',[100 100 700 500]);
    ax  = axes('Parent',fig,'Color','w','XColor','k','YColor','k', ...
               'GridColor',[0.3 0.3 0.3]); hold(ax,'on');
    plot(ax, snr, sc,   '-',  'LineWidth',1.8, 'Color',[0.00 0.45 0.74]);
    plot(ax, snr, afdm, '-.', 'LineWidth',1.8, 'Color',[0.49 0.18 0.56]);
    plot(ax, snr, otfs, '--', 'LineWidth',1.8, 'Color',[0.93 0.69 0.13]);
    plot(ax, snr, ofdm, '-',  'LineWidth',1.8, 'Color',[0.85 0.33 0.10]);
    if strcmp(jstyle,'kx')
        plot(ax, snr, jensen, 'kx', 'MarkerSize',9, 'LineWidth',1.5);
    else
        plot(ax, snr, jensen, '-',  'LineWidth',1.8, 'Color','k');
    end
    grid(ax,'on'); box(ax,'on');
    xlabel(ax,'SNR (dB)','Color','k');
    ylabel(ax,'CRB (m^2)','Color','k');
    legend(ax,{'SC','AFDM','OTFS','OFDM','Jensen Bound'}, ...
           'Location','northeast','TextColor','k','Color','w','EdgeColor','k');
    title(ax,ttl,'Color','k');
    xlim(ax,[snr(1) snr(end)]);
    exportgraphics(fig,fname,'Resolution',150,'BackgroundColor','white');
end

function s = gen_symbols(N, type)
    % Returns N x 1 i.i.d. symbol vector with E[|s|^2]=1,
    % E[s]=0, E[s^2]=0 (Assumptions 1-3 in the paper).
    switch upper(type)
        case '16-QAM'
            levels = [-3 -1 1 3];
            re = levels(randi(4, N, 1));   re = re(:);
            im = levels(randi(4, N, 1));   im = im(:);
            s  = (re + 1j*im) / sqrt(10);   % unit-power 16-QAM
        case '16-PSK'
            k = randi(16, N, 1) - 1;
            s = exp(1j*2*pi*k/16);
        otherwise
            error('Unknown constellation: %s', type);
    end
end

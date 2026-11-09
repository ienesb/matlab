% compute_mcrb_mc.m
%
% Monte-Carlo MISSPECIFIED lower bound (misspecified CRB at the pseudo-true
% parameter) for the ranging channel of
%   F. Liu, Y. Xiong, Y.-F. Liu, J. Yang, C. Masouros, S. Jin,
%   "CP-OFDM Achieves Lower Ranging CRB Than Frequency-Spread Waveforms
%    in the Large-Sample Regime," arXiv:2605.14451 (2026).
%
% This is the eva/compute_mcrb_mc.m machinery (pseudo-true projection +
% A,B,G,F,BB affine-in-SNR split) ported to the paper's channel model:
%   * delay-only ranging (NO Doppler, single OFDM symbol, M = 1);
%   * L targets with distinct normalized delays (Lemma 1), unit complex gains;
%   * the four waveforms SC / AFDM / OTFS / OFDM enter through the per-
%     subcarrier modulated symbol  g = Q^H s = F_N U s  (eqs. (8)-(9));
%   * parameter stacking per target eta_l = [Re(beta_l); Im(beta_l); tau_l],
%     tau in NORMALIZED samples (a_n(tau)=exp(-j2pi n tau/N), eq. (29)).
%
% Misspecification.  The eva BSC on QPSK is generalized to a constellation
% SYMBOL-ERROR rate p: each transmitted symbol is, with probability p, replaced
% by a uniformly random DIFFERENT constellation point, giving the receiver-side
% symbols shat.  The assumed (mismatched) mean is built from shat, the true mean
% from s.  This works for both 16-QAM and 16-PSK and reduces to the bit-flip
% BSC when the constellation is QPSK.
%
% For every realization (s, shat) and waveform:
%   1. g = Q^H s (true), ghat = Q^H shat (assumed).
%   2. eta_0 = argmin_eta || mu(g, eta_bar) - mu_tilde(ghat, eta) ||   (pseudo-true)
%      mu_n(g, eta) = g_n * sum_l beta_l exp(-j2pi n tau_l / N).
%   3. A(eta_0), B(eta_0), Jtilde(eta_0), b(eta_0) at sigma^2 = 1
%        B = Jtilde + b b^T,        A = -Jtilde + A_Delta,
%      with Delta = mu(g,eta_bar) - mu_tilde(ghat,eta_0).
%   4. SNR-affine lower bound (since A,B ~ 1/sigma^2):
%        LB(sigma^2) = sigma^2 * G + F + BB,        evaluated at sigma^2 = 1,
%        G  = A^-1 Jtilde A^-1     (CRB-like, scales with noise)
%        F  = A^-1 (b b^T) A^-1    (score-bias floor, sigma^2-independent)
%        BB = (eta_bar - eta_0)(eta_bar - eta_0)^T   (pseudo-true bias floor).
%
% At p = 0 (shat = s): eta_0 = eta_bar, A = -Jtilde, b = 0, F = BB = 0,
% G = Jtilde^-1, so LB = sigma^2 Jtilde^-1 and the ranging trace reproduces
% run_sim.m exactly (matched CRB of Fig. 1).
%
% Output: tr{T J^-1} over the tau block, converted to ranging MSE in m^2 via
% range_factor = (c0/(2B))^2, swept over SNR (paper axis) and the SER p.
clear; clc; close all;
warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:singularMatrix');
rng(20260514);   % same seed / target geometry as run_sim.m

%% ===== System parameters (paper, Sec. V.C) ======================
N  = 128;        % subcarriers / transmitted symbols (eq. (1))
L  = 40;         % number of targets
N1 = 32;         % OTFS Kronecker factor 1 (eq. (8c)), N1*N2 = N
N2 = 4;          % OTFS Kronecker factor 2
c1 = 1/16;       % AFDM chirp rate 1 (eq. (8d))
c2 = 1/8;        % AFDM chirp rate 2
B  = 160e6;      % effective bandwidth [Hz]
c0 = 3e8;        % speed of light
range_factor = (c0/(2*B))^2;     % normalized delta-tau^2 -> ranging m^2

SNR_dB  = -10:1:6;               % SNR = 1/sigma^2 (paper Fig. 1 axis)
SNR_lin = 10.^(SNR_dB/10);
n_snr   = numel(SNR_dB);

%% ===== Monte-Carlo / misspecification sweep =====================
nTrials = 40;                    % MC trials per (waveform, constellation, p)
p_list  = [0 0.05 0.10 0.20];    % symbol-error-rate sweep (0 = matched CRB)
np      = numel(p_list);
constellations = {'16-QAM','16-PSK'};
nConst  = numel(constellations);
p_show  = 0.10;                  % SER used for the comparison figures

%% ===== Random target geometry (distinct delays, Lemma 1) ========
min_sep = 1.0;
usable  = 127 - (L-1)*min_sep;
u       = sort(usable*rand(L,1));
tau     = u + min_sep*(0:L-1).';     % TRUE normalized delays [samples]
beta    = ones(L,1);                 % TRUE unit complex gains

nE      = 3*L;                       % parameters: [Re;Im;tau] per target
idx_tau = 3*(1:L);                   % tau positions in the stacking
eta_bar = packEta(beta, tau);        % ground-truth parameter (3L x 1)

%% ===== Modulation bases U and Q^H = F_N U (eqs. (8)-(9)) =========
n_idx  = (0:N-1).';
F_N    = fft(eye(N))  / sqrt(N);
F_N1   = fft(eye(N1)) / sqrt(N1);
Lam_c1 = diag(exp(-1j*2*pi*c1*(n_idx.^2)));
Lam_c2 = diag(exp(-1j*2*pi*c2*(n_idx.^2)));

U_sc   = eye(N);                                 % SC   (8a)
U_ofdm = F_N';                                   % OFDM (8b)
U_otfs = kron(F_N1', eye(N2));                   % OTFS (8c)
U_afdm = Lam_c1' * F_N' * Lam_c2';               % AFDM (8d)

waveforms = {'SC','AFDM','OTFS','OFDM'};
QH = { F_N*U_sc, F_N*U_afdm, F_N*U_otfs, F_N*U_ofdm };   % Q^H = F_N U
nWf = numel(waveforms);

% Constellation alphabets (unit power, E[s]=E[s^2]=0).
alphabets = cell(nConst,1);
for ci = 1:nConst, alphabets{ci} = constellation(constellations{ci}); end

%% ===== Storage (MC averages of the affine pieces) ===============
G_bar  = zeros(nE,nE,np,nWf,nConst);
F_bar  = zeros(nE,nE,np,nWf,nConst);
BB_bar = zeros(nE,nE,np,nWf,nConst);
trG    = zeros(np,nWf,nConst);   % <tr_tau{G}>  (CRB-like, /SNR slope)
trF    = zeros(np,nWf,nConst);   % <tr_tau{F}>  (score-bias floor)
trBB   = zeros(np,nWf,nConst);   % <tr_tau{BB}> (pseudo-true bias floor)

fprintf('MC misspecified CRB  LB(sigma2)=sigma2*G+F+BB  (ranging, paper channel)\n');
fprintf('N=%d L=%d, %d trials, p in [%s], SNR %d..%d dB\n', ...
        N, L, nTrials, num2str(p_list), SNR_dB(1), SNR_dB(end));

nWork = nConst*nWf*np;  doneWork = 0;  tAll = tic;

for ci = 1:nConst
    alphabet = alphabets{ci};
    for w = 1:nWf
        Qh = QH{w};
        for ip = 1:np
            pp = p_list(ip);
            G_t = zeros(nE,nE,nTrials);
            F_t = zeros(nE,nE,nTrials);
            BB_t = zeros(nE,nE,nTrials);

            parfor t = 1:nTrials
                % ---- data, symbol-error channel, modulated symbols ----
                [s, idx] = gen_symbols(alphabet, N);
                shat     = symbol_error(alphabet, idx, pp);
                g    = Qh * s;          % true per-subcarrier symbol  (eq. (9))
                ghat = Qh * shat;       % receiver-side (mismatched) symbol

                % ---- pseudo-true eta_0 (LS projection of the true mean) ----
                eta0 = pseudoTrue(g, ghat, eta_bar, n_idx, N, idx_tau);

                % ---- A, B, Jtilde, b and the affine pieces at sigma2 = 1 ----
                [A, ~, Jt, b] = build_AB(g, ghat, eta0, eta_bar, n_idx, N, 1);
                G  = stab_quad(A, Jt);             % A^-1 Jtilde A^-1
                F  = stab_quad(A, b*b.');          % A^-1 (b b^T) A^-1
                bias = eta_bar - eta0;
                BB = bias*bias.';                  % pseudo-true bias outer product

                G_t(:,:,t)  = (G+G.')/2;
                F_t(:,:,t)  = (F+F.')/2;
                BB_t(:,:,t) = BB;
            end

            G_bar(:,:,ip,w,ci)  = mean(G_t,3);
            F_bar(:,:,ip,w,ci)  = mean(F_t,3);
            BB_bar(:,:,ip,w,ci) = mean(BB_t,3);
            trG(ip,w,ci)  = sum(real(diag(G_bar(idx_tau,idx_tau,ip,w,ci))));
            trF(ip,w,ci)  = sum(real(diag(F_bar(idx_tau,idx_tau,ip,w,ci))));
            trBB(ip,w,ci) = sum(real(diag(BB_bar(idx_tau,idx_tau,ip,w,ci))));

            doneWork = doneWork + 1;
            fprintf('  %-6s %-6s p=%.2f : tr_tau{G,F,BB}=[%.2e %.2e %.2e]  (%4.0f%%, %s)\n', ...
                constellations{ci}, waveforms{w}, pp, ...
                trG(ip,w,ci), trF(ip,w,ci), trBB(ip,w,ci), ...
                100*doneWork/nWork, fmt_hms(toc(tAll)));
        end
    end
end

%% ===== Ranging MCRB (m^2): affine in 1/SNR per (p, waveform) =====
% LB_tau(sigma2) = sigma2*<tr_tau G> + <tr_tau F> + <tr_tau BB>,  sigma2 = 1/SNR.
% range MCRB [m^2] = LB_tau * range_factor.
mcrb = zeros(n_snr, np, nWf, nConst);
for ci = 1:nConst
    for w = 1:nWf
        for ip = 1:np
            mcrb(:,ip,w,ci) = ( (1./SNR_lin).' * trG(ip,w,ci) ...
                              + trF(ip,w,ci) + trBB(ip,w,ci) ) * range_factor;
        end
    end
end

save('compute_mcrb_mc_results.mat', 'SNR_dB','p_list','waveforms', ...
     'constellations','mcrb','trG','trF','trBB','G_bar','F_bar','BB_bar', ...
     'range_factor','N','L','B','nTrials','eta_bar','tau','idx_tau');
fprintf('\nSaved MC misspecified ranging CRB to compute_mcrb_mc_results.mat\n');

%% ===== Figures: matched (p=0) vs misspecified (p=p_show) =========
ip0 = find(p_list==0, 1);  ips = find(abs(p_list-p_show)<1e-9, 1);
if isempty(ip0),  ip0 = 1;        end
if isempty(ips),  ips = min(2,np); end
for ci = 1:nConst
    make_fig(SNR_dB, squeeze(mcrb(:,ip0,:,ci)), squeeze(mcrb(:,ips,:,ci)), ...
        waveforms, p_list(ips), ...
        sprintf('Ranging MCRB under %s (solid: SER=%.2f, dotted: matched)', ...
                constellations{ci}, p_list(ips)), ...
        sprintf('fig_mcrb_%s.png', lower(strrep(constellations{ci},'-',''))));
end
warning('on','MATLAB:nearlySingularMatrix');
warning('on','MATLAB:singularMatrix');


%% ================= Local functions ==============================
function alphabet = constellation(type)
% Unit-power constellation alphabet (column), E[s]=E[s^2]=0 (Assumptions 1-3).
    switch upper(type)
        case '16-QAM'
            lev = [-3 -1 1 3];
            [RE, IM] = ndgrid(lev, lev);
            alphabet = (RE(:) + 1j*IM(:)) / sqrt(10);
        case '16-PSK'
            k = (0:15).';
            alphabet = exp(1j*2*pi*k/16);
        otherwise
            error('Unknown constellation: %s', type);
    end
end

function [s, idx] = gen_symbols(alphabet, N)
% N x 1 i.i.d. symbols drawn uniformly from the alphabet; idx = chosen indices.
    K   = numel(alphabet);
    idx = randi(K, N, 1);
    s   = alphabet(idx);
end

function shat = symbol_error(alphabet, idx, p)
% Symbol-error channel: each symbol is, w.p. p, replaced by a uniformly random
% DIFFERENT alphabet point (generalizes the QPSK bit-flip BSC of eva/).
    K   = numel(alphabet);
    N   = numel(idx);
    hit = rand(N,1) < p;
    off = randi(K-1, N, 1);                 % 1..K-1 offset -> guarantees change
    idxh = idx;
    idxh(hit) = mod(idx(hit)-1 + off(hit), K) + 1;
    shat = alphabet(idxh);
end

function eta = packEta(beta, tau)
% Stack [Re(beta_l); Im(beta_l); tau_l] per target -> 3L x 1.
    eta = reshape([real(beta(:)).'; imag(beta(:)).'; tau(:).'], [], 1);
end

function [beta, tau] = unpackEta(eta)
    beta = eta(1:3:end) + 1j*eta(2:3:end);
    tau  = eta(3:3:end);
end

function mu = ranging_mean(g, eta, n_idx, N)
% mu_n = g_n * sum_l beta_l exp(-j2pi n tau_l / N)   (delay-only, M = 1).
    [beta, tau] = unpackEta(eta);
    Amat = exp(-1j*2*pi*(n_idx*(tau.'/N)));     % N x L,  a_n(tau_l)
    mu   = g .* (Amat * beta);
end

function [Gmat, mu] = ranging_grad(g, eta, n_idx, N)
% Gradient stack d mu / d eta (N x 3L), column order [Re,Im,tau] per target,
% and the mean mu. Used both as the lsqnonlin Jacobian and for A,B,Jtilde,b.
    [beta, tau] = unpackEta(eta);
    L    = numel(beta);
    Amat = exp(-1j*2*pi*(n_idx*(tau.'/N)));     % N x L
    cN   = -1j*2*pi*n_idx/N;                     % N x 1, d/dtau log a_n
    gA   = g .* Amat;                            % N x L, g_n a_n(tau_l)
    Gmat = zeros(N, 3*L);
    Gmat(:,1:3:end) =          gA;               % d/dRe(beta)
    Gmat(:,2:3:end) = 1j    .* gA;               % d/dIm(beta)
    Gmat(:,3:3:end) = (cN .* gA) .* beta.';      % d/dtau
    mu = gA * beta;
end

function eta0 = pseudoTrue(g, ghat, eta_bar, n_idx, N, idx_tau)
% eta0 = argmin_eta || mu(g,eta_bar) - mu_tilde(ghat,eta) ||  (Gauss-Newton).
% Analytic Jacobian = -[Re;Im] of the assumed-mean gradient. tau scaled by N
% so all variables are O(1) for the optimizer.
    mu_true = ranging_mean(g, eta_bar, n_idx, N);
    sc = ones(numel(eta_bar),1);  sc(idx_tau) = N;     % tau ~ N samples -> O(1)
    th0 = eta_bar ./ sc;
    function [r, J] = resfun(th)
        eta = th .* sc;
        [Gm, mu] = ranging_grad(ghat, eta, n_idx, N);
        d = mu_true - mu;
        r = [real(d); imag(d)];
        if nargout > 1
            Jc = -Gm .* sc.';                          % chain rule for scaling
            J  = [real(Jc); imag(Jc)];
        end
    end
    opts = optimoptions('lsqnonlin','Display','off','SpecifyObjectiveGradient',true, ...
        'FunctionTolerance',1e-14,'StepTolerance',1e-12,'OptimalityTolerance',1e-12, ...
        'MaxFunctionEvaluations',400*numel(eta_bar),'MaxIterations',200);
    th = lsqnonlin(@resfun, th0, [], [], opts);
    eta0 = th .* sc;
end

function [A, B, Jt, b] = build_AB(g, ghat, eta, eta_bar, n_idx, N, sigma2)
% Realization MCRB components at the evaluation point eta (delay-only model).
%   Jt = (2/s2) Re{G^H G},  b = (2/s2) Re{G^H Delta},  B = Jt + b b^T,
%   A  = -Jt + A_Delta,  Delta = mu(g,eta_bar) - mu_tilde(ghat,eta).
% G and the Hessian use the receiver-side symbols ghat at eta; the true mean
% uses g and eta_bar.  Only [Re,tau], [Im,tau], [tau,tau] Hessians are nonzero
% (the mean is linear in beta), and cross-target Hessians vanish.
    [Gmat, mu_as] = ranging_grad(ghat, eta, n_idx, N);
    mu_true = ranging_mean(g, eta_bar, n_idx, N);
    Delta   = mu_true - mu_as;

    Jt = (2/sigma2) * real(Gmat' * Gmat);
    b  = (2/sigma2) * real(Gmat' * Delta);
    B  = Jt + b*b.';

    [beta, tau] = unpackEta(eta);
    L    = numel(beta);
    Amat = exp(-1j*2*pi*(n_idx*(tau.'/N)));
    cN   = -1j*2*pi*n_idx/N;
    gA   = ghat .* Amat;                          % N x L
    Adelta = zeros(3*L, 3*L);
    cd = @(h) (2/sigma2) * real(sum(conj(Delta) .* h));   % contraction with Delta
    for l = 1:L
        iR = 3*(l-1)+1; iI = iR+1; iT = iR+2;
        gAl = gA(:,l);
        hRt = cN .* gAl;                          % d^2 mu / dRe dtau
        hIt = 1j * (cN .* gAl);                   % d^2 mu / dIm dtau
        htt = beta(l) * (cN.^2 .* gAl);           % d^2 mu / dtau^2
        Adelta(iR,iT) = cd(hRt); Adelta(iT,iR) = Adelta(iR,iT);
        Adelta(iI,iT) = cd(hIt); Adelta(iT,iI) = Adelta(iI,iT);
        Adelta(iT,iT) = cd(htt);
    end
    A = -Jt + Adelta;
end

function Y = stab_quad(A, C)
% Y = A^-1 C A^-1 with a Jacobi (diagonal, magnitude) preconditioner (real out).
    d = 1 ./ sqrt(abs(diag(A)));
    d(~isfinite(d)) = 1;
    D = diag(d);
    As = D*A*D;  Cs = D*C*D;
    Y  = real(D*(As \ Cs / As)*D);
end

function s = fmt_hms(sec)
    sec = max(0, round(sec));
    s = sprintf('%02d:%02d:%02d', floor(sec/3600), floor(mod(sec,3600)/60), mod(sec,60));
end

function make_fig(snr, crb0, crbp, waveforms, pser, ttl, fname)
% Solid: misspecified (SER=pser).  Dotted: matched (p=0).  One color/waveform.
    cols = [0.00 0.45 0.74; 0.49 0.18 0.56; 0.93 0.69 0.13; 0.85 0.33 0.10];
    fig = figure('Color','w','Position',[100 100 720 520]);
    ax  = axes('Parent',fig,'Color','w','XColor','k','YColor','k'); hold(ax,'on');
    h = gobjects(numel(waveforms),1);
    for w = 1:numel(waveforms)
        plot(ax, snr, crb0(:,w), ':',  'LineWidth',1.4, 'Color',cols(w,:));
        h(w) = plot(ax, snr, crbp(:,w), '-', 'LineWidth',1.8, 'Color',cols(w,:));
    end
    set(ax,'YScale','log'); grid(ax,'on'); box(ax,'on');
    xlabel(ax,'SNR (dB)','Color','k'); ylabel(ax,'Ranging MCRB (m^2)','Color','k');
    legend(ax, h, waveforms, 'Location','northeast', ...
           'TextColor','k','Color','w','EdgeColor','k');
    title(ax, ttl, 'Color','k'); xlim(ax,[snr(1) snr(end)]);
    exportgraphics(fig, fname, 'Resolution',150, 'BackgroundColor','white');
end

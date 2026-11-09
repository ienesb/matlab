% verify_my_derivation.m
%
% Computes the multi-target FIM and ranging CRB for the EVA/paper channel
% using the user's derivation in ofdm_fim.tex (Section 3.2.1, "No ICI
% Channel" — multi-target, single OFDM symbol, generic |X|^2 form), and
% compares it with the paper's CRB stored in crb_new_results.mat.
%
% Doppler is eliminated by setting M = 1: the inner I-kernel reduces to
% delta[n'-n], I'(.) = 0, and the parameter vector per target becomes
% [Re(alpha_p), Im(alpha_p), tau_p]. The user's QPSK simplification is
% replaced by the constellation-dependent |s_n|^2 weighting, exactly as
% the paper does for 16-QAM / 16-PSK.

clear; clc; close all;
S = load('crb_new_results.mat');     % saved paper-style curves
SNR_dB  = S.SNR_dB;
SNR_lin = 10.^(SNR_dB/10);
n_snr   = numel(SNR_dB);
n_mc    = 500;

%% ===== Identical channel setup to run_sim.m =====================
N  = 128;        L  = 40;
B  = 160e6;      c0 = 3e8;

rng(20260514);
min_sep   = 1.0;
usable    = 127 - (L-1)*min_sep;
u         = sort(usable*rand(L,1));
tau_norm  = u + min_sep*(0:L-1).';      % paper's tau (normalized samples)
beta      = ones(L,1);                  % complex amplitudes

% Map to physical OFDM grid used in the user's derivation:
Df       = B / N;                       % subcarrier spacing
tau_phys = tau_norm / B;                % physical delay (seconds)
n_idx    = (0:N-1).';
DT       = tau_phys.' - tau_phys;       % (p,q) -> tau_q - tau_p, L x L
tauPos   = 3:3:3*L;                     % tau positions in user's stacking

%% ===== Per-realization sanity check ==============================
% Build both FIMs for ONE shared realization, then verify the user's
% FIM equals the paper's FIM after (i) parameter reordering and (ii)
% scaling tau rows/columns by B. This is the strongest algebraic test.

% Paper-side Jacobian H = [H_tau, H_a, H_b]
a_mat    = exp(-1j*2*pi*(n_idx)*(tau_norm.'/N));
adot_mat = (-1j*2*pi*n_idx/N) .* a_mat;
H_paper  = [adot_mat .* beta.', a_mat, 1j*a_mat];

s_test = gen_symbols(N, '16-QAM');
sn2    = abs(s_test).^2;
J0_paper  = real(H_paper' * (sn2 .* H_paper));        % paper layout: [tau; R; I]
Jcore_usr = build_user_fim(s_test, n_idx, Df, DT, beta, L);

% Permutation from paper layout to user layout
Pperm = zeros(3*L);
for p = 1:L
    Pperm((p-1)*3 + 1, L  + p) = 1;     % R_p
    Pperm((p-1)*3 + 2, 2*L+ p) = 1;     % I_p
    Pperm((p-1)*3 + 3,       p) = 1;    % tau_p
end
% Unit rescaling for tau (chain rule, tau_user = tau_paper/B)
Dscale = ones(3*L,1);  Dscale(tauPos) = B;
J0_paper_mapped = diag(Dscale)*(Pperm*J0_paper*Pperm')*diag(Dscale);

err_rel = norm(Jcore_usr - J0_paper_mapped, 'fro') / norm(Jcore_usr, 'fro');
fprintf('Per-realization FIM relative error (after layout+unit map) = %.3e\n', err_rel);

%% ===== Monte-Carlo ranging CRB from the user's FIM ==============
constellations = {'16-QAM','16-PSK'};
crb_user = zeros(n_snr, 2);

% RNG state is exactly where run_sim.m was when its MC loops started
% (same seed, same tau-generation sequence, no other rand calls since).
for ci = 1:2
    cname = constellations{ci};
    fprintf('Constellation %s ... ', cname);
    acc = 0; cnt = 0;
    for tr = 1:n_mc
        s = gen_symbols(N, cname);
        Jc = build_user_fim(s, n_idx, Df, DT, beta, L);
        Ji = Jc \ eye(3*L);
        tr_T = sum(Ji(sub2ind(size(Ji), tauPos, tauPos)));
        if isfinite(tr_T) && tr_T > 0
            acc = acc + tr_T;  cnt = cnt + 1;
        end
    end
    fprintf('valid=%d/%d\n', cnt, n_mc);
    % Range CRB (m^2): (c/2)^2 * (sigma^2/2) * <sum [J^-1]_{tau_p,tau_p}>
    crb_user(:,ci) = (acc/cnt) ./ SNR_lin / 2 * (c0/2)^2;
end

save('crb_user_verify.mat', 'SNR_dB', 'crb_user', 'err_rel');

%% ===== Plots: user curve overlaid on paper curve =================
make_plot(SNR_dB, S.crb_ofdm(:,1), crb_user(:,1), ...
          '16-QAM: OFDM ranging CRB', 'fig_verify_16qam.png');
make_plot(SNR_dB, S.crb_ofdm(:,2), crb_user(:,2), ...
          '16-PSK: OFDM ranging CRB', 'fig_verify_16psk.png');


%% ===== Helpers ====================================================
function make_plot(snr, paper, user, ttl, fname)
    fig = figure('Color','w','Position',[100 100 700 500]);
    ax  = axes('Parent',fig,'Color','w','XColor','k','YColor','k', ...
               'GridColor',[0.3 0.3 0.3]); hold(ax,'on');
    plot(ax, snr, paper, '-',  'LineWidth', 1.8, 'Color', [0.85 0.33 0.10]);
    plot(ax, snr, user,  'k+', 'MarkerSize', 9, 'LineWidth', 1.4);
    grid(ax,'on'); box(ax,'on');
    xlabel(ax,'SNR (dB)','Color','k');
    ylabel(ax,'CRB (m^2)','Color','k');
    legend(ax,{'Paper formulation (run\_sim.m)', ...
               'User derivation (ofdm\_fim.tex)'}, ...
           'Location','northeast','TextColor','k','Color','w','EdgeColor','k');
    title(ax,ttl,'Color','k');
    xlim(ax,[snr(1) snr(end)]);
    exportgraphics(fig,fname,'Resolution',150,'BackgroundColor','white');
end

function J = build_user_fim(s, n_idx, Df, DT, beta, L)
    % Section 3.2.1 of ofdm_fim.tex with M=1, nu=0, generic |s_n|^2.
    % Returns J without the leading 2/sigma^2 (caller scales by sigma^2/2).
    N   = numel(s);
    sn2 = abs(s).^2;
    % Lambda^{(a)}_{p,q} = sum_n (n*Df)^a |s_n|^2 exp(-j2*pi*n*Df*DT(p,q))
    Ph  = exp(-1j*2*pi*Df*(n_idx * DT(:).'));    % N x L^2
    w0  = sn2;
    w1  = (n_idx*Df).* sn2;
    w2  = (n_idx*Df).^2.* sn2;
    L00 = reshape(w0.'*Ph, L, L);
    L10 = reshape(w1.'*Ph, L, L);
    L20 = reshape(w2.'*Ph, L, L);

    J = zeros(3*L);
    for p = 1:L
        ip = (p-1)*3 + (1:3);   ap = beta(p);
        for q = 1:L
            iq = (q-1)*3 + (1:3);   aq = beta(q);
            l00 = L00(p,q);  l10 = L10(p,q);  l20 = L20(p,q);

            J(ip(1), iq(1)) =  real(l00);                  % R_p R_q
            J(ip(1), iq(2)) = -imag(l00);                  % R_p I_q
            J(ip(1), iq(3)) =  2*pi*imag(aq*l10);          % R_p tau_q
            J(ip(2), iq(1)) =  imag(l00);                  % I_p R_q
            J(ip(2), iq(2)) =  real(l00);                  % I_p I_q
            J(ip(2), iq(3)) = -2*pi*real(aq*l10);          % I_p tau_q
            J(ip(3), iq(1)) = -2*pi*imag(conj(ap)*l10);    % tau_p R_q
            J(ip(3), iq(2)) = -2*pi*real(conj(ap)*l10);    % tau_p I_q
            J(ip(3), iq(3)) =  4*pi^2*real(conj(ap)*aq*l20);% tau_p tau_q
        end
    end
end

function s = gen_symbols(N, type)
    switch upper(type)
        case '16-QAM'
            levels = [-3 -1 1 3];
            re = levels(randi(4, N, 1));  re = re(:);
            im = levels(randi(4, N, 1));  im = im(:);
            s = (re + 1j*im) / sqrt(10);
        case '16-PSK'
            k = randi(16, N, 1) - 1;
            s = exp(1j*2*pi*k/16);
    end
end

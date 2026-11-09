function [A, B, Jt, b] = getAB_OFDM_core(X, Xhat, deltaf, Ts, sigma2, alpha, tau, nu, use_ici, alpha_bar, tau_bar, nu_bar)
% Core realization MCRB component matrices A and B for OFDM (ofdm_mcrb.tex).
% Handles single/multi target and no-ICI/ICI through (numel(alpha), use_ici).
% Parameter stacking: eta = [eta_1; ...; eta_P], eta_p = [aR_p, aI_p, tau_p, nu_p].
%
%   B = Jtilde + b b^T,         A = -Jtilde + A_Delta,
%   [Jtilde]_kl = (2/s2) Re{ gt_k^H gt_l },     (FIM evaluated at Xhat, eta)
%   [b]_k       = (2/s2) Re{ Delta^H gt_k },
%   [A_Delta]_kl= (2/s2) Re{ Delta^H ht_kl },
% where gt_k = d mu_tilde/d eta_k and ht_kl = d^2 mu_tilde/d eta_k d eta_l are
% built from the receiver-side symbols Xhat at the EVALUATION parameter
% eta = (alpha, tau, nu), and the mismatch is the full signal difference
%   Delta = mu(X, eta_bar) - mu_tilde(Xhat, eta),
% with the TRUE mean mu built from the transmitted symbols X and the true
% parameter eta_bar = (alpha_bar, tau_bar, nu_bar).
%
% By default eta_bar = eta (the matrices are evaluated at the ground-truth
% parameter, and Delta reduces to sum_p alpha_p * Phi_{E,p}, E = X - Xhat).  To
% evaluate the MCRB at the PSEUDO-TRUE parameter (pseudo_true.tex), pass eta as
% the pseudo-true eta_0 and eta_bar as the true parameter -- see
% getAB_OFDM_pseudoTrue.
%
% Uses I' := (1/j2pi) dI/dnu and I'' := (1/j2pi) dI'/dnu (ofdm_mcrb/ofdm_fim
% convention); get_ici_coeffs_derivative{,2} return the full dI/dnu, d^2I/dnu^2.
%
% Inputs:
%   X, Xhat : N x M transmitted and receiver-side symbol grids
%   deltaf  : subcarrier spacing [Hz]
%   Ts      : OFDM symbol duration incl. CP [s]
%   sigma2  : noise variance
%   alpha   : P x 1 complex path gains   (EVALUATION point eta)
%   tau     : P x 1 path delays   [s]
%   nu      : P x 1 path Dopplers [Hz]
%   use_ici : true for the ICI kernel, false for the no-ICI (I = delta) model
%   alpha_bar, tau_bar, nu_bar : (optional) TRUE parameter eta_bar for the true
%             mean in Delta; default = (alpha, tau, nu).
%
% Outputs:
%   A, B    : 4P x 4P MCRB component matrices.

    [N, M] = size(X);
    P  = numel(alpha);
    T  = 1/deltaf;

    if nargin < 12 || isempty(alpha_bar)
        alpha_bar = alpha;  tau_bar = tau;  nu_bar = nu;
    end

    ns = (0:N-1).';     % subcarrier index (n' inside the ICI sum)
    ms = (0:M-1);       % OFDM-symbol index

    if use_ici
        [NP, NN] = ndgrid(0:N-1, 0:N-1);
        diff_n = NP - NN;
    end

    % True noise-free mean mu(X, eta_bar): the projection target of Delta.
    mu_true = ofdm_mean(X, deltaf, Ts, alpha_bar, tau_bar, nu_bar, use_ici);

    G      = zeros(N*M, 4*P);    % gradient stack gt_k (at eta)
    mu_as  = zeros(N*M, 1);      % assumed mean mu_tilde(Xhat, eta)

    % Per-target Hessian shorthand fields, stored for the A_Delta pass.
    Psi_c = cell(P,1);  Om_c = cell(P,1);  Xi_c = cell(P,1);
    Th_c  = cell(P,1);  Om1_c = cell(P,1);

    for p = 1:P
        ap = alpha(p);
        phase = exp(-1j*2*pi*deltaf*ns*tau(p)) * exp(1j*2*pi*Ts*ms*nu(p));  % N x M
        Ahat  = Xhat .* phase;

        if use_ici
            Imat = get_ici_coeffs           (diff_n, nu(p), T);
            Ip   = get_ici_coeffs_derivative (diff_n, nu(p), T) / (1j*2*pi);    % I'
            Ipp  = get_ici_coeffs_derivative2(diff_n, nu(p), T) / (1j*2*pi)^2;  % I''

            Phi  = Imat.' * Ahat;
            Psi  = Imat.' * ((ns*deltaf)    .* Ahat);
            Xi   = Imat.' * ((ns*deltaf).^2 .* Ahat);
            D1   = Ip.'   * Ahat;                       % sum Xhat phi I'
            D1n  = Ip.'   * ((ns*deltaf)    .* Ahat);   % sum Xhat phi n'df I'
            D2   = Ipp.'  * Ahat;                       % sum Xhat phi I''

            Omega  = (ms*Ts)    .* Phi + D1;
            Omega1 = (ms*Ts)    .* Psi + D1n;
            Theta  = (ms*Ts).^2 .* Phi + 2*(ms*Ts).*D1 + D2;
        else
            Phi    = Ahat;                              % I = delta
            Psi    = (ns*deltaf)    .* Phi;
            Xi     = (ns*deltaf).^2 .* Phi;
            Omega  = (ms*Ts)        .* Phi;
            Omega1 = (ns*deltaf) .* (ms*Ts) .* Phi;
            Theta  = (ms*Ts).^2     .* Phi;
        end

        col = 4*(p-1);
        G(:, col+1) =                 Phi(:);
        G(:, col+2) = 1j           *  Phi(:);
        G(:, col+3) = -1j*2*pi*ap  *  Psi(:);
        G(:, col+4) =  1j*2*pi*ap  *  Omega(:);

        mu_as = mu_as + ap * Phi(:);

        Psi_c{p} = Psi(:);  Om_c{p} = Omega(:);  Xi_c{p} = Xi(:);
        Th_c{p}  = Theta(:); Om1_c{p} = Omega1(:);
    end

    Delta = mu_true(:) - mu_as;                         % mu(X,eta_bar) - mu_tilde(Xhat,eta)

    % Jtilde (FIM at Xhat) and bias vector.
    Jt = (2/sigma2) * real(G' * G);
    b  = (2/sigma2) * real(G' * Delta);     % 4P x 1
    B  = Jt + b*b.';

    % A_Delta: block-diagonal (cross-target Hessians vanish).
    Adelta = zeros(4*P, 4*P);
    contract = @(h) (2/sigma2) * real(sum(conj(Delta) .* h));
    for p = 1:P
        ap = alpha(p);
        i1 = 4*(p-1)+1; i2 = i1+1; i3 = i1+2; i4 = i1+3;

        hARt = -1j*2*pi      * Psi_c{p};
        hARn =  1j*2*pi      * Om_c{p};
        hAIt =  2*pi         * Psi_c{p};
        hAIn = -2*pi         * Om_c{p};
        htt  = -4*pi^2 * ap  * Xi_c{p};
        hnn  = -4*pi^2 * ap  * Th_c{p};
        htn  =  4*pi^2 * ap  * Om1_c{p};

        Adelta(i1,i3) = contract(hARt);  Adelta(i3,i1) = Adelta(i1,i3);
        Adelta(i1,i4) = contract(hARn);  Adelta(i4,i1) = Adelta(i1,i4);
        Adelta(i2,i3) = contract(hAIt);  Adelta(i3,i2) = Adelta(i2,i3);
        Adelta(i2,i4) = contract(hAIn);  Adelta(i4,i2) = Adelta(i2,i4);
        Adelta(i3,i3) = contract(htt);
        Adelta(i4,i4) = contract(hnn);
        Adelta(i3,i4) = contract(htn);   Adelta(i4,i3) = Adelta(i3,i4);
    end

    A = -Jt + Adelta;
end

% ----- noise-free OFDM mean mu(S, eta) (ofdm_mcrb.tex eqs (3)-(4)) -------
function mu = ofdm_mean(S, deltaf, Ts, alpha, tau, nu, use_ici)
    [N, M] = size(S);  P = numel(alpha);
    ns = (0:N-1).';  ms = (0:M-1);  mu = zeros(N, M);
    if use_ici
        T = 1/deltaf;  [NP, NN] = ndgrid(0:N-1, 0:N-1);  diff_n = NP - NN;
    end
    for p = 1:P
        phase = exp(-1j*2*pi*deltaf*ns*tau(p)) * exp(1j*2*pi*Ts*ms*nu(p));
        A = S .* phase;
        if use_ici
            mu = mu + alpha(p) * (get_ici_coeffs(diff_n, nu(p), T).' * A);
        else
            mu = mu + alpha(p) * A;
        end
    end
end

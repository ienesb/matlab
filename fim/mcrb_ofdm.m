function out = mcrb_ofdm(X, X_hat, deltaf, Ts, sigma2, alphas, taus, nus)
% MCRB component matrices for OFDM with ICI (ofdm_mcrb.tex Sections 1-2).
% Supports both single (P=1) and multi-target (P>1) cases.
% Parameter stacking: eta = [eta_1; ...; eta_P], eta_p = [aR_p, aI_p, tau_p, nu_p].
%
% Inputs:
%   X      : N x M true symbol grid
%   X_hat  : N x M assumed (BSC-corrupted) symbol grid
%   ...    : scalar/vector parameters
%
% Output struct fields (all 4P x 4P unless noted):
%   J_tilde : tilde-FIM at X_hat
%   b       : 4P x 1 bias vector
%   A_Delta : Hessian-mean correction
%   A       = -J_tilde + A_Delta
%   B       = J_tilde + b*b'
%   Delta   : N x M mismatch signal mu(X) - mu(X_hat)

    [N, M] = size(X);
    P = numel(alphas);
    T = 1/deltaf;

    ns = (0:N-1).';
    ms = (0:M-1);
    ds = (-(N-1):(N-1)).';
    E  = X - X_hat;

    G       = zeros(N*M, 4*P);
    Delta   = zeros(N, M);
    A_Delta = zeros(4*P, 4*P);

    Psi_c    = cell(P,1);
    Omega_c  = cell(P,1);
    Xi_c     = cell(P,1);
    Omega1_c = cell(P,1);
    Theta_c  = cell(P,1);

    for p = 1:P
        a_p   = alphas(p);
        tau_p = taus(p);
        nu_p  = nus(p);

        phase  = exp(-1j*2*pi*deltaf*ns*tau_p) .* exp(1j*2*pi*Ts*ms*nu_p);
        A_arr  = X_hat .* phase;
        E_arr  = E     .* phase;

        I_full   = get_ici_coeffs           (ds, nu_p, T);
        dI_full  = get_ici_coeffs_derivative(ds, nu_p, T);
        ddI_full = get_ici_coeffs_derivative2(ds, nu_p, T);

        I_mat   = toeplitz(I_full  (N:-1:1), I_full  (N:end));
        dI_mat  = toeplitz(dI_full (N:-1:1), dI_full (N:end));
        ddI_mat = toeplitz(ddI_full(N:-1:1), ddI_full(N:end));

        % Six shorthand quantities at X_hat
        Phi    = I_mat * A_arr;
        Psi    = I_mat * (A_arr .* (ns*deltaf));
        Xi     = I_mat * (A_arr .* (ns*deltaf).^2);
        Omega  = (ms*Ts) .* Phi + dI_mat * A_arr;
        Omega1 = (ms*Ts) .* Psi + dI_mat * (A_arr .* (ns*deltaf));
        Theta  = (ms*Ts).^2 .* Phi + 2*(ms*Ts) .* (dI_mat * A_arr) + ddI_mat * A_arr;

        % Delta contribution from E
        Phi_E = I_mat * E_arr;
        Delta = Delta + a_p * Phi_E;

        idx = 4*(p-1);
        G(:, idx+1) =  Phi(:);
        G(:, idx+2) =  1j*Phi(:);
        G(:, idx+3) = -1j*2*pi*a_p*Psi(:);
        G(:, idx+4) =  1j*2*pi*a_p*Omega(:);

        Psi_c{p}    = Psi(:);
        Omega_c{p}  = Omega(:);
        Xi_c{p}     = Xi(:);
        Omega1_c{p} = Omega1(:);
        Theta_c{p}  = Theta(:);
    end

    J_tilde = (2/sigma2) * real(G' * G);

    Delta_vec = Delta(:);
    b = (2/sigma2) * real(G' * Delta_vec);   % 4P x 1

    % A_Delta: nonzero only inside each per-target 4x4 block.
    for p = 1:P
        idx = 4*(p-1);
        a_p = alphas(p);
        Psi    = Psi_c{p};
        Omega  = Omega_c{p};
        Xi     = Xi_c{p};
        Omega1 = Omega1_c{p};
        Theta  = Theta_c{p};

        v_aRt = -1j*2*pi*Psi;
        v_aRn =  1j*2*pi*Omega;
        v_aIt =  2*pi*Psi;
        v_aIn = -2*pi*Omega;
        v_tt  = -4*pi^2 * a_p * Xi;
        v_nn  = -4*pi^2 * a_p * Theta;
        v_tn  =  4*pi^2 * a_p * Omega1;

        A_Delta(idx+1, idx+3) = (2/sigma2) * real(Delta_vec' * v_aRt);
        A_Delta(idx+1, idx+4) = (2/sigma2) * real(Delta_vec' * v_aRn);
        A_Delta(idx+2, idx+3) = (2/sigma2) * real(Delta_vec' * v_aIt);
        A_Delta(idx+2, idx+4) = (2/sigma2) * real(Delta_vec' * v_aIn);
        A_Delta(idx+3, idx+3) = (2/sigma2) * real(Delta_vec' * v_tt);
        A_Delta(idx+4, idx+4) = (2/sigma2) * real(Delta_vec' * v_nn);
        A_Delta(idx+3, idx+4) = (2/sigma2) * real(Delta_vec' * v_tn);

        % Symmetrize
        A_Delta(idx+3, idx+1) = A_Delta(idx+1, idx+3);
        A_Delta(idx+4, idx+1) = A_Delta(idx+1, idx+4);
        A_Delta(idx+3, idx+2) = A_Delta(idx+2, idx+3);
        A_Delta(idx+4, idx+2) = A_Delta(idx+2, idx+4);
        A_Delta(idx+4, idx+3) = A_Delta(idx+3, idx+4);
    end

    A = -J_tilde + A_Delta;
    B =  J_tilde + b * b.';

    out.J_tilde = J_tilde;
    out.b       = b;
    out.A_Delta = A_Delta;
    out.A       = A;
    out.B       = B;
    out.Delta   = Delta;
end

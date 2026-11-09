function out = mcrb_otfs(x, x_hat, deltaf, Ts, sigma2, alphas, taus, nus)
% MCRB component matrices for OTFS, no OFDM-side ICI form (otfs_mcrb.tex
% Sections 1-2, subsection "No OFDM-side ICI"). Supports both single (P=1)
% and multi-target (P>1) cases.
% Parameter stacking: eta = [eta_1; ...; eta_P], eta_p = [aR_p, aI_p, tau_p, nu_p].
%
% Inputs / Output: same shape as mcrb_ofdm.

    [N, M] = size(x);
    P = numel(alphas);

    ns   = (0:N-1).';
    ms   = (0:M-1).';
    ds_A = (-(N-1):(N-1)).';
    ds_B = (-(M-1):(M-1)).';
    e    = x - x_hat;

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

        % Delay kernels A_p, A'_p, A''_p
        E_A   = exp(1j*2*pi*ns*(ds_A.'/N - deltaf*tau_p));
        A_d   = sum(E_A,        1).';
        Ap_d  = sum(ns   .* E_A,1).';
        App_d = sum(ns.^2.* E_A,1).';

        A_mat   = toeplitz(A_d  (N:end), A_d  (N:-1:1));
        Ap_mat  = toeplitz(Ap_d (N:end), Ap_d (N:-1:1));
        App_mat = toeplitz(App_d(N:end), App_d(N:-1:1));

        % Doppler kernels B_p, B'_p, B''_p
        E_B   = exp(1j*2*pi*ms*(ds_B.'/M + Ts*nu_p));
        B_d   = sum(E_B,        1).';
        Bp_d  = sum(ms   .* E_B,1).';
        Bpp_d = sum(ms.^2.* E_B,1).';

        B_mat   = toeplitz(B_d  (M:end), B_d  (M:-1:1));
        Bp_mat  = toeplitz(Bp_d (M:end), Bp_d (M:-1:1));
        Bpp_mat = toeplitz(Bpp_d(M:end), Bpp_d(M:-1:1));

        scale = 1/(N*M);
        Phi    = scale * (A_mat   * x_hat * B_mat );
        Psi    = scale * (Ap_mat  * x_hat * B_mat );
        Omega  = scale * (A_mat   * x_hat * Bp_mat);
        Xi     = scale * (App_mat * x_hat * B_mat );
        Omega1 = scale * (Ap_mat  * x_hat * Bp_mat);
        Theta  = scale * (A_mat   * x_hat * Bpp_mat);

        % Delta contribution from e = x - x_hat
        Phi_e = scale * (A_mat * e * B_mat);
        Delta = Delta + a_p * Phi_e;

        idx = 4*(p-1);
        G(:, idx+1) =  Phi(:);
        G(:, idx+2) =  1j*Phi(:);
        G(:, idx+3) = -1j*2*pi*deltaf*a_p*Psi(:);
        G(:, idx+4) =  1j*2*pi*Ts    *a_p*Omega(:);

        Psi_c{p}    = Psi(:);
        Omega_c{p}  = Omega(:);
        Xi_c{p}     = Xi(:);
        Omega1_c{p} = Omega1(:);
        Theta_c{p}  = Theta(:);
    end

    J_tilde = (2/sigma2) * real(G' * G);

    Delta_vec = Delta(:);
    b = (2/sigma2) * real(G' * Delta_vec);

    for p = 1:P
        idx = 4*(p-1);
        a_p = alphas(p);
        Psi    = Psi_c{p};
        Omega  = Omega_c{p};
        Xi     = Xi_c{p};
        Omega1 = Omega1_c{p};
        Theta  = Theta_c{p};

        v_aRt = -1j*2*pi*deltaf * Psi;
        v_aRn =  1j*2*pi*Ts     * Omega;
        v_aIt =  2*pi*deltaf    * Psi;
        v_aIn = -2*pi*Ts        * Omega;
        v_tt  = -4*pi^2*(deltaf)^2 * a_p * Xi;
        v_nn  = -4*pi^2*(Ts)^2     * a_p * Theta;
        v_tn  =  4*pi^2*deltaf*Ts  * a_p * Omega1;

        A_Delta(idx+1, idx+3) = (2/sigma2) * real(Delta_vec' * v_aRt);
        A_Delta(idx+1, idx+4) = (2/sigma2) * real(Delta_vec' * v_aRn);
        A_Delta(idx+2, idx+3) = (2/sigma2) * real(Delta_vec' * v_aIt);
        A_Delta(idx+2, idx+4) = (2/sigma2) * real(Delta_vec' * v_aIn);
        A_Delta(idx+3, idx+3) = (2/sigma2) * real(Delta_vec' * v_tt);
        A_Delta(idx+4, idx+4) = (2/sigma2) * real(Delta_vec' * v_nn);
        A_Delta(idx+3, idx+4) = (2/sigma2) * real(Delta_vec' * v_tn);

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

function J = getFIM_OFDM_multi_noICI_expected(N, M, deltaf, Ts, sigma2, alpha, tau, nu)
% Multi-target EXPECTED OFDM Fisher information matrix, no-ICI channel.
% Closed form for i.i.d. unit-energy QPSK symbols, i.e. the data-averaged FIM
% of Section 3.2.1 / eqs. (234)-(255) of ofdm_fim.tex.  Parameter stacking:
%   eta = [eta_1; ...; eta_P],  eta_p = [alphaR_p, alphaI_p, tau_p, nu_p].
%
% Off-diagonal blocks depend only on the inter-target separations
% (Delta tau_pq, Delta nu_pq) through the cross-target moments
%   Lambda_{p,q}^{(a,b)} = (df)^a (Ts)^b F_a(-df*Dtau_pq, N) F_b(Ts*Dnu_pq, M),
%   F_a(z,K) = sum_{k=0}^{K-1} k^a e^{j2pi k z}        (Dirichlet moment),
% with Dtau_pq = tau_q - tau_p,  Dnu_pq = nu_q - nu_p.  The same formula gives
% the (symmetric, tau/nu-independent) diagonal blocks when p = q.
%
% Inputs:
%   N, M    : grid size (subcarriers x OFDM symbols)
%   deltaf  : subcarrier spacing [Hz]
%   Ts      : OFDM symbol duration including CP [s]
%   sigma2  : noise variance
%   alpha   : P x 1 complex path gains
%   tau     : P x 1 path delays   [s]
%   nu      : P x 1 path Dopplers [Hz]
%
% Output:
%   J       : 4P x 4P Fisher information matrix.

    P = numel(alpha);

    c2 = 2/sigma2;
    c4 = 4*pi/sigma2;
    c8 = 8*pi^2/sigma2;

    J = zeros(4*P, 4*P);

    for p = 1:P
        for q = 1:P
            ap = alpha(p);  aq = alpha(q);

            zN = -deltaf * (tau(q) - tau(p));   % argument of the N-sum
            zM =  Ts     * (nu(q)  - nu(p));    % argument of the M-sum

            F0n = dirichlet_moment(0, zN, N);
            F1n = dirichlet_moment(1, zN, N);
            F2n = dirichlet_moment(2, zN, N);
            F0m = dirichlet_moment(0, zM, M);
            F1m = dirichlet_moment(1, zM, M);
            F2m = dirichlet_moment(2, zM, M);

            L00 = F0n * F0m;
            L10 = deltaf    * F1n * F0m;
            L01 = Ts        * F0n * F1m;
            L20 = deltaf^2  * F2n * F0m;
            L02 = Ts^2      * F0n * F2m;
            L11 = deltaf*Ts * F1n * F1m;

            B = zeros(4,4);
            % alpha-alpha (eqs. 234-236)
            B(1,1) =  c2*real(L00);   B(2,2) =  c2*real(L00);
            B(1,2) = -c2*imag(L00);   B(2,1) =  c2*imag(L00);
            % alpha_p (rows) vs tau_q/nu_q (cols)  (eqs. 237-240)
            B(1,3) =  c4*imag(aq*L10);   B(1,4) = -c4*imag(aq*L01);
            B(2,3) = -c4*real(aq*L10);   B(2,4) =  c4*real(aq*L01);
            % tau_p/nu_p (rows) vs alpha_q (cols)  (eqs. 241-244)
            B(3,1) = -c4*imag(conj(ap)*L10);   B(3,2) = -c4*real(conj(ap)*L10);
            B(4,1) =  c4*imag(conj(ap)*L01);   B(4,2) =  c4*real(conj(ap)*L01);
            % tau/nu vs tau/nu  (eqs. 245-248)
            B(3,3) =  c8*real(conj(ap)*aq*L20);
            B(4,4) =  c8*real(conj(ap)*aq*L02);
            B(3,4) = -c8*real(conj(ap)*aq*L11);
            B(4,3) = -c8*real(conj(ap)*aq*L11);

            J(4*(p-1)+(1:4), 4*(q-1)+(1:4)) = B;
        end
    end
end

function F = dirichlet_moment(a, z, K)
% F_a(z,K) = sum_{k=0}^{K-1} k^a e^{j2pi k z}  (exact finite sum; a=0,1,2).
    k = 0:K-1;
    F = sum( (k.^a) .* exp(1j*2*pi*k*z) );
end

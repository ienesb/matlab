function [A, B] = getAB_OFDM_multi_noICI_expected(N, M, deltaf, Ts, sigma2, alpha, tau, nu, p)
% Multi-target EXPECTED MCRB component matrices A and B for OFDM, no-ICI
% channel, i.i.d. unit-energy QPSK with BSC bit-error rate p (ofdm_mcrb.tex
% Sec. 1.2.1 / 2.2.1).  Parameter stacking eta = [eta_1;...;eta_P],
% eta_p = [aR_p, aI_p, tau_p, nu_p].
%
%   E{A} = -E{Jtilde} + E{A_Delta},   E{B} = E{Jtilde} + E{b}E{b}^T.
% E{Jtilde} is the expected no-ICI multi FIM.  A_Delta is block-diagonal, so the
% off-diagonal blocks of E{A} are just -E{Jtilde}_{p,q}; the diagonal blocks use
%   Sigma_p^{(a,b)} = sum_r alpha_r^*  E{Gamma_{r,p}^{(a,b)}},
%   E{Gamma_{r,p}^{(a,b)}} = -2p df^a Ts^b F_a(-df*(tau_p-tau_r),N) F_b(Ts*(nu_p-nu_r),M).
%
% NOTE 1 (doc gap): exact E{B} also needs the bias covariance Cov(b), absent from
% ofdm_mcrb.tex; this returns the mean-bias E{Jtilde} + E{b}E{b}^T.
% NOTE 2 (doc typo): eqs.(359)-(366) carry explicit df,Ts factors on top of the
% (n df)^a (m Ts)^b already inside Gamma_{r,p}; that double counts.  Here the
% moments include df^a Ts^b and the block formulas omit the extra factors, which
% reduces correctly to the verified single-target case at P=1.
%
% Inputs: grid (N,M), deltaf, Ts, sigma2, alpha/tau/nu (P x 1), BSC p.
% Output: 4P x 4P A, B.

    P = numel(alpha);
    J = getFIM_OFDM_multi_noICI_expected(N, M, deltaf, Ts, sigma2, alpha, tau, nu);

    A = -J;                          % off-diagonal blocks are final (-Jtilde_{p,q})
    b = zeros(4*P, 1);

    c2 = 2/sigma2;  c4 = 4*pi/sigma2;  c8 = 8*pi^2/sigma2;

    for q = 1:P                      % target whose 4x4 diagonal block we build
        % Sigma_q^{(a,b)} = -2p sum_r alpha_r^* df^a Ts^b F_a(..) F_b(..)
        Sg = @(a,bb) sigma_moment(a, bb, q, alpha, tau, nu, deltaf, Ts, N, M, p);
        S00 = Sg(0,0); S10 = Sg(1,0); S01 = Sg(0,1);
        S20 = Sg(2,0); S02 = Sg(0,2); S11 = Sg(1,1);
        aq = alpha(q);

        Ad = zeros(4);
        Ad(1,3) =  c4*imag(S10);   Ad(1,4) = -c4*imag(S01);
        Ad(2,3) =  c4*real(S10);   Ad(2,4) = -c4*real(S01);
        Ad(3,3) = -c8*real(aq*S20);
        Ad(4,4) = -c8*real(aq*S02);
        Ad(3,4) =  c8*real(aq*S11);
        Ad(3,1)=Ad(1,3); Ad(4,1)=Ad(1,4); Ad(3,2)=Ad(2,3); Ad(4,2)=Ad(2,4); Ad(4,3)=Ad(3,4);

        idx = 4*(q-1)+(1:4);
        A(idx, idx) = A(idx, idx) + Ad;

        b(idx) = [ c2*real(S00); -c2*imag(S00); c4*imag(aq*S10); -c4*imag(aq*S01) ];
    end

    % Exact E{B} = E{Jtilde} + E{b}E{b}^T + Cov(b).
    Cb = getMCRB_covb_OFDM(N, M, deltaf, Ts, sigma2, alpha, tau, nu, p, false);
    B = J + b*b.' + Cb;
end

function S = sigma_moment(a, bb, q, alpha, tau, nu, deltaf, Ts, N, M, p)
    S = 0;
    for r = 1:numel(alpha)
        Fa = deltaf^a * dirichlet_moment(a, -deltaf*(tau(q)-tau(r)), N);
        Fb = Ts^bb    * dirichlet_moment(bb, Ts*(nu(q)-nu(r)), M);
        S = S + conj(alpha(r)) * Fa * Fb;
    end
    S = -2*p * S;
end

function F = dirichlet_moment(a, z, K)
    k = 0:K-1;
    F = sum( (k.^a) .* exp(1j*2*pi*k*z) );
end

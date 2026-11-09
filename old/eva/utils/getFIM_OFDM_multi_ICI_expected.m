function J = getFIM_OFDM_multi_ICI_expected(N, M, deltaf, Ts, sigma2, alpha, tau, nu)
% Multi-target EXPECTED OFDM Fisher information matrix, ICI channel.
% Closed form for i.i.d. unit-energy, zero-mean QPSK symbols, i.e. the
% data-averaged FIM of Section 1.2.2 / eqs. (107)-(117) of ofdm_fim.pdf.
% Parameter stacking:
%   eta = [eta_1; ...; eta_P],  eta_p = [alphaR_p, alphaI_p, tau_p, nu_p].
%
% In expectation the cross-bin data products collapse onto n'=n'', leaving
% the data-independent cross-target moments (Dtau_pq = tau_q - tau_p):
%   S_a^{pq}  = sum_{n,n'} (n')^a Ip*(n'-n) Iq(n'-n)  e^{-j2pi n' df Dtau_pq}
%   T_a^{pq}  = sum_{n,n'} (n')^a Ip*(n'-n) Iq'(n'-n) e^{-j2pi n' df Dtau_pq}
%   Tc_a^{pq} = sum_{n,n'} (n')^a Ip'*(n'-n) Iq(n'-n) e^{-j2pi n' df Dtau_pq}
%   U^{pq}    = sum_{n,n'}       Ip'*(n'-n) Iq'(n'-n) e^{-j2pi n' df Dtau_pq}
% and the Doppler Dirichlet sum F_b^{pq} = sum_m m^b e^{j2pi m Ts Dnu_pq}.
% Here I' := (1/j2pi) dI/dnu (PDF eq. 8-9); get_ici_coeffs_derivative returns
% the full dI/dnu = j2pi*I', hence the 1/(j2pi) scaling applied below.
%
% At p=q this reduces to the single-target ICI expected FIM
% (getFIM_OFDM_ICI_expected); as I->delta it reduces to the no-ICI multi FIM.
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
    T = 1/deltaf;

    [NP, NN] = ndgrid(0:N-1, 0:N-1);
    diff_n   = NP - NN;
    npcol    = (0:N-1).';

    % Per-target kernels and PDF-convention derivatives I' = (1/j2pi) dI/dnu.
    Imat  = cell(P,1);
    Ipmat = cell(P,1);
    for p = 1:P
        Imat{p}  = get_ici_coeffs(diff_n, nu(p), T);
        Ipmat{p} = get_ici_coeffs_derivative(diff_n, nu(p), T) / (1j*2*pi);
    end

    c2 = 2/sigma2;
    c4 = 4*pi/sigma2;
    c8 = 8*pi^2/sigma2;

    J = zeros(4*P, 4*P);

    for p = 1:P
        for q = 1:P
            ap = alpha(p);  aq = alpha(q);

            % n'-phase from the delay separation (rows = n').
            wph = exp(-1j*2*pi*deltaf*npcol*(tau(q)-tau(p)));   % N x 1
            w0 = wph;  w1 = npcol.*wph;  w2 = (npcol.^2).*wph;

            % Row sums over n (dim 2) of the four kernel products.
            rS  = sum(conj(Imat{p})  .* Imat{q},  2);   % I_p* I_q
            rT  = sum(conj(Imat{p})  .* Ipmat{q}, 2);   % I_p* I_q'
            rTc = sum(conj(Ipmat{p}) .* Imat{q},  2);   % I_p'* I_q
            rU  = sum(conj(Ipmat{p}) .* Ipmat{q}, 2);   % I_p'* I_q'

            S0 = sum(w0.*rS);   S1 = sum(w1.*rS);   S2 = sum(w2.*rS);
            T0 = sum(w0.*rT);   T1 = sum(w1.*rT);
            Tc0 = sum(w0.*rTc); Tc1 = sum(w1.*rTc);
            Uu = sum(w0.*rU);

            % Doppler Dirichlet moments F_b^{pq}.
            zM  = Ts*(nu(q)-nu(p));
            F0 = dirichlet_moment(0, zM, M);
            F1 = dirichlet_moment(1, zM, M);
            F2 = dirichlet_moment(2, zM, M);

            % Expected inner products (eqs. 111-115).
            PhiPhi = S0 * F0;
            PhiPsi = deltaf   * S1 * F0;      % E<Phi_p,Psi_q>
            PsiPhi = deltaf   * S1 * F0;      % E<Psi_p,Phi_q>
            PsiPsi = deltaf^2 * S2 * F0;
            PhiOm  = Ts*S0*F1 + T0*F0;        % E<Phi_p,Omega_q>
            OmPhi  = Ts*S0*F1 + Tc0*F0;       % E<Omega_p,Phi_q>
            PsiOm  = deltaf*( Ts*S1*F1 + T1*F0 );
            OmPsi  = deltaf*( Ts*S1*F1 + Tc1*F0 );
            OmOm   = Ts^2*S0*F2 + Ts*(T0+Tc0)*F1 + Uu*F0;

            % Block assembly (eqs. 99-106).
            B = zeros(4,4);
            B(1,1) =  c2*real(PhiPhi);   B(2,2) =  c2*real(PhiPhi);
            B(1,2) = -c2*imag(PhiPhi);   B(2,1) =  c2*imag(PhiPhi);

            B(1,3) =  c4*imag(aq*PhiPsi);   B(1,4) = -c4*imag(aq*PhiOm);
            B(2,3) = -c4*real(aq*PhiPsi);   B(2,4) =  c4*real(aq*PhiOm);

            B(3,1) = -c4*imag(conj(ap)*PsiPhi);   B(3,2) = -c4*real(conj(ap)*PsiPhi);
            B(4,1) =  c4*imag(conj(ap)*OmPhi);    B(4,2) =  c4*real(conj(ap)*OmPhi);

            B(3,3) =  c8*real(conj(ap)*aq*PsiPsi);
            B(4,4) =  c8*real(conj(ap)*aq*OmOm);
            B(3,4) = -c8*real(conj(ap)*aq*PsiOm);
            B(4,3) = -c8*real(conj(ap)*aq*OmPsi);

            J(4*(p-1)+(1:4), 4*(q-1)+(1:4)) = B;
        end
    end
end

function F = dirichlet_moment(a, z, K)
% F_a(z,K) = sum_{k=0}^{K-1} k^a e^{j2pi k z}  (exact finite sum; a=0,1,2).
    k = 0:K-1;
    F = sum( (k.^a) .* exp(1j*2*pi*k*z) );
end

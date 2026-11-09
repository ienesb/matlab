function Jbar = fim_expected(wf, g, cr, use_ici, sigma2)
% FIM_EXPECTED  Closed-form data-averaged FIM under i.i.d. unit-energy,
% zero-mean, zero-pseudo-variance symbols (QPSK / PSK / QAM after power
% normalization).
%
%   Jbar = fim_expected(wf, g, cr, use_ici, sigma2)      4P x 4P
%
% 'cp-ofdm', 'cp-otfs' : identical (ISFFT/SFFT unitarity), single and multi
%     target, no-ICI (Dirichlet cross moments Lambda^{(a,b)}, cp_ofdm_fim.tex
%     eqs. 66-80) and ICI (cross-target S_a^{pq}, T_a^{pq}, Tcheck_a^{pq}, U^{pq}
%     with the Doppler Dirichlet sums F_b, cp_otfs_fim.tex multi-target ICI).
% 'rcp-otfs'           : SAMPLED-branch single-target closed form
%     (rcp_otfs_fim.tex eq. RCP_Obar_moment): the nu entries carry the
%     positive corrections 2 Ts T P1 Q1 (O(1/M)) and Ts^2 M P2 (O(1/M^2)) when
%     use_ici = true (intra-frame phase kept); use_ici = false gives the CP
%     form with TS -> T.  Requires P = 1 and a sample-spaced delay.
% 'rcp-ofdm', 'rcp-otfs-cont' : no closed form (Monte Carlo only) -> error.
%
% Ported from old/eva/utils/getFIM_OFDM_multi_{noICI,ICI}_expected.m.
    if nargin < 5 || isempty(sigma2), sigma2 = 1; end
    wf = lower(wf);
    N = g.N;  M = g.M;  df = g.deltaf;
    P = numel(cr.alpha);
    c2 = 2/sigma2;  c4 = 4*pi/sigma2;  c8 = 8*pi^2/sigma2;

    switch wf
        case {'cp-ofdm', 'cp-otfs'}
            TS = g.TS;
            if use_ici
                [NP, NN] = ndgrid(0:N-1, 0:N-1);
                Imat = cell(P,1);  Ipmat = cell(P,1);
                for p = 1:P
                    [Imat{p}, Ipmat{p}] = ici_kernel(NP - NN, cr.nu(p), g.T);
                end
                npcol = (0:N-1).';
            end
            Jbar = zeros(4*P);
            for p = 1:P
                for q = 1:P
                    ap = cr.alpha(p);  aq = cr.alpha(q);
                    zM = TS*(cr.nu(q) - cr.nu(p));
                    F0m = dirichlet_moment(0, zM, M);  F1m = dirichlet_moment(1, zM, M);  F2m = dirichlet_moment(2, zM, M);
                    if ~use_ici
                        zN = -df*(cr.tau(q) - cr.tau(p));
                        F0n = dirichlet_moment(0, zN, N);  F1n = dirichlet_moment(1, zN, N);  F2n = dirichlet_moment(2, zN, N);
                        PhiPhi = F0n*F0m;
                        PhiPsi = df*F1n*F0m;      PsiPhi = PhiPsi;
                        PsiPsi = df^2*F2n*F0m;
                        PhiOm  = TS*F0n*F1m;      OmPhi  = PhiOm;
                        PsiOm  = df*TS*F1n*F1m;   OmPsi  = PsiOm;
                        OmOm   = TS^2*F0n*F2m;
                    else
                        wph = exp(-1j*2*pi*df*npcol*(cr.tau(q) - cr.tau(p)));
                        w0 = wph;  w1 = npcol.*wph;  w2 = (npcol.^2).*wph;
                        rS  = sum(conj(Imat{p})  .* Imat{q},  2);
                        rT  = sum(conj(Imat{p})  .* Ipmat{q}, 2);
                        rTc = sum(conj(Ipmat{p}) .* Imat{q},  2);
                        rU  = sum(conj(Ipmat{p}) .* Ipmat{q}, 2);
                        S0 = sum(w0.*rS);  S1 = sum(w1.*rS);  S2 = sum(w2.*rS);
                        T0 = sum(w0.*rT);  T1 = sum(w1.*rT);
                        Tc0 = sum(w0.*rTc); Tc1 = sum(w1.*rTc);
                        Uu = sum(w0.*rU);
                        PhiPhi = S0*F0m;
                        PhiPsi = df*S1*F0m;       PsiPhi = PhiPsi;
                        PsiPsi = df^2*S2*F0m;
                        PhiOm  = TS*S0*F1m + T0*F0m;
                        OmPhi  = TS*S0*F1m + Tc0*F0m;
                        PsiOm  = df*(TS*S1*F1m + T1*F0m);
                        OmPsi  = df*(TS*S1*F1m + Tc1*F0m);
                        OmOm   = TS^2*S0*F2m + TS*(T0 + Tc0)*F1m + Uu*F0m;
                    end
                    B = zeros(4);
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
                    Jbar(4*(p-1)+(1:4), 4*(q-1)+(1:4)) = B;
                end
            end

        case 'rcp-otfs'
            assert(P == 1, 'fim_expected:rcp', 'RCP-OTFS closed form is single-target only (P = %d).', P);
            assert(abs(cr.tau/g.Ts - round(cr.tau/g.Ts)) < 1e-6, 'fim_expected:rcp', 'RCP-OTFS closed form needs tau = l_p Ts.');
            T = g.T;  Ts = g.Ts;
            if ~use_ici, Ts = 0; end          % intra-frame phase dropped -> CP form with TS -> T
            a = cr.alpha;  aR = real(a);  aI = imag(a);  a2 = abs(a)^2;
            P1 = power_sums(N,1);  P2 = power_sums(N,2);  Q1 = power_sums(M,1);  Q2 = power_sums(M,2);
            Jbar = zeros(4);
            Jbar(1,1) = c2*N*M;  Jbar(2,2) = c2*N*M;
            Jbar(3,3) = c8*a2*df^2*M*P2;
            Jbar(4,4) = c8*a2*(T^2*N*Q2 + 2*Ts*T*P1*Q1 + Ts^2*M*P2);
            Jbar(3,4) = -c8*a2*df*(T*P1*Q1 + Ts*M*P1^2/N);   Jbar(4,3) = Jbar(3,4);
            Jbar(1,3) =  c4*aI*df*M*P1;   Jbar(3,1) = Jbar(1,3);
            Jbar(2,3) = -c4*aR*df*M*P1;   Jbar(3,2) = Jbar(2,3);
            Jbar(1,4) = -c4*aI*(T*N*Q1 + Ts*M*P1);   Jbar(4,1) = Jbar(1,4);
            Jbar(2,4) =  c4*aR*(T*N*Q1 + Ts*M*P1);   Jbar(4,2) = Jbar(2,4);

        case {'rcp-ofdm', 'rcp-otfs-cont'}
            error('fim_expected:none', 'No closed-form expected FIM for %s (rcp_ofdm_fim.tex prescribes Monte Carlo).', wf);
        otherwise
            error('fim_expected:wf', 'unknown waveform %s', wf);
    end
end

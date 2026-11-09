function [A, B, Jbar, bbar, Cb] = mcrb_AB_expected(wf, g, cr, use_ici, sigma2, p)
% MCRB_AB_EXPECTED  Closed-form expected MCRB matrices under i.i.d. unit-energy
% Gray-QPSK data with a BSC of rate p (latex/mcrb/cp_ofdm_mcrb.tex,
% cp_otfs_mcrb.tex; review Finding 5 included):
%   E{A} = -E{Jt} + E{A_Delta},   E{B} = E{Jt} + E{b}E{b}^T + Cov(b)
% 'cp-ofdm' and 'cp-otfs' share E{Jt}, E{b}, E{A_Delta} (unitarity); only Cov(b)
% differs (covb_cp_ofdm / covb_cp_otfs).  No-ICI diagonal blocks use the Dirichlet
% moments Sigma_q^{(a,b)}; ICI blocks use S_a^{rq}, T_a^{rq}, R^{rq}, F_b^{rq}
% with the (nu,nu) entry built from R = sum I* I'' (not U).  Oracle only.
% Ported from old/eva/utils/getAB_OFDM_multi_{noICI,ICI}_expected.m.
    N = g.N;  M = g.M;  df = g.deltaf;  TS = g.TS;
    alpha = cr.alpha(:);  tau = cr.tau(:);  nu = cr.nu(:);
    P = numel(alpha);
    Jbar = fim_expected('cp-ofdm', g, cr, use_ici, sigma2);
    A = -Jbar;  bbar = zeros(4*P, 1);
    c2 = 2/sigma2;  c4 = 4*pi/sigma2;  c8 = 8*pi^2/sigma2;
    if use_ici
        [NP, NN] = ndgrid(0:N-1, 0:N-1);  npcol = (0:N-1).';
        Imat = cell(P,1); Ipmat = cell(P,1); Ippmat = cell(P,1);
        for t = 1:P, [Imat{t}, Ipmat{t}, Ippmat{t}] = ici_kernel(NP - NN, nu(t), g.T); end
    end
    for q = 1:P
        aq = alpha(q);
        if ~use_ici
            Sg = @(a, bb) sigma_moment(a, bb, q, alpha, tau, nu, df, TS, N, M, p);
            S00 = Sg(0,0); S10 = Sg(1,0); S01 = Sg(0,1); S20 = Sg(2,0); S02 = Sg(0,2); S11 = Sg(1,1);
            IPphi = S00; IPpsi = S10; IPom = S01; IPxi = S20; IPth = S02; IPom1 = S11;
        else
            IPphi = 0; IPpsi = 0; IPxi = 0; IPom = 0; IPom1 = 0; IPth = 0;
            for r = 1:P
                wph = exp(-1j*2*pi*df*npcol*(tau(q) - tau(r)));
                rS = sum(conj(Imat{r}).*Imat{q}, 2);  rT = sum(conj(Imat{r}).*Ipmat{q}, 2);  rR = sum(conj(Imat{r}).*Ippmat{q}, 2);
                S0 = sum(wph.*rS);  S1 = sum(npcol.*wph.*rS);  S2 = sum(npcol.^2.*wph.*rS);
                T0 = sum(wph.*rT);  T1 = sum(npcol.*wph.*rT);  R0 = sum(wph.*rR);
                zM = TS*(nu(q) - nu(r));
                F0 = dirichlet_moment(0, zM, M);  F1 = dirichlet_moment(1, zM, M);  F2 = dirichlet_moment(2, zM, M);
                w = conj(alpha(r))*(-2*p);
                IPphi = IPphi + w*S0*F0;
                IPpsi = IPpsi + w*df*S1*F0;
                IPxi  = IPxi  + w*df^2*S2*F0;
                IPom  = IPom  + w*(TS*S0*F1 + T0*F0);
                IPom1 = IPom1 + w*df*(TS*S1*F1 + T1*F0);
                IPth  = IPth  + w*(TS^2*S0*F2 + 2*TS*T0*F1 + R0*F0);
            end
        end
        Ad = zeros(4);
        Ad(1,3) =  c4*imag(IPpsi);   Ad(1,4) = -c4*imag(IPom);
        Ad(2,3) =  c4*real(IPpsi);   Ad(2,4) = -c4*real(IPom);
        Ad(3,3) = -c8*real(aq*IPxi);
        Ad(4,4) = -c8*real(aq*IPth);
        Ad(3,4) =  c8*real(aq*IPom1);
        Ad(3,1) = Ad(1,3); Ad(4,1) = Ad(1,4); Ad(3,2) = Ad(2,3); Ad(4,2) = Ad(2,4); Ad(4,3) = Ad(3,4);
        idx = 4*(q-1)+(1:4);
        A(idx, idx) = A(idx, idx) + Ad;
        bbar(idx) = [c2*real(IPphi); -c2*imag(IPphi); c4*imag(aq*IPpsi); -c4*imag(aq*IPom)];
    end
    switch lower(wf)
        case 'cp-ofdm', Cb = covb_cp_ofdm(g, cr, use_ici, sigma2, p);
        case 'cp-otfs', Cb = covb_cp_otfs(g, cr, use_ici, sigma2, p);
        otherwise, error('mcrb_AB_expected:wf', 'expected MCRB exists for cp-ofdm / cp-otfs only');
    end
    B = Jbar + bbar*bbar.' + Cb;
end

function S = sigma_moment(a, bb, q, alpha, tau, nu, df, TS, N, M, p)
    S = 0;
    for r = 1:numel(alpha)
        Fa = df^a * dirichlet_moment(a, -df*(tau(q) - tau(r)), N);
        Fb = TS^bb * dirichlet_moment(bb, TS*(nu(q) - nu(r)), M);
        S = S + conj(alpha(r))*Fa*Fb;
    end
    S = -2*p*S;
end

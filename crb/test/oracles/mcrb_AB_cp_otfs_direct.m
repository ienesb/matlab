function [A, B] = mcrb_AB_cp_otfs_direct(x, xhat, g, cr, sigma2)
% MCRB_AB_CP_OTFS_DIRECT  Independent delay-Doppler-domain CP-OTFS MCRB A, B
% (no OFDM-side ICI) built from the Dirichlet kernels A_p, A_p', A_p'', B_p,
% B_p', B_p'' of latex/mcrb/cp_otfs_mcrb.tex (gradients lines 52-57, Hessian
% lines 72-78), evaluated at cr (eta_eval = eta_true).  Test oracle only.
    [N, M] = size(x);  P = numel(cr.alpha);
    ns = (0:N-1).';  ms = (0:M-1).';
    [Lp, L] = ndgrid(0:N-1, 0:N-1);  [K, Kp] = ndgrid(0:M-1, 0:M-1);
    G = zeros(N*M, 4*P);  mu_true = zeros(N*M, 1);  mu_as = zeros(N*M, 1);
    Hf = cell(P, 1);
    for p = 1:P
        a = cr.alpha(p);  tau = cr.tau(p);  nu = cr.nu(p);
        kerA = @(w) arrayfun(@(d) sum(w(ns) .* exp(1j*2*pi*ns*(d/N - g.deltaf*tau))), Lp - L);
        kerB = @(w) arrayfun(@(d) sum(w(ms) .* exp(1j*2*pi*ms*(d/M + g.TS*nu))), K - Kp);
        Am = kerA(@(n) ones(size(n)));  A1 = kerA(@(n) n);  A2 = kerA(@(n) n.^2);
        Bm = kerB(@(m) ones(size(m)));  B1 = kerB(@(m) m);  B2 = kerB(@(m) m.^2);
        f = @(S, Ad, Bd) Ad*S*Bd/(N*M);
        Phi = f(xhat, Am, Bm);  Psi = f(xhat, A1, Bm);  Om = f(xhat, Am, B1);
        Xi = f(xhat, A2, Bm);   Om1 = f(xhat, A1, B1);  Th = f(xhat, Am, B2);
        col = 4*(p-1);
        G(:, col+1) = Phi(:);  G(:, col+2) = 1j*Phi(:);
        G(:, col+3) = -1j*2*pi*g.deltaf*a*Psi(:);  G(:, col+4) = 1j*2*pi*g.TS*a*Om(:);
        mu_as = mu_as + a*Phi(:);
        mu_true = mu_true + a*reshape(f(x, Am, Bm), [], 1);
        Hf{p} = struct('Psi', Psi, 'Om', Om, 'Xi', Xi, 'Om1', Om1, 'Th', Th);
    end
    Delta = mu_true - mu_as;
    Jt = (2/sigma2)*real(G'*G);  b = (2/sigma2)*real(G'*Delta);
    B = Jt + b*b.';
    Ad = zeros(4*P);
    ct = @(h) (2/sigma2)*real(Delta'*h(:));
    df = g.deltaf;  TS = g.TS;
    for p = 1:P
        a = cr.alpha(p);  H = Hf{p};  i1 = 4*(p-1)+1; i2 = i1+1; i3 = i1+2; i4 = i1+3;
        Ad(i1,i3) = ct(-1j*2*pi*df*H.Psi);   Ad(i3,i1) = Ad(i1,i3);
        Ad(i1,i4) = ct( 1j*2*pi*TS*H.Om);    Ad(i4,i1) = Ad(i1,i4);
        Ad(i2,i3) = ct( 2*pi*df*H.Psi);      Ad(i3,i2) = Ad(i2,i3);
        Ad(i2,i4) = ct(-2*pi*TS*H.Om);       Ad(i4,i2) = Ad(i2,i4);
        Ad(i3,i3) = ct(-4*pi^2*df^2*a*H.Xi);
        Ad(i4,i4) = ct(-4*pi^2*TS^2*a*H.Th);
        Ad(i3,i4) = ct( 4*pi^2*df*TS*a*H.Om1);  Ad(i4,i3) = Ad(i3,i4);
    end
    A = -Jt + Ad;
end

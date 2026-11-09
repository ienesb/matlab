function Cb = covb_cp_ofdm(g, cr, use_ici, sigma2, p)
% COVB_CP_OFDM  Bias-vector covariance Cov(b) = E{b b^T} - E{b}E{b}^T for the
% CP-OFDM MCRB under i.i.d. unit-energy Gray-QPSK with BSC rate p
% (latex/mcrb/cp_ofdm_mcrb.tex Sec. "Bias covariance").  Column-wise bilinear
% forms S_k = sum_m e_m^H K_k^{(m)} xhat_m with data-independent K; with the
% per-symbol moments E[e* xhat] = -2p, E|e|^2 = 4p, E[(e* xhat)^2] = 4p^2:
%   Cov(S_k,S_l)  = 4p <K_k,K_l>_F - 4p^2 sum_n K_k(n,n) conj(K_l(n,n))
%   Cov~(S_k,S_l) = 4p^2 ( tr(K_k K_l) - sum_n K_k(n,n) K_l(n,n) )
%   [Cov b]_kl    = (2/sigma2^2) sum_m Re{ Cov + Cov~ }
% No-ICI single target collapses to (4p(1-p)|alpha|^2/sigma2) J on every block.
% Oracle only (16 x 8 grids).  Port of old/eva/utils/getMCRB_covb_OFDM.m.
    N = g.N;  M = g.M;  df = g.deltaf;  TS = g.TS;
    alpha = cr.alpha(:);  tau = cr.tau(:);  nu = cr.nu(:);
    P = numel(alpha);
    Acel = cell(P,1); Bcel = cell(P,1); Dcel = cell(P,1);
    if use_ici, [NP, NN] = ndgrid(0:N-1, 0:N-1); end
    for q = 1:P
        rowph = exp(-1j*2*pi*df*(0:N-1)*tau(q));                 % 1 x N over n'
        if use_ici
            [I, Ip] = ici_kernel(NP - NN, nu(q), g.T);             % (n',n) = I(n'-n)
            Aq = I.'  .* rowph;                                    % (n,n') = I(n'-n) e^{-j2pi n' df tau}
            Dq = Ip.' .* rowph;
        else
            Aq = diag(exp(-1j*2*pi*df*(0:N-1).'*tau(q)));
            Dq = zeros(N);
        end
        Bq = Aq .* ((0:N-1)*df);
        Acel{q} = Aq;  Bcel{q} = Bq;  Dcel{q} = Dq;
    end
    GAA = cell(P,P); GAB = cell(P,P); GAD = cell(P,P);
    for pp = 1:P, for q = 1:P
        GAA{pp,q} = Acel{pp}' * Acel{q};  GAB{pp,q} = Acel{pp}' * Bcel{q};  GAD{pp,q} = Acel{pp}' * Dcel{q};
    end, end
    c2pi = 1j*2*pi;
    Cb = zeros(4*P);
    for m = 0:M-1
        K = cell(1, 4*P);
        for q = 1:P
            KaR = zeros(N); Ktau = zeros(N); Knu = zeros(N);
            for pp = 1:P
                w = conj(alpha(pp)) * exp(c2pi*m*TS*(nu(q) - nu(pp)));
                KaR  = KaR  + w*GAA{pp,q};
                Ktau = Ktau + w*GAB{pp,q};
                Knu  = Knu  + w*(m*TS*GAA{pp,q} + GAD{pp,q});
            end
            col = 4*(q-1);
            K{col+1} = KaR;  K{col+2} = 1j*KaR;
            K{col+3} = -c2pi*alpha(q)*Ktau;  K{col+4} = c2pi*alpha(q)*Knu;
        end
        dg = cellfun(@diag, K, 'UniformOutput', false);
        for k = 1:4*P
            for l = k:4*P
                CovH = 4*p*sum(sum(K{k}.*conj(K{l}))) - 4*p^2*sum(dg{k}.*conj(dg{l}));
                CovP = 4*p^2*(sum(sum(K{k}.*(K{l}.'))) - sum(dg{k}.*dg{l}));
                val = (2/sigma2^2)*real(CovH + CovP);
                Cb(k,l) = Cb(k,l) + val;
                if l > k, Cb(l,k) = Cb(l,k) + val; end
            end
        end
    end
end

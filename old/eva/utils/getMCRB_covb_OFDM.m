function Cb = getMCRB_covb_OFDM(N, M, deltaf, Ts, sigma2, alpha, tau, nu, p, use_ici)
% Bias-vector covariance Cov(b) = E{b b^T} - E{b}E{b}^T for the OFDM MCRB
% under i.i.d. unit-energy QPSK with BSC bit-error rate p (ofdm_mcrb.tex
% "Bias covariance" section).  Returns the 4P x 4P real symmetric Cov(b).
%
% Bilinear fourth-moment result: writing S_k = Delta^H gtilde_k = sum_m e_m^H
% K_k^{(m)} x_hat_m (data-independent K), with per-symbol moments
%   E[e* xhat]=-2p, E|e|^2=4p, E[e^2]=E[xhat^2]=E[e xhat]=0, E|xhat|^2=1,
% the covariance of b_k=(2/s2)Re{S_k} is
%   Cov(b_k,b_l) = (2/s4) sum_m Re{ CovH + CovP },
%   CovH = 4p sum_ij K_k(i,j) conj(K_l(i,j)) - 4p^2 sum_n K_k(n,n) conj(K_l(n,n)),
%   CovP = 4p^2 ( tr(K_k K_l) - sum_n K_k(n,n) K_l(n,n) ).
% For no-ICI single target this reduces to Cov(b) = (4p(1-p)|alpha|^2/s2) J.
%
% I' = (1/j2pi) dI/dnu convention.  Inputs mirror getAB_OFDM_*_expected.

    P  = numel(alpha);
    T  = 1/deltaf;
    % Per-target receiver-side kernel matrices A,B,D (rows n, cols n').
    Acel = cell(P,1); Bcel = cell(P,1); Dcel = cell(P,1);
    if use_ici
        [NP, NN] = ndgrid(0:N-1, 0:N-1);
        diff_n = NP - NN;
    end
    for q = 1:P
        rowph = exp(-1j*2*pi*deltaf*(0:N-1)*tau(q));     % 1 x N over n'
        if use_ici
            Im  = get_ici_coeffs(diff_n, nu(q), T).';                    % (n,n') = I(n'-n)
            Ipm = (get_ici_coeffs_derivative(diff_n, nu(q), T)/(1j*2*pi)).';
            Aq = Im  .* rowph;
            Dq = Ipm .* rowph;
        else
            Aq = diag(exp(-1j*2*pi*deltaf*(0:N-1).'*tau(q)));            % I = delta
            Dq = zeros(N);
        end
        Bq = Aq .* ((0:N-1)*deltaf);                                     % column weight n'df
        Acel{q} = Aq;  Bcel{q} = Bq;  Dcel{q} = Dq;
    end

    % m-independent Gram blocks G^{AA}_{pq}=A_p^H A_q, etc.
    GAA = cell(P,P); GAB = cell(P,P); GAD = cell(P,P);
    for pp = 1:P
        for q = 1:P
            GAA{pp,q} = Acel{pp}' * Acel{q};
            GAB{pp,q} = Acel{pp}' * Bcel{q};
            GAD{pp,q} = Acel{pp}' * Dcel{q};
        end
    end

    c2pi = 1j*2*pi;
    Cb = zeros(4*P);

    for m = 0:M-1
        % Build the 4P data-independent K matrices for this column.
        K = cell(1, 4*P);
        for q = 1:P
            KaR = zeros(N); Ktau = zeros(N); Knu = zeros(N);
            for pp = 1:P
                ph = exp(c2pi*m*Ts*(nu(q)-nu(pp)));     % e^{j2pi m Ts (nu_q-nu_p)}
                w  = conj(alpha(pp)) * ph;
                KaR  = KaR  + w * GAA{pp,q};
                Ktau = Ktau + w * GAB{pp,q};
                Knu  = Knu  + w * (m*Ts*GAA{pp,q} + GAD{pp,q});
            end
            col = 4*(q-1);
            K{col+1} = KaR;
            K{col+2} = 1j * KaR;
            K{col+3} = -c2pi*alpha(q) * Ktau;
            K{col+4} =  c2pi*alpha(q) * Knu;
        end

        % Accumulate covariance contributions for this column.
        dg = cell(1,4*P);
        for k = 1:4*P, dg{k} = diag(K{k}); end
        for k = 1:4*P
            for l = k:4*P
                CovH = 4*p   * sum(sum(K{k}.*conj(K{l}))) ...
                     - 4*p^2 * sum(dg{k}.*conj(dg{l}));
                CovP = 4*p^2 * ( sum(sum(K{k}.*(K{l}.'))) ...   % tr(K_k K_l)
                               - sum(dg{k}.*dg{l}) );
                val = (2/sigma2^2) * real(CovH + CovP);
                Cb(k,l) = Cb(k,l) + val;
                if l > k, Cb(l,k) = Cb(l,k) + val; end
            end
        end
    end
end

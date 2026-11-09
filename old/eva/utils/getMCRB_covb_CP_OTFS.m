function Cb = getMCRB_covb_CP_OTFS(N, M, deltaf, Ts, sigma2, alpha, tau, nu, p, use_ici)
% Bias-vector covariance Cov(b) for the CP-OTFS MCRB under i.i.d. unit-energy
% QPSK with a delay-Doppler BSC of rate p.  Direct implementation of
% cp_otfs_mcrb.tex Section 4 ("Bias Covariance and Exact Matrix B").
%
% Each score is one bilinear form over the whole delay-Doppler grid,
%   S_k = e^H N_k xhat,   N_k = sum_r conj(alpha_r) c_{q,k}/(NM)^2 *
%         sum_t ( A_r^H Kdelay_{q,k,t} ) (x) ( B_r^H Kdopp_{q,k,t} ),
% a sum of scaled Kronecker products (delay (x) Doppler).  With the per-symbol
% QPSK+BSC moments (E[e* xhat]=-2p, E|e|^2=4p, E[(e* xhat)^2]=4p^2, ...),
%   Cov(S_k,S_l)      = 4p <N_k,N_l>_F        - 4p^2 sum_i [N_k]ii conj([N_l]ii)
%   Cov~(S_k,S_l)     = 4p^2 ( tr(N_k N_l)    - sum_i [N_k]ii [N_l]ii )
%   [Cov(b)]_kl       = (2/s4) Re{ Cov(S_k,S_l) + Cov~(S_k,S_l) }.
% Every term factorises into a delay-side (N x N) and a Doppler-side (M x M)
% Gram trace, so no NM x NM matrix is ever formed.
%
% Delay kernels: no-ICI A_p(d)=sum_n e^{j2pi n(d/N - df tau_p)} (circulant),
% A_p'(d) its n-weighted version; ICI K_p, K_p^{(1)}, K_p^{I'} (two-index).
% Doppler kernels: B_p(d)=sum_m e^{j2pi m(d/M + Ts nu_p)} (circulant), B_p'.
% For ICI the nu-field Lambda carries two Kronecker terms (Ts K (x) B' and
% K^{I'} (x) B).  I' = (1/j2pi) dI/dnu.

    P = numel(alpha);
    T = 1/deltaf;

    % ---- per-target delay (N x N) and Doppler (M x M) kernel matrices ----
    K0d = cell(P,1); K1d = cell(P,1); KIp = cell(P,1);
    B0  = cell(P,1); B1  = cell(P,1);
    nN = (0:N-1).';  nM = (0:M-1).';
    if use_ici
        [Nrow,Npc] = ndgrid(0:N-1,0:N-1);     % Nrow=n, Npc=n'
        [Npr ,Lc ] = ndgrid(0:N-1,0:N-1);     % Npr=n', Lc=l
        [Lpr ,Nc ] = ndgrid(0:N-1,0:N-1);     % Lpr=l', Nc=n
        Flp = exp(1j*2*pi*Nc.*Lpr/N);         % [l',n]
    end
    for q = 1:P
        if use_ici
            Im  = get_ici_coeffs(Npc-Nrow, nu(q), T);                       % [n,n']=I(n'-n)
            Imp = get_ici_coeffs_derivative(Npc-Nrow, nu(q), T)/(1j*2*pi);  % I'(n'-n)
            h   = exp(-1j*2*pi*Npr.*(deltaf*tau(q) + Lc/N));                % [n',l]
            K0d{q} = Flp * Im  * h;
            K1d{q} = Flp * Im  * ((0:N-1).'.*h);
            KIp{q} = Flp * Imp * h;
        else
            cdel = exp(-1j*2*pi*nN*deltaf*tau(q));
            a0 = N*ifft(cdel);  a1 = N*ifft(nN.*cdel);
            [Lp,Ll] = ndgrid(0:N-1,0:N-1);  Dd = mod(Lp-Ll,N);             % [l',l], d=l'-l
            K0d{q} = a0(Dd+1);  K1d{q} = a1(Dd+1);
        end
        g = exp(1j*2*pi*nM*Ts*nu(q));
        b0 = M*ifft(g);  b1 = M*ifft(nM.*g);
        [Kp,Kk] = ndgrid(0:M-1,0:M-1);  De = mod(Kk-Kp,M);                 % [k',k], d=k-k'
        B0{q} = b0(De+1);  B1{q} = b1(De+1);
    end

    % ---- build N_k term-lists: {weight, delayGram (NxN), doppGram (MxM)} ----
    nq = 4*P;
    Nterms = cell(1, nq);
    sc = 1/(N*M)^2;
    for q = 1:P
        % per-parameter coefficient c and the (delayKer, doppKer, extra) terms
        for kk = 1:4
            switch kk
                case 1, c = 1;                          T_ = {K0d{q}, B0{q}, 1};
                case 2, c = 1j;                         T_ = {K0d{q}, B0{q}, 1};
                case 3, c = -1j*2*pi*deltaf*alpha(q);   T_ = {K1d{q}, B0{q}, 1};
                case 4
                    if use_ici
                        c = 1j*2*pi*alpha(q);
                        T_ = {K0d{q}, B1{q}, Ts ;  KIp{q}, B0{q}, 1};   % Lambda: Ts K(x)B' + K^{I'}(x)B
                    else
                        c = 1j*2*pi*Ts*alpha(q);
                        T_ = {K0d{q}, B1{q}, 1};   % Omega: A (x) B'
                    end
            end
            terms = {};
            for ti = 1:size(T_,1)
                delayKer = T_{ti,1};  doppKer = T_{ti,2};  extra = T_{ti,3};
                for r = 1:P
                    w  = conj(alpha(r)) * c * extra * sc;
                    Pd = K0d{r}' * delayKer;     % A_r^H * delay
                    Pe = B0{r}' * doppKer;       % B_r^H * dopp
                    terms(end+1,:) = {w, Pd, Pe}; %#ok<AGROW>
                end
            end
            Nterms{4*(q-1)+kk} = terms;
        end
    end

    % ---- assemble Cov(b) via factorised trace formula ----
    Cb = zeros(nq);
    for a = 1:nq
        Ta = Nterms{a};
        for b = a:nq
            Tb = Nterms{b};
            CovH = 0; CovP = 0;
            for ia = 1:size(Ta,1)
                wa = Ta{ia,1}; Pda = Ta{ia,2}; Pea = Ta{ia,3};
                dDa = diag(Pda); dEa = diag(Pea);
                for ib = 1:size(Tb,1)
                    wb = Tb{ib,1}; Pdb = Tb{ib,2}; Peb = Tb{ib,3};
                    % Frobenius (Hermitian) part: <Pd_a,Pd_b>_F * <Pe_a,Pe_b>_F
                    frob = sum(sum(Pda.*conj(Pdb))) * sum(sum(Pea.*conj(Peb)));
                    dgH  = sum(dDa.*conj(diag(Pdb))) * sum(dEa.*conj(diag(Peb)));
                    CovH = CovH + 4*p*(wa*conj(wb))*frob - 4*p^2*(wa*conj(wb))*dgH;
                    % complementary (pseudo) part: tr(Pd_a Pd_b) * tr(Pe_a Pe_b)
                    trp  = sum(sum(Pda.*(Pdb.'))) * sum(sum(Pea.*(Peb.')));
                    dgP  = sum(dDa.*diag(Pdb)) * sum(dEa.*diag(Peb));
                    CovP = CovP + 4*p^2*(wa*wb)*(trp - dgP);
                end
            end
            val = (2/sigma2^2) * real(CovH + CovP);
            Cb(a,b) = val;  if b>a, Cb(b,a) = val; end
        end
    end
end

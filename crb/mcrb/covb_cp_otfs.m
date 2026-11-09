function Cb = covb_cp_otfs(g, cr, use_ici, sigma2, p)
% COVB_CP_OTFS  Bias-vector covariance Cov(b) for the CP-OTFS MCRB under i.i.d.
% unit-energy Gray-QPSK with a delay-Doppler BSC of rate p, direct implementation
% of latex/mcrb/cp_otfs_mcrb.tex Sec. 4 ("Bias covariance and exact matrix B"):
% each score is one bilinear form S_k = e^H N_k xhat over the whole DD grid with
% N_k a sum of Kronecker products (delay Gram (x) Doppler Gram), so every trace
% factorizes and no NM x NM matrix is formed.  Oracle only.
% Port of old/eva/utils/getMCRB_covb_CP_OTFS.m with ici_kernel and TS.
    N = g.N;  M = g.M;  df = g.deltaf;  TS = g.TS;  T = g.T;
    alpha = cr.alpha(:);  tau = cr.tau(:);  nu = cr.nu(:);
    P = numel(alpha);
    K0d = cell(P,1); K1d = cell(P,1); KIp = cell(P,1);  B0 = cell(P,1); B1 = cell(P,1);
    nN = (0:N-1).';  nM = (0:M-1).';
    if use_ici
        [Nrow, Npc] = ndgrid(0:N-1, 0:N-1);    % Nrow = n, Npc = n'
        [Npr,  Lc]  = ndgrid(0:N-1, 0:N-1);    % Npr = n', Lc = l
        [Lpr,  Nc]  = ndgrid(0:N-1, 0:N-1);    % Lpr = l', Nc = n
        Flp = exp(1j*2*pi*Nc.*Lpr/N);          % [l', n]
    end
    for q = 1:P
        if use_ici
            [Im, Imp] = ici_kernel(Npc - Nrow, nu(q), T);        % [n,n'] = I(n'-n)
            h = exp(-1j*2*pi*Npr.*(df*tau(q) + Lc/N));            % [n', l]
            K0d{q} = Flp*Im*h;  K1d{q} = Flp*Im*((0:N-1).'.*h);  KIp{q} = Flp*Imp*h;
        else
            cdel = exp(-1j*2*pi*nN*df*tau(q));
            a0 = N*ifft(cdel);  a1 = N*ifft(nN.*cdel);
            [Lp, Ll] = ndgrid(0:N-1, 0:N-1);  Dd = mod(Lp - Ll, N);
            K0d{q} = a0(Dd+1);  K1d{q} = a1(Dd+1);
        end
        gm = exp(1j*2*pi*nM*TS*nu(q));
        b0 = M*ifft(gm);  b1 = M*ifft(nM.*gm);
        [Kp, Kk] = ndgrid(0:M-1, 0:M-1);  De = mod(Kk - Kp, M);
        B0{q} = b0(De+1);  B1{q} = b1(De+1);
    end
    nq = 4*P;  Nterms = cell(1, nq);  sc = 1/(N*M)^2;
    for q = 1:P
        for kk = 1:4
            switch kk
                case 1, c = 1;                          T_ = {K0d{q}, B0{q}, 1};
                case 2, c = 1j;                         T_ = {K0d{q}, B0{q}, 1};
                case 3, c = -1j*2*pi*df*alpha(q);       T_ = {K1d{q}, B0{q}, 1};
                case 4
                    if use_ici, c = 1j*2*pi*alpha(q);   T_ = {K0d{q}, B1{q}, TS; KIp{q}, B0{q}, 1};
                    else,       c = 1j*2*pi*TS*alpha(q); T_ = {K0d{q}, B1{q}, 1}; end
            end
            terms = {};
            for ti = 1:size(T_, 1)
                for r = 1:P
                    w  = conj(alpha(r))*c*T_{ti,3}*sc;
                    terms(end+1, :) = {w, K0d{r}'*T_{ti,1}, B0{r}'*T_{ti,2}}; %#ok<AGROW>
                end
            end
            Nterms{4*(q-1)+kk} = terms;
        end
    end
    Cb = zeros(nq);
    for a = 1:nq
        Ta = Nterms{a};
        for b = a:nq
            Tb = Nterms{b};  CovH = 0;  CovP = 0;
            for ia = 1:size(Ta, 1)
                wa = Ta{ia,1}; Pda = Ta{ia,2}; Pea = Ta{ia,3};  dDa = diag(Pda);  dEa = diag(Pea);
                for ib = 1:size(Tb, 1)
                    wb = Tb{ib,1}; Pdb = Tb{ib,2}; Peb = Tb{ib,3};
                    frob = sum(sum(Pda.*conj(Pdb))) * sum(sum(Pea.*conj(Peb)));
                    dgH  = sum(dDa.*conj(diag(Pdb))) * sum(dEa.*conj(diag(Peb)));
                    CovH = CovH + 4*p*(wa*conj(wb))*frob - 4*p^2*(wa*conj(wb))*dgH;
                    trp  = sum(sum(Pda.*(Pdb.'))) * sum(sum(Pea.*(Peb.')));
                    dgP  = sum(dDa.*diag(Pdb)) * sum(dEa.*diag(Peb));
                    CovP = CovP + 4*p^2*(wa*wb)*(trp - dgP);
                end
            end
            val = (2/sigma2^2)*real(CovH + CovP);
            Cb(a,b) = val;  if b > a, Cb(b,a) = val; end
        end
    end
end

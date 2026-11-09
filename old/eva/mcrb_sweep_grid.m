function mcrb_sweep_grid(grid, SNR_dbs, nTrials, outname, Emax_dd)
% MCRB_SWEEP_GRID  Low-SNR sweep over (N, M, deltaf) to find a grid where OFDM
% beats CP-OTFS on the LB for BOTH delay AND Doppler, despite OFDM's higher BER.
% LB = MCRB + BB (= N0*G + F + BB) from the realistic detector (mcrb_ce_ici core:
% OFDM dftdd-CE + 1-tap MMSE; CP-OTFS threshold-CE + MP).  Reports the ratio
% r = LB_OTFS / LB_OFDM per parameter -- r>1 means OFDM wins.  "OFDM wins both"
% = delay ratio>1 AND Doppler ratio>1 at low SNR.
%
% Hypothesis: small M (short Doppler aperture) starves OTFS's MP of Doppler
% diversity -> more detector failures -> higher OTFS BB on both params -> OFDM
% wins both.  (Grids keep N*deltaf >= ~7.68 MHz so the [0,150,370] ns delays stay
% on distinct taps.)
    thisdir=fileparts(mfilename('fullpath'));
    addpath(fullfile(thisdir,'utils'));
    warning('off','MATLAB:nearlySingularMatrix');

    seed_base=7700;
    if nargin<2||isempty(SNR_dbs), SNR_dbs=0:2:6; end
    if nargin<3||isempty(nTrials), nTrials=200; end
    if nargin<4||isempty(outname), outname='mcrb_sweep_grid_results.mat'; end
    if nargin<5||isempty(Emax_dd), Emax_dd=Inf; end       % DD Doppler-spread truncation (Inf=exact)
    %                              N     M    deltaf   label
    if nargin<1||isempty(grid), grid={ 512,  32, 15e3, 'N512 M32  df15'
           512,  64, 15e3, 'N512 M64  df15'
           512, 128, 15e3, 'N512 M128 df15 (baseline)'
           512, 256, 15e3, 'N512 M256 df15'
           512,  64, 30e3, 'N512 M64  df30 (B=15.4M)'
           512, 128, 30e3, 'N512 M128 df30 (B=15.4M)'
           256, 128, 30e3, 'N256 M128 df30 (B=7.68M)'
           256,  64, 30e3, 'N256 M64  df30 (B=7.68M)' }; end
    nG=size(grid,1);
    res=struct('label',{},'N',{},'M',{},'deltaf',{},'SNR',{}, ...
               'dOF',{},'dOT',{},'nOF',{},'nOT',{},'dr',{},'nr',{},'ser',{});

    fprintf('Low-SNR LB sweep (SNR %s dB, %d trials).  ratio r = LB_OTFS/LB_OFDM (>1 => OFDM wins)\n', ...
            num2str(SNR_dbs), nTrials);
    fprintf('%-28s | %4s | %8s %8s | %8s %8s | %s\n', ...
            'config','SNR','delay_r','dopp_r','ser_OF','ser_OT','OFDM wins both?');
    tAll=tic;
    for gi=1:nG
        N=grid{gi,1}; M=grid{gi,2}; deltaf=grid{gi,3};
        if size(grid,2)>=5, guard_extra=grid{gi,4}; lab=grid{gi,5}; else, guard_extra=[]; lab=grid{gi,4}; end
        [dLB,nLB,ser]=compute_lb(N,M,deltaf,SNR_dbs,nTrials,seed_base,guard_extra,Emax_dd);
        for si=1:numel(SNR_dbs)
            dr=nLB(si,2)/nLB(si,1);  % note: cols are [OFDM OTFS]; ratio OTFS/OFDM
            drt=dLB(si,2)/dLB(si,1);
            win=tern(drt>1 && dr>1,'YES','no');
            fprintf('%-28s | %4d | %8.3f %8.3f | %8.1e %8.1e | %s\n', ...
                    lab, SNR_dbs(si), drt, dr, ser(si,1), ser(si,2), win);
            k=numel(res)+1;
            res(k).label=lab; res(k).N=N; res(k).M=M; res(k).deltaf=deltaf; res(k).SNR=SNR_dbs(si);
            res(k).dOF=dLB(si,1); res(k).dOT=dLB(si,2); res(k).nOF=nLB(si,1); res(k).nOT=nLB(si,2);
            res(k).dr=drt; res(k).nr=dr; res(k).ser=ser(si,:);
        end
        fprintf('  (%.0fs)\n',toc(tAll));
    end
    save(fullfile(thisdir,'results',outname),'res','grid','SNR_dbs','nTrials');
    fprintf('\nSaved results/%s  (Done %.0fs)\n',outname,toc(tAll));

    % summary: configs where OFDM wins both at the lowest SNR
    fprintf('\n=== configs where OFDM wins BOTH (delay & Doppler) at SNR=%d dB ===\n',SNR_dbs(1));
    any_win=false;
    for gi=1:nG
        k=(gi-1)*numel(SNR_dbs)+1;                 % SNR_dbs(1) row for this config
        if res(k).dr>1 && res(k).nr>1
            fprintf('  %-28s  delay_r=%.2f  dopp_r=%.2f\n',res(k).label,res(k).dr,res(k).nr); any_win=true;
        end
    end
    if ~any_win, fprintf('  (none at SNR=%d dB)\n',SNR_dbs(1)); end
end

% ---------------------------------------------------------------------------
function [dLB,nLB,ser]=compute_lb(N,M,deltaf,SNR_dbs,nTrials,seed_base,guard_extra,Emax_dd)
    if nargin<8||isempty(Emax_dd), Emax_dd=Inf; end
% Realistic-detector LB (delay & Doppler) for OFDM and CP-OTFS on the reduced-EVA
% channel at grid (N,M,deltaf).  Returns dLB,nLB,ser as nSNR x 2 ([OFDM OTFS]).
    Ts=1/deltaf; fc=4e9; c=3e8; v=120/3.6; nu_max=v*fc/c;
    M_order=4; nbits=log2(M_order); constellation=qammod((0:M_order-1).',M_order,UnitAveragePower=true);
    eva_delays_ns=[0 150 370]; eva_pow_db=[0 -1.4 -0.6]; P=numel(eva_delays_ns);
    delay_taps=round(eva_delays_ns*1e-9*N*deltaf);
    assert(numel(unique(delay_taps))==P,'delays collide onto same tap at N=%d df=%g',N,deltaf);
    tau=(delay_taps/(N*deltaf)).'; amp=(10.^(eva_pow_db/20)).'; amp=amp/norm(amp);
    l_tau=ceil(max(eva_delays_ns*1e-9)*N*deltaf); k_nu=max(1,ceil(nu_max*M/deltaf));
    nE=4*P; idx_tau=4*(0:P-1)+3; idx_nu=4*(0:P-1)+4;
    n_iter_mp=30; damping=0.3; ce_mode='dftdd'; n_dd=2;

    % OTFS delay guard: G = l_tau + 1 + guard_extra bins.  guard_extra=l_tau -> the
    % full 2*l_tau+1 guard (default); guard_extra=0 -> minimal (l_tau+1) so data sits
    % right after the pilot's channel spread.  G-1>=l_tau always, so the largest
    % delay always fits the guard (pilot response never wraps).
    if nargin<7||isempty(guard_extra), guard_extra=l_tau; end
    guard_extra=max(0,min(round(guard_extra),l_tau));
    G=l_tau+1+guard_extra; l_p=guard_extra; k_p=M/2;
    guard_rows=1:G; ce_rows=(l_p+1):G;
    data_mask_otfs=true(N,M); data_mask_otfs(guard_rows,:)=false;
    otfs_nondata=numel(guard_rows)*M;
    [dN_p,dM_p,n_pilot_cells]=match_pilot_spacing(N,M,otfs_nondata,l_tau,k_nu);
    n_pilots=1:dN_p:N; m_pilots=1:dM_p:M;
    pilot_mask=false(N,M); pilot_mask(n_pilots,m_pilots)=true; data_mask_ofdm=~pilot_mask;
    a_p=1; x_p=sqrt(n_pilot_cells);

    nSNR=numel(SNR_dbs); dLB=zeros(nSNR,2); nLB=zeros(nSNR,2); ser=zeros(nSNR,2);
    for ci=1:2
        modn=tern(ci==1,'ofdm','otfs');
        for si=1:nSNR
            N0=1/db2pow(SNR_dbs(si));
            G_t=zeros(nE,nE,nTrials); F_t=zeros(nE,nE,nTrials); BB_t=zeros(nE,nE,nTrials); be=zeros(1,nTrials);
            parfor t=1:nTrials
                warning('off','MATLAB:nearlySingularMatrix');
                rs=RandStream('Threefry','Seed',seed_base); rs.Substream=10000*si+t; RandStream.setGlobalStream(rs);
                nu=(nu_max*cos(2*pi*rand(1,P)-pi)).'; alpha=amp.*exp(1j*2*pi*rand(P,1));
                eta_bar=packEta(alpha,tau,nu);
                sym=randi(M_order,N,M)-1; Xd=qammod(sym,M_order,UnitAveragePower=true);
                Z=sqrt(N0/2)*(randn(N,M)+1j*randn(N,M));
                if strcmp(modn,'ofdm')
                    d=Xd; d(pilot_mask)=a_p;
                    Y=ofdm_mean(d,deltaf,Ts,alpha,tau,nu,true)+Z;
                    H_hat=ce_ofdm(Y,d,pilot_mask,n_pilots,m_pilots,N,M,l_tau,N0,M_order,ce_mode,n_dd);
                    Xs=conj(H_hat).*Y./(abs(H_hat).^2+N0);
                    dh=d; bh=qamdemod(Xs(data_mask_ofdm),M_order,OutputType="bit",UnitAveragePower=true);
                    dh(data_mask_ofdm)=qammod(bh,M_order,InputType="bit",UnitAveragePower=true);
                    errmask=data_mask_ofdm; Xtf=d; Xhtf=dh;
                    eta0=pseudoTrue_OFDM(d,dh,deltaf,Ts,alpha,tau,nu,true,[],true);
                else
                    d=zeros(N,M); d(data_mask_otfs)=Xd(data_mask_otfs); d(l_p+1,k_p+1)=x_p;
                    Gtrue=build_effective_dd_eva(alpha,delay_taps,nu,N,M,deltaf,Emax_dd);
                    Y=apply_dd_channel_eva(d,Gtrue)+Z;
                    h_ce=zeros(N,M); h_ce(ce_rows,:)=Y(ce_rows,:)/x_p; h_ce(abs(h_ce)<3*sqrt(N0)/x_p)=0;
                    [rc,cc,gc]=find(h_ce); dh=d;
                    if ~isempty(gc)
                        dl_c=mod(rc-(l_p+1),N); dk_c=mod(cc-(k_p+1),M); known=zeros(N,M); known(l_p+1,k_p+1)=x_p;
                        xhat=mp_detector(Y,dl_c,dk_c,gc,N0,constellation,n_iter_mp,data_mask_otfs,known,damping);
                        dh(data_mask_otfs)=xhat(data_mask_otfs);
                    end
                    errmask=data_mask_otfs; Xtf=isfft_u(d); Xhtf=isfft_u(dh);
                    eta0=pseudoTrue_CP_OTFS(d,dh,deltaf,Ts,alpha,tau,nu,true,[],true);
                end
                be(t)=mean(dh(errmask)~=d(errmask));
                [a0,t0,n0]=unpackEta(eta0);
                [A,~,Jt,b]=getAB_OFDM_core(Xtf,Xhtf,deltaf,Ts,1,a0,t0,n0,true,alpha,tau,nu);
                G_t(:,:,t)=stab_quad(A,Jt); F_t(:,:,t)=stab_quad(A,b*b.');
                bias=eta_bar-eta0; BB_t(:,:,t)=bias*bias.';
            end
            LB=N0*mean(G_t,3)+mean(F_t,3)+mean(BB_t,3);
            dLB(si,ci)=sqrt(sum(real(diag(LB(idx_tau,idx_tau)))));
            nLB(si,ci)=sqrt(sum(real(diag(LB(idx_nu, idx_nu )))));
            ser(si,ci)=mean(be);
            fprintf('    compute_lb %s SNR=%2d | dLB=%.3e nLB=%.3e ser=%.2e\n',modn,SNR_dbs(si),dLB(si,ci),nLB(si,ci),ser(si,ci));
            try, save(fullfile(fileparts(mfilename('fullpath')),'results','CKPT_compute_lb.mat'), ...
                      'dLB','nLB','ser','SNR_dbs','N','M','deltaf','ci','si'); catch, end   % crash-safe
        end
    end
end

% ---- local helpers (copied from mcrb_ce_ici.m) ----
function eta=packEta(alpha,tau,nu)
    eta=reshape([real(alpha(:)).'; imag(alpha(:)).'; tau(:).'; nu(:).'],[],1);
end
function Y=stab_quad(A,C)
    d=1./sqrt(abs(diag(A))); d(~isfinite(d))=1; D=diag(d);
    As=D*A*D; Cs=D*C*D; As=(As+As.')/2; Ai=pinv(As); Y=real(D*(Ai*Cs*Ai)*D);
end
function H_hat=ce_ofdm(Y,d,pilot_mask,n_pilots,m_pilots,N,M,l_tau,N0,M_order,ce_mode,n_dd)
    H_ls=Y(n_pilots,m_pilots)./d(n_pilots,m_pilots); [Nq,Mq]=ndgrid(1:N,1:M);
    Afull=exp(-1j*2*pi*(0:N-1).'*(0:l_tau)/N);
    Abas=exp(-1j*2*pi*(n_pilots(:)-1)*(0:l_tau)/N); Hp=Afull*(Abas\H_ls);
    F_H=griddedInterpolant({(1:N).',m_pilots},Hp,'linear','nearest'); H_hat=F_H(Nq,Mq);
    if strcmp(ce_mode,'dftdd')
        for it=1:n_dd
            Xs=conj(H_hat).*Y./(abs(H_hat).^2+N0);
            dhat=qammod(qamdemod(Xs,M_order,UnitAveragePower=true),M_order,UnitAveragePower=true);
            dhat(pilot_mask)=d(pilot_mask); H_hat=Afull*(Afull\(Y./dhat));
        end
    end
end
function s=tern(c,a,b), if c, s=a; else, s=b; end, end

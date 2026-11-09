function diag_marginal_bsc()
% DIAG_MARGINAL_BSC  Decisive test of the delay/Doppler MCRB asymmetry for OFDM.
% Compares the pseudo-true bias floor tr(F+BB) on the delay block and the Doppler
% block under four error PLACEMENTS with the SAME per-trial error count:
%   real    : the actual detector errors (fade-aligned in frequency)
%   iid     : i.i.d. uniform scatter over data cells (the abstract BSC)
%   freqBSC : preserves the per-SUBCARRIER error count r_n, random within each
%             subcarrier over time  (keeps the frequency marginal, kills time struct)
%   timeBSC : preserves the per-SYMBOL error count c_m, random within each symbol
%             over frequency (keeps the time marginal, kills frequency struct)
% Hypothesis: delay bias is set by the FREQUENCY marginal (delay = phase across n),
% so freqBSC ~ real >> iid for delay; Doppler bias is set by the TIME marginal,
% so timeBSC ~ real for Doppler, and since real errors spread ~uniformly in time,
% iid ~ real for Doppler too (small gap).
    thisdir=fileparts(mfilename('fullpath')); addpath(fullfile(thisdir,'utils'));
    warning('off','MATLAB:nearlySingularMatrix');

    N=512; M=128; deltaf=15e3; Ts=1/deltaf; fc=4e9; c=3e8; v=120/3.6;
    nu_max=v*fc/c; M_order=4;
    eva_delays_ns=[0 150 370]; eva_pow_db=[0 -1.4 -0.6]; P=numel(eva_delays_ns);
    delay_taps=round(eva_delays_ns*1e-9*N*deltaf);
    tau=(delay_taps/(N*deltaf)).'; amp=(10.^(eva_pow_db/20)).'; amp=amp/norm(amp);
    l_tau=ceil(max(eva_delays_ns*1e-9)*N*deltaf); k_nu=ceil(nu_max*M/deltaf);
    idx_tau=4*(0:P-1)+3; idx_nu=4*(0:P-1)+4; nE=4*P;

    guard_rows=1:(2*l_tau+1); otfs_nondata=numel(guard_rows)*M;
    [dN_p,dM_p,~]=match_pilot_spacing(N,M,otfs_nondata,l_tau,k_nu);
    n_pilots=1:dN_p:N; m_pilots=1:dM_p:M;
    pilot_mask=false(N,M); pilot_mask(n_pilots,m_pilots)=true; data_mask=~pilot_mask; a_p=1;
    ce_mode='dftdd'; n_dd=2;

    SNRs=[3 6 12]; nTrials=24; seed_base=9000;
    tags={'real','iid','freqBSC','timeBSC'};

    fprintf('OFDM tr(F+BB) delay & Doppler blocks vs error placement (%d trials)\n',nTrials);
    fprintf('%-5s | %-34s | %-34s\n','SNR','DELAY  real  iid  freqBSC timeBSC','DOPPLER real iid  freqBSC timeBSC');
    for si=1:numel(SNRs)
        SNR=SNRs(si); N0=1/db2pow(SNR);
        Dt=nan(nTrials,4); Dn=nan(nTrials,4);
        parfor t=1:nTrials
            warning('off','MATLAB:nearlySingularMatrix');
            rs=RandStream('Threefry','Seed',seed_base); rs.Substream=1000*si+t;
            RandStream.setGlobalStream(rs);
            nu=(nu_max*cos(2*pi*rand(1,P)-pi)).'; alpha=amp.*exp(1j*2*pi*rand(P,1));
            eta_bar=reshape([real(alpha)';imag(alpha)';tau';nu'],[],1);
            sym=randi(M_order,N,M)-1; Xd=qammod(sym,M_order,UnitAveragePower=true);
            Z=sqrt(N0/2)*(randn(N,M)+1j*randn(N,M));
            d=Xd; d(pilot_mask)=a_p;
            Y=ofdm_mean(d,deltaf,Ts,alpha,tau,nu,true)+Z;
            H_hat=ce_ofdm(Y,d,pilot_mask,n_pilots,m_pilots,N,M,l_tau,N0,M_order,ce_mode,n_dd);
            Xs=conj(H_hat).*Y./(abs(H_hat).^2+N0);
            dh=d; dh(data_mask)=qammod(qamdemod(Xs(data_mask),M_order,UnitAveragePower=true),M_order,UnitAveragePower=true);
            E=false(N,M); E(data_mask)=dh(data_mask)~=d(data_mask);
            ne=nnz(E);
            if ne<2, continue; end
            dts=nan(1,4); dns=nan(1,4);
            for pl=1:4
                switch pl
                    case 1, dhx=dh;
                    case 2, dhx=scatter_iid(d,sym,data_mask,ne,M_order);
                    case 3, dhx=scatter_marg(d,sym,data_mask,E,M_order,'freq');
                    case 4, dhx=scatter_marg(d,sym,data_mask,E,M_order,'time');
                end
                eta0=pseudoTrue_OFDM(d,dhx,deltaf,Ts,alpha,tau,nu,true,[],true);
                [a0,t0,n0]=unpackEta(eta0);
                [A,~,~,b]=getAB_OFDM_core(d,dhx,deltaf,Ts,1,a0,t0,n0,true,alpha,tau,nu);
                F=stab_quad(A,b*b.'); bias=eta_bar-eta0; BB=bias*bias.';
                FB=F+BB;
                dts(pl)=sum(real(diag(FB(idx_tau,idx_tau))));
                dns(pl)=sum(real(diag(FB(idx_nu,idx_nu))));
            end
            Dt(t,:)=dts; Dn(t,:)=dns;
        end
        mt=sqrt(nanmedian(Dt)); mn=sqrt(nanmedian(Dn));   % RMSE-like (sqrt of trace)
        fprintf('%3ddB | %.3e %.3e %.3e %.3e | %.3e %.3e %.3e %.3e\n',...
            SNR, mt(1),mt(2),mt(3),mt(4), mn(1),mn(2),mn(3),mn(4));
    end
    fprintf(['\ntags: ' strjoin(tags,'  ') '\n' ...
      'If freqBSC~real and iid<<real for DELAY, the delay bias is a FREQUENCY-marginal\n' ...
      '(delay-signature) effect. If timeBSC~real~iid for DOPPLER, the time axis is\n' ...
      'unstructured, so BSC already captures the Doppler bias.\n']);
end

function dhx=scatter_iid(d,sym,data_mask,ne,M_order)
    dhx=d; di=find(data_mask); pe=di(randperm(numel(di),ne));
    w=mod(sym(pe)+randi(M_order-1,ne,1),M_order);
    dhx(pe)=qammod(w,M_order,UnitAveragePower=true);
end

function dhx=scatter_marg(d,sym,data_mask,E,M_order,dim)
% Keep the per-line error COUNT (freq: per subcarrier row; time: per symbol col)
% but randomize the position within that line and the wrong symbol.
    [N,M]=size(d); dhx=d;
    if strcmp(dim,'freq')
        for n=1:N
            k=nnz(E(n,:)); if k==0, continue; end
            cells=find(data_mask(n,:)); pos=cells(randperm(numel(cells),k));
            for j=pos
                w=mod(sym(n,j)+randi(M_order-1),M_order);
                dhx(n,j)=qammod(w,M_order,UnitAveragePower=true);
            end
        end
    else
        for m=1:M
            k=nnz(E(:,m)); if k==0, continue; end
            cells=find(data_mask(:,m)); pos=cells(randperm(numel(cells),k));
            for i=pos
                w=mod(sym(i,m)+randi(M_order-1),M_order);
                dhx(i,m)=qammod(w,M_order,UnitAveragePower=true);
            end
        end
    end
end

function H_hat=ce_ofdm(Y,d,pilot_mask,n_pilots,m_pilots,N,M,l_tau,N0,M_order,ce_mode,n_dd)
    H_ls=Y(n_pilots,m_pilots)./d(n_pilots,m_pilots); [Nq,Mq]=ndgrid(1:N,1:M);
    Afull=exp(-1j*2*pi*(0:N-1).'*(0:l_tau)/N);
    Abas=exp(-1j*2*pi*(n_pilots(:)-1)*(0:l_tau)/N);
    Hp=Afull*(Abas\H_ls);
    F_H=griddedInterpolant({(1:N).',m_pilots},Hp,'linear','nearest'); H_hat=F_H(Nq,Mq);
    if strcmp(ce_mode,'dftdd')
        for it=1:n_dd
            Xs=conj(H_hat).*Y./(abs(H_hat).^2+N0);
            dhat=qammod(qamdemod(Xs,M_order,UnitAveragePower=true),M_order,UnitAveragePower=true);
            dhat(pilot_mask)=d(pilot_mask); H_hat=Afull*(Afull\(Y./dhat));
        end
    end
end

function Y=stab_quad(A,C)
    d=1./sqrt(abs(diag(A))); d(~isfinite(d))=1; D=diag(d);
    As=D*A*D; Cs=D*C*D; As=(As+As.')/2; Ai=pinv(As); Y=real(D*(Ai*Cs*Ai)*D);
end
function y=nanmedian(x), y=median(x(~isnan(x(:,1)),:),1); end

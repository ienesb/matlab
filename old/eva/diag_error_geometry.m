function diag_error_geometry()
% DIAG_ERROR_GEOMETRY  Where do the OFDM detector's errors land in the (n,m) grid?
% Tests the hypothesis behind the delay/Doppler MCRB asymmetry:
%   delay  is read from the phase across SUBCARRIERS n (frequency),
%   Doppler is read from the phase across SYMBOLS   m (time).
% If the real errors CLUSTER in frequency (few subcarriers, at the delay-set deep
% fades) but SPREAD in time (all symbols alike), then the delay bias is strongly
% error-STRUCTURED (BSC can't mimic it) while the Doppler bias looks BSC-like.
%
% Per trial we measure the normalized participation ratio of the error marginals:
%   PR_n = (sum_n r_n)^2 / (N_data_sc  * sum_n r_n^2),   r_n = #errors at subcarrier n
%   PR_m = (sum_m c_m)^2 / (M_data_sym * sum_m c_m^2),   c_m = #errors at symbol m
% PR = 1 means uniform over that axis; PR -> small means concentrated.  We compare
% the REAL detector to a same-count i.i.d. scatter (BSC null).  We also report the
% mean |H|^2 at error cells relative to the frame average (deep-fade indicator).
    thisdir = fileparts(mfilename('fullpath'));
    addpath(fullfile(thisdir,'utils'));
    warning('off','MATLAB:nearlySingularMatrix');

    N=512; M=128; deltaf=15e3; Ts=1/deltaf; fc=4e9; c=3e8; v=120/3.6;
    nu_max=v*fc/c; M_order=4;
    eva_delays_ns=[0 150 370]; eva_pow_db=[0 -1.4 -0.6]; P=numel(eva_delays_ns);
    delay_taps=round(eva_delays_ns*1e-9*N*deltaf);
    tau=(delay_taps/(N*deltaf)).'; amp=(10.^(eva_pow_db/20)).'; amp=amp/norm(amp);
    l_tau=ceil(max(eva_delays_ns*1e-9)*N*deltaf); k_nu=ceil(nu_max*M/deltaf);

    guard_rows=1:(2*l_tau+1); otfs_nondata=numel(guard_rows)*M;
    [dN_p,dM_p,n_pilot_cells]=match_pilot_spacing(N,M,otfs_nondata,l_tau,k_nu);
    n_pilots=1:dN_p:N; m_pilots=1:dM_p:M;
    pilot_mask=false(N,M); pilot_mask(n_pilots,m_pilots)=true;
    data_mask=~pilot_mask; a_p=1;
    ce_mode='dftdd'; n_dd=2;

    SNRs=[3 6 12]; nTrials=40; seed_base=9000;
    Nsc_data=sum(any(data_mask,2)); Msy_data=sum(any(data_mask,1));

    fprintf('OFDM error geometry (reduced EVA 512x128, %s CE, %d trials)\n',ce_mode,nTrials);
    fprintf('%-6s | %-22s | %-22s | %-10s\n','SNR','PR_n (freq)  real|BSC','PR_m (time)  real|BSC','|H|^2@err/avg');
    for si=1:numel(SNRs)
        SNR=SNRs(si); N0=1/db2pow(SNR);
        prn_r=zeros(nTrials,1); prm_r=zeros(nTrials,1);
        prn_b=zeros(nTrials,1); prm_b=zeros(nTrials,1); hfade=zeros(nTrials,1);
        for t=1:nTrials
            rs=RandStream('Threefry','Seed',seed_base); rs.Substream=1000*si+t;
            RandStream.setGlobalStream(rs);
            nu=(nu_max*cos(2*pi*rand(1,P)-pi)).'; alpha=amp.*exp(1j*2*pi*rand(P,1));
            sym=randi(M_order,N,M)-1; Xd=qammod(sym,M_order,UnitAveragePower=true);
            Z=sqrt(N0/2)*(randn(N,M)+1j*randn(N,M));
            d=Xd; d(pilot_mask)=a_p;
            [Y,H]=ofdm_mean_H(d,deltaf,Ts,alpha,tau,nu);
            Y=Y+Z;
            H_hat=ce_ofdm(Y,d,pilot_mask,n_pilots,m_pilots,N,M,l_tau,N0,M_order,ce_mode,n_dd);
            Xs=conj(H_hat).*Y./(abs(H_hat).^2+N0);
            dh=d;
            dh(data_mask)=qammod(qamdemod(Xs(data_mask),M_order,UnitAveragePower=true),M_order,UnitAveragePower=true);
            E=false(N,M); E(data_mask)=dh(data_mask)~=d(data_mask);
            ne=nnz(E);
            if ne<2, prn_r(t)=NaN; prm_r(t)=NaN; prn_b(t)=NaN; prm_b(t)=NaN; hfade(t)=NaN; continue; end
            [prn_r(t),prm_r(t)]=prnm(E,Nsc_data,Msy_data);
            % |H|^2 at error cells vs frame average
            H2=abs(H).^2; hfade(t)=mean(H2(E))/mean(H2(data_mask));
            % BSC null: same count scattered i.i.d. over data cells
            di=find(data_mask); pe=di(randperm(numel(di),ne));
            Eb=false(N,M); Eb(pe)=true;
            [prn_b(t),prm_b(t)]=prnm(Eb,Nsc_data,Msy_data);
        end
        fprintf('%3ddB  |   %.3f  |  %.3f       |   %.3f  |  %.3f       |   %.3f\n',...
            SNR, nanmedian(prn_r),nanmedian(prn_b), nanmedian(prm_r),nanmedian(prm_b), nanmedian(hfade));
    end
    fprintf(['\nRead: PR=1 uniform over that axis, smaller=more clustered.\n' ...
             'Frequency clustering shows as PR_n(real) << PR_n(BSC); if PR_m(real)~PR_m(BSC)\n' ...
             'the time axis is unstructured (BSC-like) -> delay bias diverges from BSC, Doppler does not.\n']);
end

function [prn,prm]=prnm(E,Nsc,Msy)
    rn=sum(E,2); cm=sum(E,1);
    prn=(sum(rn)^2)/(Nsc*sum(rn.^2));
    prm=(sum(cm)^2)/(Msy*sum(cm.^2));
end

function [Y,H]=ofdm_mean_H(d,deltaf,Ts,alpha,tau,nu)
% Same mean (ICI) TF channel as ofdm_mean, but also return the per-cell gain H
% that multiplies each data symbol (diagonal term), for the deep-fade metric.
    Y=ofdm_mean(d,deltaf,Ts,alpha,tau,nu,true);
    [N,M]=size(d); n=(0:N-1).'; m=(0:M-1);
    H=zeros(N,M);
    for p=1:numel(alpha)
        H=H + alpha(p).*exp(-1j*2*pi*n*deltaf*tau(p)).*exp(1j*2*pi*m*Ts*nu(p));
    end
end

function H_hat=ce_ofdm(Y,d,pilot_mask,n_pilots,m_pilots,N,M,l_tau,N0,M_order,ce_mode,n_dd)
    H_ls=Y(n_pilots,m_pilots)./d(n_pilots,m_pilots);
    [Nq,Mq]=ndgrid(1:N,1:M);
    Afull=exp(-1j*2*pi*(0:N-1).'*(0:l_tau)/N);
    Abas=exp(-1j*2*pi*(n_pilots(:)-1)*(0:l_tau)/N);
    Hp=Afull*(Abas\H_ls);
    F_H=griddedInterpolant({(1:N).',m_pilots},Hp,'linear','nearest');
    H_hat=F_H(Nq,Mq);
    if strcmp(ce_mode,'dftdd')
        for it=1:n_dd
            Xs=conj(H_hat).*Y./(abs(H_hat).^2+N0);
            dhat=qammod(qamdemod(Xs,M_order,UnitAveragePower=true),M_order,UnitAveragePower=true);
            dhat(pilot_mask)=d(pilot_mask);
            H_hat=Afull*(Afull\(Y./dhat));
        end
    end
end

function y=nanmedian(x), y=median(x(~isnan(x))); end

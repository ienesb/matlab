function ber_awgn()
% BER_AWGN  BER vs SNR over a pure AWGN channel (single unit tap, tau=0, nu=0 ->
% H = 1, flat, no fading, no ICI) for OFDM with all channel-estimation methods
% (perfect CSI / LS+interp / DFT / DFT+DD) and CP-OTFS (MP), plus the theoretical
% 4-QAM AWGN reference BER = Q(sqrt(SNR)).  Validation case: with a flat unit
% channel the OTFS ISFFT is a unitary transform, so OFDM and OTFS must both land
% on the theory curve; the CE methods sit slightly above it by their estimation
% noise (perfect CSI is exact).
    thisdir=fileparts(mfilename('fullpath'));
    addpath(fullfile(thisdir,'utils'));
    resdir=fullfile(thisdir,'results');
    warning('off','MATLAB:nearlySingularMatrix');

    %% ---- AWGN channel: one unit-gain tap, no delay, no Doppler --------------
    N=512; M=128; deltaf=15e3; Ts=1/deltaf; M_order=4; nbits=log2(M_order);
    constellation=qammod((0:M_order-1).',M_order,UnitAveragePower=true);
    P=1; delay_taps=0; amp=1; tau=0; l_tau=0; k_nu=1;   % flat AWGN

    SNR_dbs=-6:1:12; nSNR=numel(SNR_dbs); nTrials=500; seed_base=9300;
    n_iter_mp=30; damping=0.3; n_dd=2;

    %% ---- OTFS embedded pilot + guard, OFDM overhead-matched lattice ---------
    l_p=l_tau; k_p=M/2;
    guard_rows=1:(2*l_tau+1); ce_rows=(l_p+1):(2*l_tau+1);
    data_mask_otfs=true(N,M); data_mask_otfs(guard_rows,:)=false;
    otfs_nondata=numel(guard_rows)*M;
    [dN_p,dM_p,n_pilot_cells]=match_pilot_spacing(N,M,otfs_nondata,l_tau,k_nu);
    n_pilots=1:dN_p:N; m_pilots=1:dM_p:M;
    pilot_mask=false(N,M); pilot_mask(n_pilots,m_pilots)=true; data_mask=~pilot_mask;
    E_pil=n_pilot_cells; a_p=1; x_p=sqrt(E_pil);

    configs={ 'OFDM (perfect CSI)','ofdm','perfect'
              'OFDM (LS+interp)',  'ofdm','interp'
              'OFDM (DFT CE)',     'ofdm','dft'
              'OFDM (DFT+DD CE)',  'ofdm','dftdd'
              'CP-OTFS (MP)',      'otfs','' };
    nC=size(configs,1); ber=nan(nSNR,nC);

    fprintf('AWGN BER (flat unit channel), N=%d M=%d, %d trials\n',N,M,nTrials);
    tAll=tic;
    for ci=1:nC
        modn=configs{ci,2}; ce_mode=configs{ci,3};
        for si=1:nSNR
            N0=1/db2pow(SNR_dbs(si)); be_t=zeros(1,nTrials);
            parfor t=1:nTrials
                warning('off','MATLAB:nearlySingularMatrix');
                rs=RandStream('Threefry','Seed',seed_base); rs.Substream=1000*si+t;
                RandStream.setGlobalStream(rs);
                nu=zeros(P,1); alpha=amp.*exp(1j*2*pi*rand(P,1));   % |H|=1, random phase
                sym=randi(M_order,N,M)-1; Xd=qammod(sym,M_order,UnitAveragePower=true);
                Z=sqrt(N0/2)*(randn(N,M)+1j*randn(N,M));
                if strcmp(modn,'ofdm')
                    d=Xd; d(pilot_mask)=a_p;
                    Y=ofdm_mean(d,deltaf,Ts,alpha,tau,nu,true)+Z;
                    if strcmp(ce_mode,'perfect')
                        H_hat=true_diag_gain(N,M,deltaf,Ts,alpha,tau,nu);
                    else
                        H_hat=ce_ofdm(Y,d,pilot_mask,n_pilots,m_pilots,N,M,l_tau,N0,M_order,ce_mode,n_dd);
                    end
                    Xsoft=conj(H_hat).*Y./(abs(H_hat).^2+N0);
                    sh=qamdemod(Xsoft(data_mask),M_order,UnitAveragePower=true);
                    be_t(t)=mean(int2bit(sym(data_mask),nbits)~=int2bit(sh,nbits),'all');
                else
                    d=zeros(N,M); d(data_mask_otfs)=Xd(data_mask_otfs); d(l_p+1,k_p+1)=x_p;
                    Gtrue=build_effective_dd_eva(alpha,delay_taps,nu,N,M,deltaf,Inf);
                    Y=apply_dd_channel_eva(d,Gtrue)+Z;
                    h_ce=zeros(N,M); h_ce(ce_rows,:)=Y(ce_rows,:)/x_p;
                    h_ce(abs(h_ce)<3*sqrt(N0)/x_p)=0;
                    [rc,cc,gc]=find(h_ce);
                    if isempty(gc)
                        be_t(t)=0.5;
                    else
                        dl_c=mod(rc-(l_p+1),N); dk_c=mod(cc-(k_p+1),M);
                        known=zeros(N,M); known(l_p+1,k_p+1)=x_p;
                        xhat=mp_detector(Y,dl_c,dk_c,gc,N0,constellation,n_iter_mp,data_mask_otfs,known,damping);
                        sh=qamdemod(xhat(data_mask_otfs),M_order,UnitAveragePower=true);
                        be_t(t)=mean(int2bit(sym(data_mask_otfs),nbits)~=int2bit(sh,nbits),'all');
                    end
                end
            end
            ber(si,ci)=mean(be_t);
            if mod(si,6)==1||si==nSNR
                fprintf('  %-18s SNR=%3d dB | BER=%.3e (%.0fs)\n',configs{ci,1},SNR_dbs(si),ber(si,ci),toc(tAll));
            end
        end
    end

    ber_theory=0.5*erfc(sqrt(db2pow(SNR_dbs))/sqrt(2));    % 4-QAM Gray AWGN: Q(sqrt(SNR))

    outpath=fullfile(resdir,'ber_awgn_results.mat');
    save(outpath,'SNR_dbs','ber','ber_theory','configs','N','M','nTrials');

    %% ---- plot --------------------------------------------------------------
    cols=[0 0 0; 0.95 0.55 0.35; 0.85 0.33 0.10; 0.93 0.69 0.13; 0 0.45 0.74];
    mrk={'o','o','o','o','d'};
    hf=figure('Position',[200 160 860 580],'Color','w'); hold on; grid on; box on;
    semilogy(SNR_dbs,ber_theory,'k--','LineWidth',1.6);        % theory reference
    lg={'theory Q(\surdSNR)'};
    for ci=1:nC
        y=ber(:,ci); y(y==0)=NaN;
        semilogy(SNR_dbs,y,['-' mrk{ci}],'Color',cols(ci,:),'LineWidth',1.9,'MarkerSize',5);
        lg=[lg,configs(ci,1)]; %#ok<AGROW>
    end
    set(gca,'YScale','log'); xlabel('SNR = E_s/N_0 [dB]'); ylabel('BER');
    title(sprintf(['AWGN (flat unit channel): OFDM (all CE) & CP-OTFS vs theory (%d-QAM, %d trials)\n' ...
                   'all curves should collapse onto Q(\\surdSNR); CE noise lifts the estimated arms'],M_order,nTrials));
    legend(lg,'Location','southwest','FontSize',9); set(gca,'FontSize',11);
    if exist('theme','file'), theme(hf,"light"); end
    png=fullfile(resdir,'ber_awgn.png');
    exportgraphics(hf,png,'Resolution',150,'BackgroundColor','white');
    fprintf('\nSaved %s and %s\nDone in %.0fs\n',png,outpath,toc(tAll));
end

function H=true_diag_gain(N,M,deltaf,Ts,alpha,tau,nu)
    ns=(0:N-1).'; ms=(0:M-1); H=zeros(N,M); T=1/deltaf;
    for p=1:numel(alpha)
        I0=get_ici_coeffs(0,nu(p),T);
        H=H+alpha(p)*I0*(exp(-1j*2*pi*deltaf*ns*tau(p))*exp(1j*2*pi*Ts*ms*nu(p)));
    end
end

function H_hat=ce_ofdm(Y,d,pilot_mask,n_pilots,m_pilots,N,M,l_tau,N0,M_order,ce_mode,n_dd)
    H_ls=Y(n_pilots,m_pilots)./d(n_pilots,m_pilots); [Nq,Mq]=ndgrid(1:N,1:M);
    Afull=exp(-1j*2*pi*(0:N-1).'*(0:l_tau)/N);
    switch ce_mode
        case 'interp'
            F_H=griddedInterpolant({n_pilots,m_pilots},H_ls,'linear','nearest'); H_hat=F_H(Nq,Mq); return
        case {'dft','dftdd'}
            Abas=exp(-1j*2*pi*(n_pilots(:)-1)*(0:l_tau)/N); Hp=Afull*(Abas\H_ls);
            F_H=griddedInterpolant({(1:N).',m_pilots},Hp,'linear','nearest'); H_hat=F_H(Nq,Mq);
    end
    if strcmp(ce_mode,'dftdd')
        for it=1:n_dd
            Xs=conj(H_hat).*Y./(abs(H_hat).^2+N0);
            dhat=qammod(qamdemod(Xs,M_order,UnitAveragePower=true),M_order,UnitAveragePower=true);
            dhat(pilot_mask)=d(pilot_mask); H_hat=Afull*(Afull\(Y./dhat));
        end
    end
end

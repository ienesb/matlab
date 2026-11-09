function ber_rayleigh_flat()
% BER_RAYLEIGH_FLAT  BER vs SNR over a SINGLE-PATH flat Rayleigh channel (one tap,
% tau=0, nu=0 -> H = alpha, alpha ~ CN(0,1) per frame; genuine Rayleigh magnitude,
% no frequency selectivity, no Doppler, no ICI) for OFDM (perfect CSI) and CP-OTFS
% (MP), vs the theoretical flat-Rayleigh 4-QAM BER 0.5(1-sqrt(g/(1+g))), g=Eb/N0.
%
% A single flat fade gives NO diversity to spread over, so OTFS's ISFFT (unitary)
% leaves the per-symbol SNR |alpha|^2*SNR unchanged -> OFDM and OTFS must coincide,
% both on the order-1 Rayleigh curve.  (Contrast the reduced-EVA 3-tap channel,
% where OTFS harvests delay diversity and pulls away from OFDM.)
    thisdir=fileparts(mfilename('fullpath'));
    addpath(fullfile(thisdir,'utils'));
    resdir=fullfile(thisdir,'results'); figdir=fullfile(thisdir,'figures');
    if ~exist(figdir,'dir'), mkdir(figdir); end
    warning('off','MATLAB:nearlySingularMatrix');

    N=512; M=128; deltaf=15e3; Ts=1/deltaf; M_order=4; nbits=log2(M_order);
    constellation=qammod((0:M_order-1).',M_order,UnitAveragePower=true);
    P=1; delay_taps=0; tau=0; l_tau=0; k_nu=1;          % single flat tap

    SNR_dbs=0:2:34; nSNR=numel(SNR_dbs); nTrials=1000; seed_base=9500;
    n_iter_mp=30; damping=0.3;

    % OTFS embedded pilot + guard; OFDM overhead-matched lattice
    l_p=l_tau; k_p=M/2;
    guard_rows=1:(2*l_tau+1); ce_rows=(l_p+1):(2*l_tau+1);
    data_mask_otfs=true(N,M); data_mask_otfs(guard_rows,:)=false;
    otfs_nondata=numel(guard_rows)*M;
    [dN_p,dM_p,n_pilot_cells]=match_pilot_spacing(N,M,otfs_nondata,l_tau,k_nu);
    n_pilots=1:dN_p:N; m_pilots=1:dM_p:M;
    pilot_mask=false(N,M); pilot_mask(n_pilots,m_pilots)=true; data_mask=~pilot_mask;
    a_p=1; x_p=sqrt(n_pilot_cells);

    configs={ 'OFDM (perfect CSI)','ofdm'; 'CP-OTFS (MP)','otfs' };
    nC=size(configs,1); ber=nan(nSNR,nC);

    fprintf('Single-path flat Rayleigh BER, N=%d M=%d, %d trials\n',N,M,nTrials);
    tAll=tic;
    for ci=1:nC
        modn=configs{ci,2};
        for si=1:nSNR
            N0=1/db2pow(SNR_dbs(si)); be_t=zeros(1,nTrials);
            parfor t=1:nTrials
                warning('off','MATLAB:nearlySingularMatrix');
                rs=RandStream('Threefry','Seed',seed_base); rs.Substream=1000*si+t;
                RandStream.setGlobalStream(rs);
                alpha=(randn+1j*randn)/sqrt(2);                 % CN(0,1): flat Rayleigh
                nu=0;
                sym=randi(M_order,N,M)-1; Xd=qammod(sym,M_order,UnitAveragePower=true);
                Z=sqrt(N0/2)*(randn(N,M)+1j*randn(N,M));
                if strcmp(modn,'ofdm')
                    d=Xd; d(pilot_mask)=a_p;
                    Y=alpha*d+Z;                                % flat channel (H=alpha everywhere)
                    H_hat=alpha*ones(N,M);                      % perfect CSI
                    Xs=conj(H_hat).*Y./(abs(H_hat).^2+N0);
                    sh=qamdemod(Xs(data_mask),M_order,UnitAveragePower=true);
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
            if mod(si,4)==1||si==nSNR
                fprintf('  %-18s SNR=%3d dB | BER=%.3e (%.0fs)\n',configs{ci,1},SNR_dbs(si),ber(si,ci),toc(tAll));
            end
        end
    end

    gb=db2pow(SNR_dbs)/2;                                 % Eb/N0 (QPSK: Es/N0 = SNR, 2 bits)
    ber_theory=0.5*(1-sqrt(gb./(1+gb)));                  % flat-Rayleigh 4-QAM Gray

    save(fullfile(resdir,'ber_rayleigh_flat_results.mat'),'SNR_dbs','ber','ber_theory','configs','N','M','nTrials');

    orange=[0.85 0.33 0.10]; blue=[0 0.45 0.74];
    hf=figure('Name','BER: single-path flat Rayleigh','Position',[200 160 840 560],'Color','w');
    hold on; grid on; box on;
    semilogy(SNR_dbs,ber_theory,'k--','LineWidth',1.6);
    semilogy(SNR_dbs,nz(ber(:,1)),'-o','Color',orange,'LineWidth',2,'MarkerSize',6,'MarkerFaceColor',orange);
    semilogy(SNR_dbs,nz(ber(:,2)),'-d','Color',blue,  'LineWidth',2,'MarkerSize',6,'MarkerFaceColor',blue);
    set(gca,'YScale','log'); xlabel('SNR = E_s/N_0 [dB]'); ylabel('BER');
    title(sprintf('Single-path FLAT Rayleigh (no Doppler, no ICI): OFDM vs CP-OTFS (%d-QAM, %d trials)',M_order,nTrials));
    legend({'flat-Rayleigh theory','OFDM (perfect CSI)','CP-OTFS (MP)'},'Location','southwest','FontSize',10);
    set(gca,'FontSize',11); if exist('theme','file'), theme(hf,"light"); end
    savefig(hf,fullfile(figdir,'ber_rayleigh_flat.fig'));
    exportgraphics(hf,fullfile(resdir,'ber_rayleigh_flat.png'),'Resolution',150,'BackgroundColor','white');
    fprintf('Saved %s and .png\n',fullfile(figdir,'ber_rayleigh_flat.fig'));

    fprintf('\n%4s | %10s %10s %10s\n','SNR','OFDM','CP-OTFS','theory');
    for si=1:nSNR
        fprintf('%4d | %10.3e %10.3e %10.3e\n',SNR_dbs(si),ber(si,1),ber(si,2),ber_theory(si));
    end
end

function y=nz(y), y(y==0)=NaN; end

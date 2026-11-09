function ber_confirm_config(guard_extra, Emax_dd, nTrials, SNR_dbs)
% BER_CONFIRM_CONFIG(guard_extra, Emax_dd, nTrials, SNR_dbs)
% BER vs SNR for the confirmation setup: reduced EVA, N=512, M=128, deltaf=30 kHz
% (B=15.4 MHz), Doppler ~444 Hz + ICI.  Same realistic receivers as the LB sweep
% (OFDM: dftdd CE + 1-tap MMSE ; CP-OTFS: threshold CE + MP).
%   guard_extra : OTFS delay guard G = l_tau+1+guard_extra.  Default l_tau -> full
%                 guard G=2*l_tau+1=13 (OTFS keeps its BER cliff).  =2 -> G=9 (shrunk;
%                 OTFS CE degrades -> floors).
%   Emax_dd     : truncate the DD Doppler spread to |q|<=Emax in build_effective_dd_eva
%                 (default 8; Inf = exact but the MP gets very slow at high SNR).
% Per-SNR checkpoint written to CKPT_<outfile> so a stop never loses data.
% Companion to plot_mcrb_confirm: full guard -> OTFS lower BER (cliff) while OFDM
% has the lower delay+Doppler LB at low SNR ("BER != sensing bound").
    thisdir=fileparts(mfilename('fullpath')); addpath(fullfile(thisdir,'utils'));
    resdir=fullfile(thisdir,'results'); figdir=fullfile(thisdir,'figures');
    if ~exist(figdir,'dir'), mkdir(figdir); end
    warning('off','MATLAB:nearlySingularMatrix');

    N=512; M=128; deltaf=30e3; Ts=1/deltaf; fc=4e9; c=3e8; v=120/3.6; nu_max=v*fc/c;
    M_order=4; nbits=log2(M_order); constellation=qammod((0:M_order-1).',M_order,UnitAveragePower=true);
    eva_delays_ns=[0 150 370]; eva_pow_db=[0 -1.4 -0.6]; P=numel(eva_delays_ns);
    delay_taps=round(eva_delays_ns*1e-9*N*deltaf); assert(numel(unique(delay_taps))==P);   % [0 2 6]
    tau=(delay_taps/(N*deltaf)).'; amp=(10.^(eva_pow_db/20)).'; amp=amp/norm(amp);
    l_tau=max(delay_taps); k_nu=max(1,ceil(nu_max*M/deltaf));
    n_iter_mp=30; damping=0.3; ce_mode='dftdd'; n_dd=2;

    if nargin<1||isempty(guard_extra), guard_extra=l_tau; end   % default: full guard G=2*l_tau+1
    if nargin<2||isempty(Emax_dd), Emax_dd=8; end
    guard_extra=max(0,min(round(guard_extra),l_tau));
    G=l_tau+1+guard_extra;                                % l_tau -> G=13 (full); 2 -> G=9 (shrunk)
    l_p=guard_extra; k_p=M/2; guard_rows=1:G; ce_rows=(l_p+1):G;
    data_mask_otfs=true(N,M); data_mask_otfs(guard_rows,:)=false;
    otfs_nondata=numel(guard_rows)*M;
    [dN_p,dM_p,n_pilot_cells]=match_pilot_spacing(N,M,otfs_nondata,l_tau,k_nu);
    n_pilots=1:dN_p:N; m_pilots=1:dM_p:M;
    pilot_mask=false(N,M); pilot_mask(n_pilots,m_pilots)=true; data_mask_ofdm=~pilot_mask;
    a_p=1; x_p=sqrt(n_pilot_cells);

    if nargin<3||isempty(nTrials), nTrials=100; end
    if nargin<4||isempty(SNR_dbs), SNR_dbs=0:1:25; end
    nSNR=numel(SNR_dbs); seed_base=8300;
    outfile=sprintf('ber_confirm_config_g%d_results.mat',G);   % new name per guard size
    configs={'OFDM (dftdd CE)','ofdm';'CP-OTFS (MP)','otfs'};
    ber=nan(nSNR,2);
    fprintf('BER, confirm config: N=%d M=%d df=%gk B=%.1fMHz, Doppler+ICI, OTFS G=%d, Emax=%g, %d trials\n', ...
            N,M,deltaf/1e3,N*deltaf/1e6,G,Emax_dd,nTrials);
    tAll=tic;
    for ci=1:2
        modn=configs{ci,2};
        for si=1:nSNR
            N0=1/db2pow(SNR_dbs(si)); bb=zeros(1,nTrials);
            parfor t=1:nTrials
                warning('off','MATLAB:nearlySingularMatrix');
                rs=RandStream('Threefry','Seed',seed_base); rs.Substream=1000*si+t; RandStream.setGlobalStream(rs);
                nu=(nu_max*cos(2*pi*rand(1,P)-pi)).'; alpha=amp.*exp(1j*2*pi*rand(P,1));
                sym=randi(M_order,N,M)-1; Xd=qammod(sym,M_order,UnitAveragePower=true);
                Z=sqrt(N0/2)*(randn(N,M)+1j*randn(N,M));
                if strcmp(modn,'ofdm')
                    d=Xd; d(pilot_mask)=a_p;
                    Y=ofdm_mean(d,deltaf,Ts,alpha,tau,nu,true)+Z;
                    H_hat=ce_ofdm(Y,d,pilot_mask,n_pilots,m_pilots,N,M,l_tau,N0,M_order,ce_mode,n_dd);
                    Xs=conj(H_hat).*Y./(abs(H_hat).^2+N0);
                    sh=qamdemod(Xs(data_mask_ofdm),M_order,UnitAveragePower=true);
                    bb(t)=mean(int2bit(sym(data_mask_ofdm),nbits)~=int2bit(sh,nbits),'all');
                else
                    d=zeros(N,M); d(data_mask_otfs)=Xd(data_mask_otfs); d(l_p+1,k_p+1)=x_p;
                    Gtrue=build_effective_dd_eva(alpha,delay_taps,nu,N,M,deltaf,Emax_dd);
                    Y=apply_dd_channel_eva(d,Gtrue)+Z;
                    h_ce=zeros(N,M); h_ce(ce_rows,:)=Y(ce_rows,:)/x_p; h_ce(abs(h_ce)<3*sqrt(N0)/x_p)=0;
                    [rc,cc,gc]=find(h_ce);
                    if isempty(gc)
                        bb(t)=0.5;
                    else
                        dl_c=mod(rc-(l_p+1),N); dk_c=mod(cc-(k_p+1),M); known=zeros(N,M); known(l_p+1,k_p+1)=x_p;
                        xhat=mp_detector(Y,dl_c,dk_c,gc,N0,constellation,n_iter_mp,data_mask_otfs,known,damping);
                        sh=qamdemod(xhat(data_mask_otfs),M_order,UnitAveragePower=true);
                        bb(t)=mean(int2bit(sym(data_mask_otfs),nbits)~=int2bit(sh,nbits),'all');
                    end
                end
            end
            ber(si,ci)=mean(bb);
            fprintf('  %-16s SNR=%2d | BER=%.3e (%.0fs)\n',configs{ci,1},SNR_dbs(si),ber(si,ci),toc(tAll));
            save(fullfile(resdir,['CKPT_' outfile]),'SNR_dbs','ber','configs','N','M','deltaf','G','nTrials','ci','si');  % checkpoint
        end
    end

    save(fullfile(resdir,outfile),'SNR_dbs','ber','configs','N','M','deltaf','G','Emax_dd','nTrials');
    stem=sprintf('ber_confirm_config_g%d',G);
    orange=[0.85 0.33 0.10]; blue=[0 0.45 0.74];
    hf=figure('Position',[200 160 840 560],'Color','w'); hold on; grid on; box on;
    y=ber(:,1); y(y==0)=NaN; semilogy(SNR_dbs,y,'-o','Color',orange,'LineWidth',2,'MarkerSize',6,'MarkerFaceColor',orange);
    y=ber(:,2); y(y==0)=NaN; semilogy(SNR_dbs,y,'-d','Color',blue,'LineWidth',2,'MarkerSize',6,'MarkerFaceColor',blue);
    set(gca,'YScale','log'); xlabel('SNR [dB]'); ylabel('BER');
    title(sprintf(['BER, "OFDM wins both LB" config: reduced EVA, B=15.4 MHz, Doppler+ICI, OTFS G=%d (%d trials)\n' ...
                   'OTFS has the lower BER (its cliff); OFDM has the lower delay+Doppler LB at 0-4 dB'],G,nTrials));
    legend({'OFDM (dftdd CE)','CP-OTFS (MP)'},'Location','southwest','FontSize',10); set(gca,'FontSize',11);
    if exist('theme','file'), theme(hf,"light"); end
    savefig(hf,fullfile(figdir,[stem '.fig']));
    exportgraphics(hf,fullfile(resdir,[stem '.png']),'Resolution',150,'BackgroundColor','white');
    fprintf('Saved %s and .png\nDone %.0fs\n',fullfile(figdir,[stem '.fig']),toc(tAll));

    fprintf('\n%4s | %10s %10s\n','SNR','OFDM','CP-OTFS');
    for si=1:nSNR, if mod(SNR_dbs(si),3)==0, fprintf('%4d | %10.3e %10.3e\n',SNR_dbs(si),ber(si,1),ber(si,2)); end, end
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

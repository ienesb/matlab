function ber_rayleigh_eva(fading, doppler_on)
% BER_RAYLEIGH_EVA  BER vs SNR for the reduced-EVA 3-tap delay profile.
% OFDM (perfect CSI) vs CP-OTFS (MP), vs the flat-Rayleigh (order-1) reference.
% Two toggles:
%   fading      : 'rayleigh' (default) alpha_p ~ CN(0,power_p)  -- true Rayleigh
%                 'constmag'           alpha_p = sqrt(power_p)*e^{j*phase} (usual EVA)
%   doppler_on  : false (default) NO Doppler -> NO ICI (nu=0).
%                 true            nominal EVA Doppler (~444 Hz) per tap -> ICI on
%                 (OFDM channel via the ICI kernel; per-cell perfect CSI still uses
%                 the diagonal gain, so residual ICI floors OFDM's BER).
    if nargin<1||isempty(fading), fading='rayleigh'; end
    if nargin<2||isempty(doppler_on), doppler_on=false; end
    assert(any(strcmp(fading,{'rayleigh','constmag'})),'fading must be rayleigh|constmag');
%
% Expect: OFDM stays on the order-1 Rayleigh curve (each subcarrier is a sum of
% Gaussians -> Rayleigh -> diversity 1, independent of the tap count), while
% CP-OTFS harvests the 3-path DELAY diversity and drops below it (steeper slope).
    thisdir=fileparts(mfilename('fullpath'));
    addpath(fullfile(thisdir,'utils'));
    resdir=fullfile(thisdir,'results'); figdir=fullfile(thisdir,'figures');
    if ~exist(figdir,'dir'), mkdir(figdir); end
    warning('off','MATLAB:nearlySingularMatrix');

    N=512; M=128; deltaf=15e3; Ts=1/deltaf; M_order=4; nbits=log2(M_order);
    constellation=qammod((0:M_order-1).',M_order,UnitAveragePower=true);

    eva_delays_ns=[0 150 370]; eva_pow_db=[0 -1.4 -0.6]; P=numel(eva_delays_ns);
    delay_taps=round(eva_delays_ns*1e-9*N*deltaf);            % [0 1 3]
    tau=(delay_taps/(N*deltaf)).';
    plin=10.^(eva_pow_db/20); amp=(plin/norm(plin)).';       % per-tap std (sum|amp|^2=1)
    l_tau=max(delay_taps);
    fc=4e9; c=3e8; v=120/3.6;  nu_max=(v*fc/c)*doppler_on;   % ~444 Hz if doppler_on, else 0
    k_nu=max(1,ceil(nu_max*M/deltaf));

    SNR_dbs=0:2:34; nSNR=numel(SNR_dbs); nTrials=1000; seed_base=9600;
    n_iter_mp=30; damping=0.3;

    l_p=l_tau; k_p=M/2;
    guard_rows=1:(2*l_tau+1); ce_rows=(l_p+1):(2*l_tau+1);
    data_mask_otfs=true(N,M); data_mask_otfs(guard_rows,:)=false;
    otfs_nondata=numel(guard_rows)*M;
    [dN_p,dM_p,n_pilot_cells]=match_pilot_spacing(N,M,otfs_nondata,l_tau,k_nu);
    n_pilots=1:dN_p:N; m_pilots=1:dM_p:M;
    pilot_mask=false(N,M); pilot_mask(n_pilots,m_pilots)=true; data_mask=~pilot_mask;
    a_p=1; x_p=sqrt(n_pilot_cells);
    ns=(0:N-1).'; ms=(0:M-1);

    configs={ 'OFDM (perfect CSI)','ofdm'; 'CP-OTFS (MP)','otfs' };
    nC=size(configs,1); ber=nan(nSNR,nC);

    dtxt=tern(doppler_on,'Doppler ~444 Hz + ICI','no Doppler / no ICI');
    fprintf('Reduced-EVA 3-tap BER, fading=%s, %s, N=%d M=%d, %d trials\n',fading,dtxt,N,M,nTrials);
    tAll=tic;
    for ci=1:nC
        modn=configs{ci,2};
        for si=1:nSNR
            N0=1/db2pow(SNR_dbs(si)); be_t=zeros(1,nTrials);
            parfor t=1:nTrials
                warning('off','MATLAB:nearlySingularMatrix');
                rs=RandStream('Threefry','Seed',seed_base); rs.Substream=1000*si+t;
                RandStream.setGlobalStream(rs);
                if strcmp(fading,'rayleigh')
                    alpha=(randn(P,1)+1j*randn(P,1))/sqrt(2).*amp;   % CN(0,power_p): Rayleigh
                else
                    alpha=amp.*exp(1j*2*pi*rand(P,1));               % fixed |amp|, random phase
                end
                if doppler_on, nu=(nu_max*cos(2*pi*rand(P,1)-pi)); else, nu=zeros(P,1); end
                sym=randi(M_order,N,M)-1; Xd=qammod(sym,M_order,UnitAveragePower=true);
                Z=sqrt(N0/2)*(randn(N,M)+1j*randn(N,M));
                if strcmp(modn,'ofdm')
                    d=Xd; d(pilot_mask)=a_p;
                    Y=ofdm_mean(d,deltaf,Ts,alpha,tau,nu,doppler_on)+Z;  % doppler_on -> full ICI kernel
                    H_hat=zeros(N,M);                                    % perfect CSI (diagonal gain)
                    for p=1:P
                        I0=get_ici_coeffs(0,nu(p),Ts);                   % ICI self-coeff (=1 if nu=0)
                        H_hat=H_hat+alpha(p)*I0*(exp(-1j*2*pi*deltaf*ns*tau(p))*exp(1j*2*pi*Ts*ms*nu(p)));
                    end
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

    gb=db2pow(SNR_dbs)/2;  ber_theory=0.5*(1-sqrt(gb./(1+gb)));   % flat-Rayleigh order-1 reference

    stem=['ber_' fading '_eva' tern(doppler_on,'_dopp','')];
    ftxt=tern(strcmp(fading,'rayleigh'),'Rayleigh','const-|amp|, random phase');
    save(fullfile(resdir,[stem '_results.mat']),'SNR_dbs','ber','ber_theory','configs','N','M','nTrials','eva_delays_ns','fading','doppler_on','nu_max');

    orange=[0.85 0.33 0.10]; blue=[0 0.45 0.74];
    hf=figure('Name',['BER: 3-tap EVA ' fading],'Position',[200 160 840 560],'Color','w');
    hold on; grid on; box on;
    semilogy(SNR_dbs,ber_theory,'k--','LineWidth',1.6);
    semilogy(SNR_dbs,nz(ber(:,1)),'-o','Color',orange,'LineWidth',2,'MarkerSize',6,'MarkerFaceColor',orange);
    semilogy(SNR_dbs,nz(ber(:,2)),'-d','Color',blue,  'LineWidth',2,'MarkerSize',6,'MarkerFaceColor',blue);
    set(gca,'YScale','log'); xlabel('SNR = E_s/N_0 [dB]'); ylabel('BER');
    title(sprintf('Reduced-EVA 3-tap, %s, %s: OFDM vs CP-OTFS (%d-QAM, %d trials)',ftxt,dtxt,M_order,nTrials));
    legend({'order-1 (flat-Rayleigh) reference','OFDM (perfect CSI)','CP-OTFS (MP)'},'Location','southwest','FontSize',10);
    set(gca,'FontSize',11); if exist('theme','file'), theme(hf,"light"); end
    savefig(hf,fullfile(figdir,[stem '.fig']));
    exportgraphics(hf,fullfile(resdir,[stem '.png']),'Resolution',150,'BackgroundColor','white');
    fprintf('Saved %s and .png\n',fullfile(figdir,[stem '.fig']));

    fprintf('\n%4s | %10s %10s %10s\n','SNR','OFDM','CP-OTFS','order1ref');
    for si=1:nSNR
        fprintf('%4d | %10.3e %10.3e %10.3e\n',SNR_dbs(si),ber(si,1),ber(si,2),ber_theory(si));
    end
end

function y=nz(y), y(y==0)=NaN; end
function s=tern(c,a,b), if c, s=a; else, s=b; end, end

function ber_ofdm_ce_vs_perfect(SNR_dbs, nu_scale)
% ber_ofdm_ce_vs_perfect(SNR_dbs, nu_scale)
%   SNR_dbs  : optional SNR vector; those points are computed and MERGED into the
%              matching results file (OFDM only). No argument -> -10:1:35.
%   nu_scale : Doppler scale (default 1 = nominal EVA ~444 Hz). nu_scale=0 -> no
%              Doppler (ICI kernel collapses to delta), written to a separate
%              *_nodopp results file so the nominal run is preserved.  Use this to
%              test whether the perfect-CSI floor is caused by ICI.
% BER_OFDM_CE_VS_PERFECT  OFDM BER vs SNR on the reduced-EVA channel with the
% three practical channel-estimation methods vs PERFECT channel knowledge at the
% receiver.  Isolates CE-estimation error from the equalizer/diversity limit:
% the perfect-CSI curve uses the TRUE per-cell diagonal gain in the same 1-tap
% MMSE, so any gap above it is CE error, and its own floor is the residual-ICI +
% order-1 diversity wall that no CE can beat.
%
%   OFDM (perfect CSI) : H_hat = true diagonal channel gain (genie)
%   OFDM (LS+interp)   : LS at pilots -> linear interpolation
%   OFDM (DFT CE)      : LS delay-tap fit across the pilot comb
%   OFDM (DFT+DD CE)   : + decision-directed tap refit from all subcarriers
%
% The existing CP-OTFS (MP) curve for reduced EVA is loaded from
% ber_ce_channels_results.mat and overlaid unchanged (not recomputed).
    thisdir = fileparts(mfilename('fullpath'));
    addpath(fullfile(thisdir,'utils'));
    resdir = fullfile(thisdir,'results');
    warning('off','MATLAB:nearlySingularMatrix');

    %% ---- reduced EVA, nominal (vehicular) Doppler --------------------------
    if nargin<2 || isempty(nu_scale), nu_scale=1; end
    N=512; M=128; deltaf=15e3; Ts=1/deltaf; fc=4e9; c=3e8; v=120/3.6;
    nu_max=v*fc/c * nu_scale;                         % ~444 Hz at nu_scale=1; 0 at nu_scale=0
    dtag = ''; if nu_scale==0, dtag='_nodopp'; elseif nu_scale~=1, dtag=sprintf('_nu%gx',nu_scale); end
    delay_taps=round([0 150 370]*1e-9*N*deltaf);
    amp=(10.^([0 -1.4 -0.6]/20)).'; amp=amp/norm(amp);
    tau=(delay_taps/(N*deltaf)).'; P=numel(delay_taps); l_tau=max(delay_taps);
    k_nu=max(1,ceil(nu_max*M/deltaf));
    M_order=4; nbits=log2(M_order);

    if nargin<1 || isempty(SNR_dbs), SNR_dbs=-10:1:35; end
    SNR_dbs=SNR_dbs(:).'; nSNR=numel(SNR_dbs); nTrials=500; seed_base=9100;
    n_dd=2;

    % ---- overhead-matched, data-referenced pilot lattice (as in mcrb_ce_ici) --
    guard_rows=1:(2*l_tau+1); otfs_nondata=numel(guard_rows)*M;
    [dN_p,dM_p,n_pilot_cells]=match_pilot_spacing(N,M,otfs_nondata,l_tau,k_nu);
    n_pilots=1:dN_p:N; m_pilots=1:dM_p:M;
    pilot_mask=false(N,M); pilot_mask(n_pilots,m_pilots)=true; data_mask=~pilot_mask;
    a_p=1;                                           % E_pil = n_pilot_cells, boost 0 dB -> unit pilots

    configs={ 'OFDM (perfect CSI)','perfect'
              'OFDM (LS+interp)',  'interp'
              'OFDM (DFT CE)',     'dft'
              'OFDM (DFT+DD CE)',  'dftdd' };
    nC=size(configs,1); ber=nan(nSNR,nC);

    fprintf('OFDM CE vs perfect CSI, reduced EVA (N=%d M=%d P=%d nu_max=%.0f Hz), %d trials\n',...
            N,M,P,nu_max,nTrials);
    fprintf('pilot lattice (%d,%d), %d pilot cells (%.2f%%); delay-Nyquist l_tau=%d < N/dN_p=%.1f\n',...
            dN_p,dM_p,n_pilot_cells,100*n_pilot_cells/(N*M),l_tau,N/dN_p);
    tAll=tic;
    for ci=1:nC
        ce_mode=configs{ci,2};
        for si=1:nSNR
            N0=1/db2pow(SNR_dbs(si)); be_t=zeros(1,nTrials);
            parfor t=1:nTrials
                warning('off','MATLAB:nearlySingularMatrix');
                rs=RandStream('Threefry','Seed',seed_base); rs.Substream=1000*si+t;
                RandStream.setGlobalStream(rs);
                nu=(nu_max*cos(2*pi*rand(1,P)-pi)).'; alpha=amp.*exp(1j*2*pi*rand(P,1));
                sym=randi(M_order,N,M)-1; Xd=qammod(sym,M_order,UnitAveragePower=true);
                Z=sqrt(N0/2)*(randn(N,M)+1j*randn(N,M));
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
            end
            ber(si,ci)=mean(be_t);
            if mod(si,12)==1||si==nSNR
                fprintf('  %-18s SNR=%3d dB | BER=%.3e (%.0fs)\n',configs{ci,1},SNR_dbs(si),ber(si,ci),toc(tAll));
            end
        end
    end

    %% ---- overlay existing CP-OTFS (MP) reduced-EVA curve (not recomputed) ---
    % Match the Doppler regime: nominal -> ber_ce_channels; no-Doppler -> _nodopp.
    otfs_ber=[]; otfs_snr=[];
    ref=fullfile(resdir,['ber_ce_channels' dtag '_results.mat']);
    if isfile(ref)
        R=load(ref); ei=find(strcmpi(cellstr(R.channel_set),'eva'),1);
        oi=find(strcmp(R.configs(:,1),'CP-OTFS (MP)'),1);
        if ~isempty(ei)&&~isempty(oi), otfs_ber=R.ber(:,oi,ei); otfs_snr=R.SNR_dbs(:); end
    end

    %% ---- merge with existing results (OFDM only): union of SNR points -------
    outpath=fullfile(resdir,['ber_ofdm_ce_vs_perfect' dtag '_results.mat']);
    if isfile(outpath)
        E=load(outpath);
        if isequal(E.configs(:,1),configs(:,1))
            allSNR=union(E.SNR_dbs(:),SNR_dbs(:)).';           % sorted, unique
            merged=nan(numel(allSNR),nC);
            [~,loc]=ismember(E.SNR_dbs(:),allSNR); merged(loc,:)=E.ber;   % keep old
            [~,loc]=ismember(SNR_dbs(:),allSNR);   merged(loc,:)=ber;     % new wins on overlap
            SNR_dbs=allSNR; ber=merged; nSNR=numel(allSNR);
            fprintf('Merged with existing: SNR now %d..%d dB (%d pts)\n',allSNR(1),allSNR(end),nSNR);
        end
    end
    save(outpath,'SNR_dbs','ber','configs','N','M','P','nu_max','nTrials','otfs_ber','otfs_snr');

    %% ---- plot --------------------------------------------------------------
    cols=[0 0 0; 0.95 0.55 0.35; 0.85 0.33 0.10; 0.93 0.69 0.13];  % perfect=black
    hf=figure('Position',[200 160 860 580],'Color','w'); hold on; grid on; box on;
    lg={};
    for ci=1:nC
        y=ber(:,ci); y(y==0)=NaN;
        semilogy(SNR_dbs,y,'-o','Color',cols(ci,:),'LineWidth',1.9,'MarkerSize',5);
        lg=[lg,configs(ci,1)]; %#ok<AGROW>
    end
    if ~isempty(otfs_ber)
        y=otfs_ber; y(y==0)=NaN;
        semilogy(otfs_snr,y,'-d','Color',[0 0.45 0.74],'LineWidth',2.1,'MarkerSize',6);
        lg=[lg,{'CP-OTFS (MP), existing'}];
    end
    set(gca,'YScale','log'); xlabel('SNR [dB]'); ylabel('BER');
    dopptxt = 'nominal Doppler ~444 Hz'; if nu_scale==0, dopptxt='NO Doppler (nu=0, ICI-free)'; elseif nu_scale~=1, dopptxt=sprintf('Doppler x%g',nu_scale); end
    title(sprintf(['OFDM: channel estimation vs perfect CSI, reduced EVA, %s (%d-QAM, %d trials)\n' ...
                   'gap above black = CE error; a floor here would be residual ICI'],dopptxt,M_order,nTrials));
    legend(lg,'Location','southwest','FontSize',9); set(gca,'FontSize',11);
    if exist('theme','file'), theme(hf,"light"); end
    png=fullfile(resdir,['ber_ofdm_ce_vs_perfect' dtag '.png']);
    exportgraphics(hf,png,'Resolution',150,'BackgroundColor','white');
    fprintf('\nSaved %s and %s\n',png,outpath); fprintf('Done in %.0fs\n',toc(tAll));
end

function H=true_diag_gain(N,M,deltaf,Ts,alpha,tau,nu)
% True per-cell diagonal channel gain that multiplies X[n,m] (includes the ICI
% self-coefficient I0(nu); the off-diagonal ICI leakage stays in Y as residual).
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

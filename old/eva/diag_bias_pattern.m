% diag_bias_pattern.m
% Why does the real-detector MCRB differ from the BSC-matched one, and why does
% the sign flip with SNR (and hold for OTFS too)?  Compare the linearized
% score-bias F = A^-1 b b^T A^-1 (delay trace, evaluated at eta_bar -- NO
% pseudo-true solve) for the REAL error pattern vs a BSC-matched one (same count,
% i.i.d. placement).  F ~ ||A^-1 b||^2 is the first-order proxy for BB.
%   ratio F_real/F_bsc > 1  -> real errors bias MORE  (coherent)
%                       < 1  -> real errors bias LESS  (deep-fade nulling)
clear; clc;
thisdir = fileparts(mfilename('fullpath'));  addpath(fullfile(thisdir,'utils'));
warning('off','MATLAB:nearlySingularMatrix');

N=512; M=128; deltaf=15e3; Ts=1/deltaf; fc=4e9; c=3e8; v=120/3.6; nu_max=v*fc/c;
dt=round([0 150 370]*1e-9*N*deltaf); tau=(dt/(N*deltaf)).'; amp=(10.^([0 -1.4 -0.6]/20)).'; amp=amp/norm(amp);
P=numel(tau); l_tau=max(dt); nE=4*P; idx_tau=4*(0:P-1)+3;
M_order=4; constellation=qammod((0:3).',4,UnitAveragePower=true); n_iter_mp=30; damping=0.3;
l_p=l_tau; k_p=M/2; gr=1:(2*l_tau+1); cer=(l_p+1):(2*l_tau+1); dmO=true(N,M); dmO(gr,:)=false;
otfs_nd=numel(gr)*M; [dNp,dMp,ncp]=match_pilot_spacing(N,M,otfs_nd,l_tau,ceil(nu_max*M/deltaf));
np=1:dNp:N; mp=1:dMp:M; pm=false(N,M); pm(np,mp)=true; dmF=~pm;
E_pil=ncp; a_p=sqrt(E_pil/ncp); x_p=sqrt(E_pil);
Af=exp(-1j*2*pi*(0:N-1).'*(0:l_tau)/N); Ab=exp(-1j*2*pi*(np(:)-1)*(0:l_tau)/N);

Ftr = @(A,b) sum(real(diag( real(pinv(A)*(b*b.')*pinv(A)) )),1);   % crude; refined below

fprintf('   config   SNR |  F_tau(real)  F_tau(BSC)  ratio | coh_real coh_bsc | |H|@err\n');
for cfg = {'ofdm','otfs'}
  modn = cfg{1};
  for snr = [6 12 20]
    if strcmp(modn,'otfs') && snr>12, continue; end     % OTFS has ~0 errors above
    N0=1/db2pow(snr); nt=40;
    Fr=0; Fb=0; cohr=0; cohb=0; hAtErr=0; rng(11);
    for t=1:nt
      nu=(nu_max*cos(2*pi*rand(1,P)-pi)).'; alpha=amp.*exp(1j*2*pi*rand(P,1)); eb=packEta(alpha,tau,nu);
      sym=randi(4,N,M)-1; Xd=qammod(sym,4,UnitAveragePower=true); Z=sqrt(N0/2)*(randn(N,M)+1j*randn(N,M));
      Hcell = ofdm_mean(ones(N,M),deltaf,Ts,alpha,tau,nu,true);
      if strcmp(modn,'ofdm')
        d=Xd; d(pm)=a_p; Y=ofdm_mean(d,deltaf,Ts,alpha,tau,nu,true)+Z;
        Hls=Y(np,mp)./d(np,mp); Hh=Af*(Ab\Hls); FH=griddedInterpolant({(1:N).',mp},Hh,'linear','nearest');
        [Nq,Mq]=ndgrid(1:N,1:M); Hhat=FH(Nq,Mq);
        for it=1:2, Xs=conj(Hhat).*Y./(abs(Hhat).^2+N0); dh=qammod(qamdemod(Xs,4,UnitAveragePower=true),4,UnitAveragePower=true); dh(pm)=d(pm); Hhat=Af*(Af\(Y./dh)); end
        Xs=conj(Hhat).*Y./(abs(Hhat).^2+N0); dh=d; dh(dmF)=qammod(qamdemod(Xs(dmF),4,UnitAveragePower=true),4,UnitAveragePower=true);
        em=dmF; Xtf=d;
      else
        d=zeros(N,M); d(dmO)=Xd(dmO); d(l_p+1,k_p+1)=x_p; Gt=build_effective_dd_eva(alpha,dt,nu,N,M,deltaf,Inf);
        Y=apply_dd_channel_eva(d,Gt)+Z; hce=zeros(N,M); hce(cer,:)=Y(cer,:)/x_p; hce(abs(hce)<3*sqrt(N0)/x_p)=0;
        [rc,cc,gc]=find(hce); dh=d;
        if ~isempty(gc), dl=mod(rc-(l_p+1),N); dk=mod(cc-(k_p+1),M); kn=zeros(N,M); kn(l_p+1,k_p+1)=x_p;
          xh=mp_detector(Y,dl,dk,gc,N0,constellation,n_iter_mp,dmO,kn,damping); dh(dmO)=xh(dmO); end
        em=dmO; Xtf=isfft_u(d);
      end
      dh_b = bsc_match(d,dh,sym,em,M_order);
      if strcmp(modn,'ofdm'), Xh=dh; Xhb=dh_b; else, Xh=isfft_u(dh); Xhb=isfft_u(dh_b); end
      [A ,~,~,b ]=getAB_OFDM_core(Xtf,Xh ,deltaf,Ts,1,alpha,tau,nu,true,alpha,tau,nu);
      [Ab2,~,~,bb]=getAB_OFDM_core(Xtf,Xhb,deltaf,Ts,1,alpha,tau,nu,true,alpha,tau,nu);
      Cr=sq(A,b*b.'); Cb=sq(Ab2,bb*bb.');
      Fr=Fr+sum(real(diag(Cr(idx_tau,idx_tau)))); Fb=Fb+sum(real(diag(Cb(idx_tau,idx_tau))));
      % coherence proxy: |sum E*H|^2 / (Nerr * mean|E*H|^2)  (1 = incoherent, >1 coherent)
      er=(dh~=d)&em; eb2=(dh_b~=d)&em;
      cohr=cohr+coh((d(er)-dh(er)).*Hcell(er)); cohb=cohb+coh((d(eb2)-dh_b(eb2)).*Hcell(eb2));
      hAtErr=hAtErr+mean(abs(Hcell(er)).^2)/mean(abs(Hcell(em)).^2);
    end
    fprintf('  %-6s  %3d | %11.3e %11.3e  %5.2f | %7.2f %7.2f | %5.3f\n', ...
            modn, snr, Fr/nt, Fb/nt, (Fr)/(Fb), cohr/nt, cohb/nt, hAtErr/nt);
  end
end

function C = sq(A,BB)   % A^-1 BB A^-1, Jacobi-preconditioned pinv
  d=1./sqrt(abs(diag(A))); d(~isfinite(d))=1; D=diag(d); As=(D*A*D+(D*A*D).')/2; Ai=pinv(As);
  C=real(D*(Ai*(D*BB*D)*Ai)*D);
end
function r = coh(v)     % |sum v|^2 / (n * mean|v|^2): 1=random walk, ->n=coherent
  n=numel(v); if n==0, r=NaN; return; end; r = abs(sum(v))^2 / (n*mean(abs(v).^2)+eps);
end
function dh_b = bsc_match(d,dh,sym,em,Mo)
  ne=nnz(dh(em)~=d(em)); dh_b=d; if ne>0, di=find(em); pe=di(randperm(numel(di),ne));
    w=mod(sym(pe)+randi(Mo-1,ne,1),Mo); dh_b(pe)=qammod(w,Mo,UnitAveragePower=true); end
end
function eta=packEta(a,t,n), eta=reshape([real(a(:)).';imag(a(:)).';t(:).';n(:).'],[],1); end

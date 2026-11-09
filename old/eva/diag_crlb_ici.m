function diag_crlb_ici()
% DIAG_CRLB_ICI  Clean-case (genie, no bit errors) CRLB with the ICI kernel
% MODELLED vs IGNORED, for the reduced-EVA OFDM sensing model.  Answers: does
% modelling ICI raise or lower the delay / Doppler CRLB, and how does it scale
% with the time aperture M?  (CP-OFDM here is ISI-free by construction, so only
% the ICI half is testable numerically; ISI is its time<->frequency dual.)
%
% Genie FIM: Jt = (2/sigma2) Re{G^H G} from getAB_OFDM_core(X,X,...,use_ici).
% CRB = Jt^-1.  Averaged over random QPSK data frames (expected FIM).  Same fixed
% channel (alpha,tau,nu) for ici on/off so only the ICI signature differs.
    thisdir = fileparts(mfilename('fullpath')); addpath(fullfile(thisdir,'utils'));
    warning('off','MATLAB:nearlySingularMatrix'); warning('off','MATLAB:singularMatrix');

    N = 512; deltaf = 15e3; Ts = 1/deltaf; fc = 4e9; c = 3e8; v = 120/3.6;
    nu_nom = v*fc/c;                                   % ~444 Hz
    eva_ns = [0 150 370]; eva_db = [0 -1.4 -0.6]; P = numel(eva_ns);
    tau = (round(eva_ns*1e-9*N*deltaf)/(N*deltaf)).';
    amp = (10.^(eva_db/20)).'; amp = amp/norm(amp);
    idx_tau = 4*(0:P-1)+3; idx_nu = 4*(0:P-1)+4;
    nAvg = 24; seed = 12345;

    % fixed channel (drawn once): real Doppler so ICI is physically present
    rs = RandStream('Threefry','Seed',seed); RandStream.setGlobalStream(rs);
    nu_unit = cos(2*pi*rand(1,P)-pi).';               % in [-1,1], scaled below
    alpha   = amp .* exp(1j*2*pi*rand(P,1));

    fprintf('Clean genie CRLB: ICI modelled vs ignored (reduced EVA, N=%d, %d frames)\n', N, nAvg);
    fprintf('%-5s %-9s | %-26s | %-26s\n','M','nu_max', ...
            'DELAY RMSE [s]  ici|noici|ratio','DOPPLER RMSE [Hz] ici|noici|ratio');
    for nu_scale = [1 10]
        nu = nu_scale * nu_nom * nu_unit;
        for M = [128 8 1]
            [rt_i, rn_i] = crb_case(N,M,deltaf,Ts,alpha,tau,nu,true ,idx_tau,idx_nu,nAvg,seed);
            [rt_n, rn_n] = crb_case(N,M,deltaf,Ts,alpha,tau,nu,false,idx_tau,idx_nu,nAvg,seed);
            fprintf('%-5d %6.0fHz | %9.3e %9.3e %5.2f | %9.3e %9.3e %s\n', ...
                M, nu_scale*nu_nom, rt_i, rt_n, rt_i/rt_n, ...
                rn_i, rn_n, ratstr(rn_i, rn_n));
        end
    end
    fprintf(['\nRatio = ici/noici.  <1: modelling ICI LOWERS that CRLB (ICI adds info);\n' ...
             '>1: ICI raises it (coupling nuisance).  Inf noici Doppler at M=1 => Doppler\n' ...
             'is observable ONLY through ICI when the inter-symbol aperture is gone.\n']);
end

function [rt, rn] = crb_case(N,M,deltaf,Ts,alpha,tau,nu,use_ici,idx_tau,idx_nu,nAvg,seed)
    P = numel(alpha); nE = 4*P; Jt = zeros(nE);
    rs = RandStream('Threefry','Seed',seed+7); RandStream.setGlobalStream(rs);
    for k = 1:nAvg
        X = qammod(randi(4,N,M)-1, 4, UnitAveragePower=true);
        [~,~,Jk,~] = getAB_OFDM_core(X, X, deltaf, Ts, 1, alpha, tau, nu, use_ici);
        Jt = Jt + Jk;
    end
    Jt = Jt/nAvg; Jt = (Jt+Jt.')/2;
    % Jacobi-equilibrated inverse (same hardening as stab_quad); detect unobservable.
    d = 1./sqrt(abs(diag(Jt))); d(~isfinite(d)) = 1; D = diag(d);
    Js = D*Jt*D;
    if rcond(Js) < 1e-14              % rank-deficient (e.g. Doppler unobservable)
        CRB = D*pinv(Js)*D;
        obs_nu = abs(diag(Jt(idx_nu,idx_nu))) > 1e-9*max(abs(diag(Jt)));
    else
        CRB = D*(Js\eye(nE))*D;
        obs_nu = true(numel(idx_nu),1);
    end
    rt = sqrt(sum(real(diag(CRB(idx_tau,idx_tau)))));
    if all(obs_nu)
        rn = sqrt(sum(real(diag(CRB(idx_nu,idx_nu)))));
    else
        rn = Inf;                     % at least one Doppler direction unobservable
    end
end

function s = ratstr(a,b)
    if isinf(b) && isfinite(a), s = 'Inf(enabled)';
    elseif isinf(a), s = 'unobs';
    else, s = sprintf('%.2f', a/b); end
end

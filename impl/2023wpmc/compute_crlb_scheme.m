function out = compute_crlb_scheme(scheme, paths, sys, SNRdB, nMC, seed)
% COMPUTE_CRLB_SCHEME  CRLB of delay/Doppler for one signal model, using the
% paper's closed-form channel matrices Psi^p (eq.(7),(11),(26),(40)).
%
%   scheme : 'cpfree_ofdm' | 'cp_ofdm' | 'zak_otfs' | 'twosteps_otfs'
%   paths  : struct array, fields h (complex gain), tau (s), nu (Hz)
%   sys    : from sys_params
%   SNRdB  : vector of SNR points (dB).  Noise model y = Psi x + n,
%            n ~ CN(0, sigma2 I) in the estimation domain (eq.(46)),
%            sigma2 = 10^(-SNRdB/10).
%   nMC    : number of random unit-power 4-QAM data realisations to average the
%            FIM over (1 -> single fixed realisation).  Default 1.
%   seed   : RNG seed.  Default 2023.
%
% Per path the unknowns are theta = [|h|, angle(h), tau, nu] (Sec. III).
% FIM (Slepian-Bangs, parameter-independent covariance, eq.(48)):
%   J_ij = (2/sigma2) Re{ (dmu/dtheta_i)^H (dmu/dtheta_j) },  mu = Psi x.
% d/d|h| = Psi/|h|, d/dphi = j Psi (analytic); d/dtau, d/dnu by central finite
% differences of the exact closed-form Psi (identical to eq.(30)-(36),(43)-(44),
% (51)-(52); tau stays inside one I1/I2 piece for the step used here).

if nargin < 5 || isempty(nMC),  nMC = 1;    end
if nargin < 6 || isempty(seed), seed = 2023; end

P  = numel(paths);
NM = sys.NM;
dtau = 1e-9;        % finite-difference step for d/dtau  (s)
dnu  = 1e-3;        % finite-difference step for d/dnu    (Hz)

% ---- derivative matrices dPsi/dtheta_i (data independent) ---------------
dMat = cell(1,4*P);
pmap = zeros(4*P,2);                 % [path_index, type]  3=tau, 4=nu
col = 0;
for p = 1:P
    h = paths(p).h; tau = paths(p).tau; nu = paths(p).nu;
    ha = abs(h);
    Psi = psi_build(scheme,h,tau,nu,sys);
    % amplitude / phase (Psi linear in h = |h| e^{j phi})
    col=col+1; dMat{col}=Psi/ha; pmap(col,:)=[p 1];
    col=col+1; dMat{col}=1j*Psi;  pmap(col,:)=[p 2];
    % delay (central difference)
    Ptp = psi_build(scheme,h,tau+dtau,nu,sys);
    Ptm = psi_build(scheme,h,tau-dtau,nu,sys);
    col=col+1; dMat{col}=(Ptp-Ptm)/(2*dtau); pmap(col,:)=[p 3];
    % Doppler (central difference)
    Pnp = psi_build(scheme,h,tau,nu+dnu,sys);
    Pnm = psi_build(scheme,h,tau,nu-dnu,sys);
    col=col+1; dMat{col}=(Pnp-Pnm)/(2*dnu); pmap(col,:)=[p 4];
end
nP = 4*P;

% ---- FIM kernel G (averaged over data) ----------------------------------
%   G_ij = Re{ (D_i x)^H (D_j x) }, averaged over nMC unit-power 4-QAM blocks
%   (nMC -> inf gives the E{xx^H}=I form Re{trace(D_i^H D_j)}).
rng(seed);
qam = [1+1j, 1-1j, -1+1j, -1-1j]/sqrt(2);
G = zeros(nP);
for mc = 1:nMC
    x = qam(randi(4,NM,1)).';
    V = zeros(NM,nP);
    for a = 1:nP, V(:,a) = dMat{a}*x; end
    G = G + real(V'*V);
end
G = G/nMC;

% ---- Jacobi (diagonal) preconditioning + inversion ----------------------
% theta and h-parameters differ by ~20 orders of magnitude; precondition so
% Gt = (d d').*G has unit diagonal, with [G^{-1}]_ii = d_i^2 [Gt^{-1}]_ii.
d = 1./sqrt(diag(G));
Gt = (d*d').*G;
GinvDiag = (d.^2).*real(diag(inv(Gt)));

out.tau    = zeros(P,numel(SNRdB));
out.nu     = zeros(P,numel(SNRdB));
out.condGt = cond(Gt);
for si = 1:numel(SNRdB)
    sigma2 = 10^(-SNRdB(si)/10);
    for p = 1:P
        rt = pmap(:,1)==p & pmap(:,2)==3;
        rn = pmap(:,1)==p & pmap(:,2)==4;
        out.tau(p,si) = (sigma2/2)*GinvDiag(rt);
        out.nu(p,si)  = (sigma2/2)*GinvDiag(rn);
    end
end
end

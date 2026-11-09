% verify_varpro.m
% Sanity check for the variable-projection pseudo-true solve: it must reach the
% same optimum as the full 4P lsqnonlin, at far lower cost.  Compares eta0 and
% times both, across the four {OFDM,OTFS} x {no-ICI,ICI} configs at a few p.
clear; clc;
thisdir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisdir, 'utils'));
rng(7);

deltaf = 15e3;  Ts = 1/deltaf;
N = 128;  M = 64;
fc = 2e9;  v = 120/3.6;  c = 3e8;  numax = v*fc/c;

% reduced EVA, 3 targets (same style as compute_mcrb_mc.m)
tau   = [0; 370; 150]*1e-9;
nu    = numax*cos([0.5; 2.0; 3.5]);
alpha = 10.^([0; -0.6; -1.4]/20) .* exp(1j*[0.7; 1.1; -0.5]);

configs = {'OFDM no-ICI','ofdm',false; 'OFDM ICI','ofdm',true; ...
           'OTFS no-ICI','otfs',false; 'OTFS ICI','otfs',true};
p_list = [0, 0.1, 0.3];

qpsk = @(N,M) (1/sqrt(2))*((2*bitand(randi([0 3],N,M),1)-1) + 1j*(2*floor(randi([0 3],N,M)/2)-1));
bsc  = @(X,p) ( sign(real(X)).*(1-2*(rand(size(X))<p)) + 1j*sign(imag(X)).*(1-2*(rand(size(X))<p)) )/sqrt(2);

% Compare in scale-aware units: gains as-is, delays/Dopplers in resolution bins
% (tau bin = 1/(N*deltaf), nu bin = 1/(M*Ts)), so a single number is meaningful
% across the tau~1e-7 / nu~1e2 / alpha~1 disparity.
s_tau = 1/(N*deltaf);  s_nu = 1/(M*Ts);
fprintf('%-12s %5s | %9s %9s %9s | %9s %9s | %6s\n', ...
        'config','p','d_gain','d_tau[bin]','d_nu[bin]','t_full','t_vp','speedup');
maxdev = 0;
for ci = 1:size(configs,1)
    md = configs{ci,2};  ic = configs{ci,3};
    for p = p_list
        d  = qpsk(N,M);  dh = bsc(d, p);
        if strcmp(md,'ofdm')
            f_full = @() pseudoTrue_OFDM(d,dh,deltaf,Ts,alpha,tau,nu,ic,[],false);
            f_vp   = @() pseudoTrue_OFDM(d,dh,deltaf,Ts,alpha,tau,nu,ic,[],true);
        else
            f_full = @() pseudoTrue_CP_OTFS(d,dh,deltaf,Ts,alpha,tau,nu,ic,[],false);
            f_vp   = @() pseudoTrue_CP_OTFS(d,dh,deltaf,Ts,alpha,tau,nu,ic,[],true);
        end
        t1=tic; e_full = f_full(); tf=toc(t1);
        t2=tic; e_vp   = f_vp();   tv=toc(t2);

        de = e_full - e_vp;
        dg = norm([de(1:4:end); de(2:4:end)]);   % gain (aR,aI)
        dt = norm(de(3:4:end)) / s_tau;          % delay, in bins
        dn = norm(de(4:4:end)) / s_nu;           % Doppler, in bins
        maxdev = max(maxdev, max([dg dt dn]));
        fprintf('%-12s %5.2f | %9.2e %9.2e %9.2e | %8.3fs %8.3fs | %5.1fx\n', ...
                configs{ci,1}, p, dg, dt, dn, tf, tv, tf/max(tv,eps));
    end
end
fprintf('\nmax scale-aware deviation (gain / tau-bin / nu-bin) = %.3e\n', maxdev);
if maxdev < 1e-3
    fprintf('PASS: VARPRO matches the full 4P optimum (well within one resolution bin).\n');
else
    fprintf('WARN: deviation larger than expected -- inspect.\n');
end

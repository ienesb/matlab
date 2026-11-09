% run_mcrb_all.m
% Master driver: computes MCRB component matrices (A, B, J_tilde, b, A_Delta),
% pseudo-true parameter eta_0, and MCRB(tau) / MCRB(nu) for all four cases:
%   {OFDM-ICI, OTFS} x {single target (first EVA path), multi-target (reduced EVA)}
% Sweeps BSC bit-error probability rho, averages over MC trials, saves
% everything to mcrb_results.mat.

clear; clc;
addpath('../crb/EVA');

% ---- Grid / channel parameters ----
N      = 512;
M      = 128;
deltaf = 15e3;
T      = 1/deltaf;
Tcp    = 0;
Ts     = T + Tcp;
sigma2 = 1;

c        = 3e8;
fc       = 4e9;
velocity = 120 * 1e3 / 3600;
nu_max   = velocity * fc / c;

% ---- Single-target: first EVA path (tau=0, |alpha|=1) ----
rng(0);
single.taus   = 0;
single.alphas = exp(1j*2*pi*rand);            % random phase, unit magnitude
single.nus    = nu_max * cos(2*pi*rand);
single.P      = 1;

% ---- Multi-target: reduced EVA (P = 3) ----
time_delays_ns  = [0,  370,  150];
relative_powers = [0,  -0.6, -1.4];
multi.P = numel(time_delays_ns);
multi.taus = time_delays_ns * 1e-9;
P_lin = 10.^(relative_powers/10);
P_lin = P_lin / sum(P_lin);
amps  = sqrt(P_lin);

rng(1);
phases = 2*pi*rand(multi.P,1);
thetas = 2*pi*rand(multi.P,1);
multi.alphas = amps(:) .* exp(1j*phases);
multi.nus    = nu_max * cos(thetas);

% ---- Sweep ----
rho_vals = 0:0.05:0.5;
nMC      = 10;

% ---- lsqnonlin options (used for eta_0) ----
opts = optimoptions('lsqnonlin','Display','off', ...
                    'Algorithm','trust-region-reflective', ...
                    'MaxIterations',60, 'MaxFunctionEvaluations',4000, ...
                    'FunctionTolerance',1e-12, 'StepTolerance',1e-12);

cases    = {'ofdm_single','ofdm_multi','otfs_single','otfs_multi'};
modtypes = {'ofdm',       'ofdm',      'otfs',       'otfs'      };
targets  = {single,       multi,       single,       multi       };

% Initialize result containers
results = struct;
for ic = 1:numel(cases)
    mn = cases{ic};
    P_local = targets{ic}.P;
    results.(mn).A       = cell(numel(rho_vals), nMC);
    results.(mn).B       = cell(numel(rho_vals), nMC);
    results.(mn).J_tilde = cell(numel(rho_vals), nMC);
    results.(mn).b       = cell(numel(rho_vals), nMC);
    results.(mn).A_Delta = cell(numel(rho_vals), nMC);
    results.(mn).eta_0   = cell(numel(rho_vals), nMC);
    results.(mn).MCRB_tau = zeros(numel(rho_vals), nMC, P_local);
    results.(mn).MCRB_nu  = zeros(numel(rho_vals), nMC, P_local);
    results.(mn).eta_bar  = build_eta(targets{ic}.alphas, targets{ic}.taus, targets{ic}.nus);
    results.(mn).target   = targets{ic};
end

t_start = tic;
for ir = 1:numel(rho_vals)
    rho = rho_vals(ir);
    for mc = 1:nMC
        rng(1000*ir + mc);   % reproducible per (rho, MC)
        bits = randi([0 3], N, M);
        X    = (1/sqrt(2)) * (2*bitand(bits,1) - 1 + 1j*(2*floor(bits/2) - 1));
        X_hat = bsc_qpsk(X, rho);

        for ic = 1:numel(cases)
            mn   = cases{ic};
            mod  = modtypes{ic};
            tgt  = targets{ic};
            P_local = tgt.P;

            % ---- MCRB component matrices ----
            if strcmp(mod, 'ofdm')
                out = mcrb_ofdm(X, X_hat, deltaf, Ts, sigma2, tgt.alphas, tgt.taus, tgt.nus);
            else
                out = mcrb_otfs(X, X_hat, deltaf, Ts, sigma2, tgt.alphas, tgt.taus, tgt.nus);
            end

            % ---- MCRB = A^-1 B A^-T  (A symmetric -> A^-T = A^-1) ----
            d = sqrt(abs(diag(out.A)));
            S = d * d.';
            invA = (inv(out.A ./ S)) ./ S;
            MCRB = invA * out.B * invA.';

            % ---- pseudo-true eta_0 ----
            eta_bar_vec    = results.(mn).eta_bar;
            scales         = repmat([1; 1; 1e9; 1], P_local, 1);
            eta_bar_scaled = eta_bar_vec .* scales;
            box            = repmat([2; 2; 1000; 2*nu_max], P_local, 1);
            lb = eta_bar_scaled - box;
            ub = eta_bar_scaled + box;

            if strcmp(mod, 'ofdm')
                mu_bar_eval = mu_ofdm_at(X, deltaf, Ts, tgt.alphas, tgt.taus, tgt.nus);
                res = @(es) residual_ofdm(es, scales, mu_bar_eval, X_hat, deltaf, Ts);
            else
                mu_bar_eval = mu_otfs_at(X, deltaf, Ts, tgt.alphas, tgt.taus, tgt.nus);
                res = @(es) residual_otfs(es, scales, mu_bar_eval, X_hat, deltaf, Ts);
            end
            eta0_scaled = lsqnonlin(res, eta_bar_scaled, lb, ub, opts);
            eta0_vec    = eta0_scaled ./ scales;

            % ---- store ----
            results.(mn).J_tilde{ir,mc} = out.J_tilde;
            results.(mn).b{ir,mc}       = out.b;
            results.(mn).A_Delta{ir,mc} = out.A_Delta;
            results.(mn).A{ir,mc}       = out.A;
            results.(mn).B{ir,mc}       = out.B;
            results.(mn).eta_0{ir,mc}   = eta0_vec;
            for p = 1:P_local
                results.(mn).MCRB_tau(ir,mc,p) = MCRB(4*(p-1)+3, 4*(p-1)+3);
                results.(mn).MCRB_nu (ir,mc,p) = MCRB(4*(p-1)+4, 4*(p-1)+4);
            end
        end
    end
    elapsed = toc(t_start);
    fprintf('rho = %.2f done   (elapsed %.1f s)\n', rho, elapsed);
end

% ---- MC-average MCRBs and eta_0 bias ----
for ic = 1:numel(cases)
    mn = cases{ic};
    results.(mn).MCRB_tau_avg = squeeze(mean(results.(mn).MCRB_tau, 2));
    results.(mn).MCRB_nu_avg  = squeeze(mean(results.(mn).MCRB_nu,  2));
    eta_bar = results.(mn).eta_bar;
    bias_norms = zeros(numel(rho_vals), nMC);
    for ir = 1:numel(rho_vals)
        for mc = 1:nMC
            bias_norms(ir,mc) = norm(results.(mn).eta_0{ir,mc} - eta_bar);
        end
    end
    results.(mn).bias_norm_avg = mean(bias_norms, 2);
end

% ---- Save ----
save('mcrb_results.mat', 'results', 'rho_vals', 'single', 'multi', ...
     'N','M','deltaf','Ts','sigma2','nu_max', '-v7.3');
fprintf('Saved mcrb_results.mat\n');

% ---- Quick console summary ----
fprintf('\n%-14s  | path |   rho   |   MCRB(tau)   |   MCRB(nu)   \n', 'case');
fprintf('%s\n', repmat('-',1,75));
for ir = [1, ceil(numel(rho_vals)/2), numel(rho_vals)]
    for ic = 1:numel(cases)
        mn = cases{ic};
        for p = 1:targets{ic}.P
            fprintf('%-14s  |  %d   |  %.2f   |  %.3e    |  %.3e\n', mn, p, ...
                rho_vals(ir), results.(mn).MCRB_tau_avg(ir,p), ...
                results.(mn).MCRB_nu_avg (ir,p));
        end
    end
end

% ================== local helpers ==================
function v = build_eta(alphas, taus, nus)
    P = numel(alphas);
    v = zeros(4*P, 1);
    for p = 1:P
        v(4*(p-1) + (1:4)) = [real(alphas(p)); imag(alphas(p)); taus(p); nus(p)];
    end
end

function mu = mu_ofdm_at(X, deltaf, Ts, alphas, taus, nus)
    [N, M] = size(X);
    T = 1/deltaf;
    P = numel(alphas);
    ns = (0:N-1).';
    ms = (0:M-1);
    ds = (-(N-1):(N-1)).';
    mu = zeros(N, M);
    for p = 1:P
        A      = X .* exp(-1j*2*pi*deltaf*ns*taus(p)) .* exp(1j*2*pi*Ts*ms*nus(p));
        I_full = get_ici_coeffs(ds, nus(p), T);
        I_mat  = toeplitz(I_full(N:-1:1), I_full(N:end));
        mu     = mu + alphas(p) * (I_mat * A);
    end
end

function mu = mu_otfs_at(x, deltaf, Ts, alphas, taus, nus)
    [N, M] = size(x);
    P = numel(alphas);
    ns   = (0:N-1).';
    ms   = (0:M-1).';
    ds_A = (-(N-1):(N-1)).';
    ds_B = (-(M-1):(M-1)).';
    mu = zeros(N, M);
    for p = 1:P
        E_A = exp(1j*2*pi*ns*(ds_A.'/N - deltaf*taus(p)));
        A_d = sum(E_A,1).';
        E_B = exp(1j*2*pi*ms*(ds_B.'/M + Ts*nus(p)));
        B_d = sum(E_B,1).';
        A_mat = toeplitz(A_d(N:end), A_d(N:-1:1));
        B_mat = toeplitz(B_d(M:end), B_d(M:-1:1));
        mu = mu + (alphas(p)/(N*M)) * (A_mat * x * B_mat);
    end
end

function r = residual_ofdm(eta_s, scales, mu_bar, X_hat, deltaf, Ts)
    eta    = eta_s ./ scales;
    alphas = complex(eta(1:4:end), eta(2:4:end));
    taus   = eta(3:4:end);
    nus    = eta(4:4:end);
    mu     = mu_ofdm_at(X_hat, deltaf, Ts, alphas, taus, nus);
    d      = mu_bar(:) - mu(:);
    r      = [real(d); imag(d)];
end

function r = residual_otfs(eta_s, scales, mu_bar, x_hat, deltaf, Ts)
    eta    = eta_s ./ scales;
    alphas = complex(eta(1:4:end), eta(2:4:end));
    taus   = eta(3:4:end);
    nus    = eta(4:4:end);
    mu     = mu_otfs_at(x_hat, deltaf, Ts, alphas, taus, nus);
    d      = mu_bar(:) - mu(:);
    r      = [real(d); imag(d)];
end

function [eta0, info] = pseudo_true(X, Xhat, g, cr_true, use_ici, opts)
% PSEUDO_TRUE  Pseudo-true parameter of the misspecified CP-OFDM / CP-OTFS model
% (latex/mcrb/pseudo_true.tex; cp_ofdm_mcrb.tex, cp_otfs_mcrb.tex Sec. pseudo-true):
%
%   eta0 = argmin_eta || mu(X, eta_bar) - mu_tilde(Xhat, eta) ||^2
%
% i.e. the least-squares projection of the true noise-free mean onto the
% manifold of assumed means built from the DETECTED symbols.
%
%   [eta0, info] = pseudo_true(X, Xhat, g, cr_true, use_ici, opts)
% opts (defaults):
%   method   'varpro' : variable projection over [tau; nu] in bin units with
%                       the gains concentrated (alpha = C \ mu) and the Kaufman
%                       analytic Jacobian (old/eva/utils/pseudoTrue_OFDM.m);
%            'full4P' : lsqnonlin over all 4P unknowns (slow, for tests);
%            'closed' : expectation closed form [(1-2p) alpha_bar; tau_bar; nu_bar]
%                       (needs opts.p);
%   domain   'tf' (X, Xhat are TF grids) | 'dd' (DD grids -> ISFFT first)
%   eta_init [] (default: the true parameter; the closed form with opts.p if given)
%   mask_nu  auto (true when all cr_true.nu == 0): nu is not estimated, held at nu_bar
%   p        BSC rate for the closed-form warm start
%   solver   'lsqnonlin' (Optimization Toolbox) | 'lm' (built-in Levenberg-Marquardt)
% info: resnorm, iterations, exitflag, method, seconds.
    % Tolerances are set for a trial-averaged bound: eta0 to ~1e-6 bins is far
    % below the pseudo-true jitter (~1e-2 bins), and the misspecified residual
    % is large, so Gauss-Newton converges linearly and over-tight tolerances only
    % add iterations (old/eva used 1e-10/1e-9/1e-8).
    d = struct('method', 'varpro', 'domain', 'tf', 'eta_init', [], 'mask_nu', [], 'p', [], 'solver', 'auto', ...
        'ftol', 1e-9, 'steptol', 1e-8, 'opttol', 1e-7, 'maxiter', 200);
    if nargin < 6 || isempty(opts), opts = struct(); end
    fn = fieldnames(d);  for i = 1:numel(fn), if ~isfield(opts, fn{i}), opts.(fn{i}) = d.(fn{i}); end, end
    t0 = tic;
    if strcmpi(opts.domain, 'dd'), X = isfft_dd(X);  Xhat = isfft_dd(Xhat); end
    [N, M] = size(X);
    P = numel(cr_true.alpha);
    alpha_bar = cr_true.alpha(:);  tau_bar = cr_true.tau(:);  nu_bar = cr_true.nu(:);
    if isempty(opts.mask_nu), opts.mask_nu = all(nu_bar == 0); end
    if strcmpi(opts.solver, 'auto')
        if exist('lsqnonlin', 'file') == 2, opts.solver = 'lsqnonlin'; else, opts.solver = 'lm'; end
    end

    if strcmpi(opts.method, 'closed')
        assert(~isempty(opts.p), 'pseudo_true:p', 'closed form needs opts.p');
        eta0 = eta_pack((1 - 2*opts.p)*alpha_bar, tau_bar, nu_bar);
        info = struct('method', 'closed', 'resnorm', NaN, 'iterations', 0, 'exitflag', 1, 'seconds', toc(t0));
        return;
    end

    mu_true = io_cp_ofdm(X, g, cr_true, use_ici);  mu_true = mu_true(:);
    s_tau = g.bin_tau;  s_nu = g.bin_nu;

    if isempty(opts.eta_init)
        if ~isempty(opts.p), a_init = (1 - 2*opts.p)*alpha_bar; else, a_init = alpha_bar; end
        eta_init = eta_pack(a_init, tau_bar, nu_bar);
    else
        eta_init = opts.eta_init(:);
    end
    [~, tinit, ninit] = eta_unpack(eta_init);

    switch lower(opts.method)
        case 'varpro'
            free = [true(P, 1); ~opts.mask_nu & true(P, 1)];        % [tau; nu] unknowns
            sc2  = [s_tau*ones(P, 1); s_nu*ones(P, 1)];
            th_full0 = [tinit(:); ninit(:)] ./ sc2;
            fun = @(th) varpro_resjac(Xhat, g, th, th_full0, free, sc2, P, use_ici, mu_true);
            th0 = th_full0(free);
            [th_opt, resnorm, exitflag, iters] = solve_ls(fun, th0, opts);
            th_full = th_full0;  th_full(free) = th_opt;
            thu = th_full .* sc2;
            tau0 = thu(1:P);  nu0 = thu(P+1:2*P);
            [C, ~, ~] = basis_jac(Xhat, g, tau0, nu0, use_ici, false);
            alpha0 = C \ mu_true;
            eta0 = eta_pack(alpha0, tau0, nu0);
        case 'full4p'
            sc = repmat([1; 1; s_tau; s_nu], P, 1);
            th0 = eta_init ./ sc;
            free = true(4*P, 1);  if opts.mask_nu, free(4:4:end) = false; end
            fun = @(th) full_resjac(Xhat, g, th, th0, free, sc, use_ici, mu_true);
            [th_opt, resnorm, exitflag, iters] = solve_ls(fun, th0(free), opts);
            thf = th0;  thf(free) = th_opt;
            eta0 = thf .* sc;
        otherwise
            error('pseudo_true:method', 'unknown method %s', opts.method);
    end
    info = struct('method', opts.method, 'resnorm', resnorm, 'iterations', iters, 'exitflag', exitflag, ...
        'seconds', toc(t0), 'solver', opts.solver, 'mask_nu', opts.mask_nu);
end

% ---- variable-projection residual and Kaufman Jacobian (scaled unknowns) ----
function [r, J] = varpro_resjac(Xhat, g, th, th_full0, free, sc2, P, use_ici, mu_true)
    th_full = th_full0;  th_full(free) = th;
    thu = th_full .* sc2;
    tau = thu(1:P);  nu = thu(P+1:2*P);
    [C, dCt, dCn] = basis_jac(Xhat, g, tau, nu, use_ici, true);
    a   = C \ mu_true;
    rho = mu_true - C*a;
    r   = [real(rho); imag(rho)];
    if nargout > 1
        Gd = [dCt .* a.', dCn .* a.'];               % NM x 2P
        Jc = -(Gd - C*(C \ Gd));                     % -(I - C C^+) Gd
        Jc = Jc .* sc2.';                            % chain rule to bin units
        Jc = Jc(:, free);
        J  = [real(Jc); imag(Jc)];
    end
end

function [r, J] = full_resjac(Xhat, g, th, th0, free, sc, use_ici, mu_true)
    thf = th0;  thf(free) = th;
    eta = thf .* sc;
    [alpha, tau, nu] = eta_unpack(eta);
    P = numel(alpha);
    [C, dCt, dCn] = basis_jac(Xhat, g, tau, nu, use_ici, nargout > 1);
    rho = mu_true - C*alpha;
    r = [real(rho); imag(rho)];
    if nargout > 1
        Jc = zeros(size(C, 1), 4*P);
        for p = 1:P
            col = 4*(p-1);
            Jc(:, col+1) = -C(:, p);
            Jc(:, col+2) = -1j*C(:, p);
            Jc(:, col+3) = -alpha(p)*dCt(:, p);
            Jc(:, col+4) = -alpha(p)*dCn(:, p);
        end
        Jc = Jc .* sc.';
        Jc = Jc(:, free);
        J = [real(Jc); imag(Jc)];
    end
end

% ---- gain-stripped basis C(:,p) = vec(Phi_p) and its (tau, nu) derivatives ----
function [C, dCt, dCn] = basis_jac(S, g, tau, nu, use_ici, want_d)
    [N, M] = size(S);  P = numel(tau);
    C = zeros(N*M, P);  dCt = zeros(N*M, P);  dCn = zeros(N*M, P);
    for p = 1:P
        F = fields_cp_ofdm(S, g, tau(p), nu(p), use_ici, double(want_d));
        C(:, p) = F.Phi(:);
        if want_d
            dCt(:, p) = -1j*2*pi * F.Psi(:);
            dCn(:, p) =  1j*2*pi * F.Omega(:);
        end
    end
end

% ---- least-squares driver: lsqnonlin if available, else built-in LM ----------
function [th, resnorm, exitflag, iters] = solve_ls(fun, th0, opts)
    if strcmpi(opts.solver, 'lsqnonlin')
        o = optimoptions('lsqnonlin', 'Display', 'off', 'SpecifyObjectiveGradient', true, ...
            'FunctionTolerance', opts.ftol, 'StepTolerance', opts.steptol, 'OptimalityTolerance', opts.opttol, ...
            'MaxFunctionEvaluations', 200*numel(th0), 'MaxIterations', opts.maxiter);
        [th, resnorm, ~, exitflag, out] = lsqnonlin(fun, th0, [], [], o);
        iters = out.iterations;
        return;
    end
    % Levenberg-Marquardt with the analytic Jacobian (no toolbox).  Stops on a
    % small gradient (relative to the initial one), a small step, or a
    % negligible relative decrease of the sum of squares.
    th = th0;  lambda = 1e-3;
    [r, J] = fun(th);  f = r.'*r;  iters = 0;  exitflag = 0;
    g0 = norm(J.'*r);
    for it = 1:opts.maxiter
        iters = it;
        H = J.'*J;  gvec = J.'*r;
        if norm(gvec) <= opts.opttol*max(g0, eps), exitflag = 1; break; end
        step = -(H + lambda*diag(diag(H))) \ gvec;
        thn = th + step;
        [rn, Jn] = fun(thn);  fnew = rn.'*rn;
        if fnew < f
            df = (f - fnew)/max(f, eps);
            th = thn;  r = rn;  J = Jn;  f = fnew;  lambda = max(lambda/10, 1e-12);
            if norm(step) < opts.steptol*(1 + norm(th)) || df < opts.ftol, exitflag = 1; break; end
        else
            lambda = lambda*10;
            % Exhausted damping is a failed/stalled solve, not convergence.
            % All MCRB drivers accept only positive exit flags.
            if lambda > 1e12, exitflag = -2; break; end
        end
    end
    resnorm = f;
end

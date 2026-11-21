function [eta_hat, info] = mmle_cp(Y, Xhat, g, P, use_ici, opts)
% MMLE_CP  Misspecified maximum-likelihood estimator for CP-OFDM / CP-OTFS
% (latex/mmle/ofdm_mmle.tex, cp_otfs_mmle.tex): the gain-concentrated
% subspace-projection search on the DETECTED symbols Xhat,
%   theta_hat = argmax_theta y^H C (C^H C)^-1 C^H y,  alpha_hat = C^+ y,
% C = [c_1 ... c_P], c_p = vec(Phi_p(Xhat; tau_p, nu_p)).  CLEAN initialization
% (zero-padded 2-D FFT periodogram + fminsearch refine per target, successive
% cancellation) then joint refinement of all (tau, nu) in bin units.
%   opts: Zpad (4), domain 'tf' | 'dd' (ISFFT first), mask_nu (auto when P
%   targets are delay-only is NOT inferred -- pass true to search tau only).
%   opts.method 'nm' (default; the port of old/eva/utils/mmle_OFDM.m: Nelder-Mead
%   per target and jointly, ~90 s per 512 x 128 frame with P = 3) | 'lsq'
%   (same objective, since argmax of the projected energy = argmin ||y - C alpha||^2:
%   periodogram peak + single-target Gauss-Newton per CLEAN step, then the joint
%   variable-projection least-squares fit of pseudo_true with y as the target;
%   ~1 pseudo-true solve).  opts.eta_init (4P x 1, 'lsq' only) skips the CLEAN
%   start, e.g. the true parameter for a locally initialized estimator.
% Returns eta_hat (4P x 1) in CLEAN (descending strength) order; associate to
% the true targets by nearest (tau, nu) in bin units before computing an MSE.
% info ('lsq'): tau0, nu0 (start), exitflag, resnorm, iterations of the joint fit.
% Ported from old/eva/utils/mmle_OFDM.m.
    d = struct('Zpad', 4, 'domain', 'tf', 'mask_nu', false, 'method', 'nm', 'eta_init', []);
    if nargin < 6 || isempty(opts), opts = struct(); end
    fn = fieldnames(d);  for i = 1:numel(fn), if ~isfield(opts, fn{i}), opts.(fn{i}) = d.(fn{i}); end, end
    if strcmpi(opts.domain, 'dd'), Y = isfft_dd(Y);  Xhat = isfft_dd(Xhat); end
    [N, M] = size(Y);
    s_tau = g.bin_tau;  s_nu = g.bin_nu;
    if strcmpi(opts.method, 'lsq'), [eta_hat, info] = mmle_lsq(Y, Xhat, g, P, use_ici, opts);  return; end
    assert(strcmpi(opts.method, 'nm'), 'mmle_cp:method', 'unknown method %s', opts.method);
    steer = @(tau, nu) reshape(getfield(fields_cp_ofdm(Xhat, g, tau, nu, use_ici, 0), 'Phi'), [], 1);
    tau0 = zeros(P, 1);  nu0 = zeros(P, 1);
    R = Y;
    fmopt = optimset('Display', 'off', 'TolX', 1e-6, 'TolFun', 1e-9, 'MaxFunEvals', 400, 'MaxIter', 400);
    for p = 1:P
        r  = R .* conj(Xhat);
        Nz = opts.Zpad*N;  Mz = opts.Zpad*M;
        Gp = fft(ifft(r, Nz, 1), Mz, 2);
        [~, idx] = max(abs(Gp(:)).^2);
        [ia, ib] = ind2sub([Nz, Mz], idx);
        a = ia - 1;  if a > Nz/2, a = a - Nz; end
        b = ib - 1;  if b > Mz/2, b = b - Mz; end
        if opts.mask_nu
            obj = @(th) -mf_metric(R, steer(th(1)*s_tau, 0));
            th = fminsearch(obj, a/opts.Zpad, fmopt);  tau0(p) = th*s_tau;  nu0(p) = 0;
        else
            obj = @(th) -mf_metric(R, steer(th(1)*s_tau, th(2)*s_nu));
            th = fminsearch(obj, [a/opts.Zpad; b/opts.Zpad], fmopt);
            tau0(p) = th(1)*s_tau;  nu0(p) = th(2)*s_nu;
        end
        c = steer(tau0(p), nu0(p));
        ap = (c'*R(:)) / (c'*c);
        R = R - reshape(ap*c, N, M);
    end
    if P > 1
        jopt = optimset('Display', 'off', 'TolX', 1e-6, 'TolFun', 1e-9, 'MaxFunEvals', 2000*P, 'MaxIter', 2000*P);
        if opts.mask_nu
            jobj = @(th) -proj_energy(Y, steer, th(:)*s_tau, zeros(P, 1));
            thJ = fminsearch(jobj, tau0/s_tau, jopt);
            tau_hat = thJ(:)*s_tau;  nu_hat = zeros(P, 1);
        else
            jobj = @(th) -proj_energy(Y, steer, th(1:P)*s_tau, th(P+1:2*P)*s_nu);
            thJ = fminsearch(jobj, [tau0/s_tau; nu0/s_nu], jopt);
            tau_hat = thJ(1:P)*s_tau;  nu_hat = thJ(P+1:2*P)*s_nu;
        end
    else
        tau_hat = tau0;  nu_hat = nu0;
    end
    C = zeros(N*M, P);
    for p = 1:P, C(:, p) = steer(tau_hat(p), nu_hat(p)); end
    alpha_hat = (C'*C) \ (C'*Y(:));
    eta_hat = eta_pack(alpha_hat, tau_hat, nu_hat);
    info = struct('tau0', tau0, 'nu0', nu0);
end

function [eta_hat, info] = mmle_lsq(Y, Xhat, g, P, use_ici, opts)
% Gauss-Newton version (Y, Xhat already in TF).  Every fit is pseudo_true's
% variable-projection solver with the observation (or the CLEAN residual) as target.
    [N, M] = size(Y);
    po = struct('domain', 'tf', 'mask_nu', opts.mask_nu);
    if isempty(opts.eta_init)
        a0 = zeros(P, 1);  tau0 = zeros(P, 1);  nu0 = zeros(P, 1);
        R = Y;  Nz = opts.Zpad*N;  Mz = opts.Zpad*M;
        for p = 1:P
            Gp = fft(ifft(R .* conj(Xhat), Nz, 1), Mz, 2);
            if opts.mask_nu, Gp = Gp(:, 1); end
            [~, idx] = max(abs(Gp(:)).^2);
            [ia, ib] = ind2sub(size(Gp), idx);
            a = ia - 1;  if a > Nz/2, a = a - Nz; end
            b = ib - 1;  if b > Mz/2, b = b - Mz; end
            c1 = struct('alpha', 1, 'tau', a/opts.Zpad*g.bin_tau, 'nu', b/opts.Zpad*g.bin_nu);
            po.target = R;
            e1 = pseudo_true(Xhat, Xhat, g, c1, use_ici, po);
            [a0(p), tau0(p), nu0(p)] = eta_unpack(e1);
            Fp = fields_cp_ofdm(Xhat, g, tau0(p), nu0(p), use_ici, 0);
            R = R - a0(p)*Fp.Phi;
        end
    else
        [a0, tau0, nu0] = eta_unpack(opts.eta_init(:));
    end
    po.target = Y;
    [eta_hat, pinfo] = pseudo_true(Xhat, Xhat, g, struct('alpha', a0, 'tau', tau0, 'nu', nu0), use_ici, po);
    info = struct('tau0', tau0, 'nu0', nu0, 'exitflag', pinfo.exitflag, 'resnorm', pinfo.resnorm, ...
        'iterations', pinfo.iterations, 'method', 'lsq');
end

function v = mf_metric(Y, c)
    v = abs(c'*Y(:))^2 / (real(c'*c) + eps);
end

function v = proj_energy(Y, steer, taus, nus)
    P = numel(taus);  C = zeros(numel(Y), P);
    for p = 1:P, C(:, p) = steer(taus(p), nus(p)); end
    Cy = C'*Y(:);
    v = real(Cy' * ((C'*C) \ Cy));
end

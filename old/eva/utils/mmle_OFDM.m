function [eta_hat, alpha_hat, tau_hat, nu_hat] = ...
        mmle_OFDM(Y, Xhat, deltaf, Ts, P, use_ici, Zpad)
% Misspecified maximum-likelihood estimator (MMLE) for the OFDM ISAC model
% (ofdm_mmle.tex).  Estimates the P-target channel parameters
% eta = [aR_p; aI_p; tau_p; nu_p] from the noisy time-frequency observation Y,
% using the DETECTED symbols Xhat as the despreading reference -- the
% substitution Xhat (not the transmitted X) is exactly what makes the model
% misspecified (ofdm_mmle.tex sec. 2).
%
% Per ofdm_mmle.tex the assumed mean is linear in the gains,
%   mu_tilde = sum_p alpha_p c_tilde_p(tau_p,nu_p),
%   c_tilde[n,m](tau,nu) = sum_{n'} Xhat[n',m] e^{-j2pi n' df tau}
%                                              e^{ j2pi m Ts nu} I(n'-n),
% (I(n)=delta[n] in the no-ICI channel).  Concentrating the gains out gives the
% subspace-projection (periodogram) search, single target eqs (8)-(9) and
% multi-target eqs (18)-(19),
%   theta_hat = argmax_theta || P_{C(theta)} y ||^2 = y^H C (C^H C)^-1 C^H y,
%   alpha_hat = C^dagger y = (C^H C)^-1 C^H y,   C = [c_tilde_1 ... c_tilde_P].
% The P=1 case reduces to the single-target periodogram
%   (tau,nu)=argmax |c_tilde^H y|^2/||c_tilde||^2,  alpha = c_tilde^H y/||c_tilde||^2.
%
% Search strategy: each target is initialized by a successive (CLEAN-type)
% interference-cancellation pass (estimate the strongest target with the
% single-target periodogram, subtract its reconstruction from the residual,
% repeat -- ofdm_mmle.tex sec. 3.3), after which all P delay-Doppler pairs are
% refined JOINTLY on the projection energy.  The joint stage is essential when
% targets are sub-bin spaced, where CLEAN alone cannot resolve the overlapping
% steering vectors.
%
% Inputs:
%   Y       : N x M noisy time-frequency observation (true mean + noise)
%   Xhat    : N x M detected symbol grid (despreading reference)
%   deltaf  : subcarrier spacing [Hz]
%   Ts      : OFDM symbol duration incl. CP [s]
%   P       : number of targets to extract
%   use_ici : true -> ICI kernel I(n'-n); false -> no-ICI (I = delta)
%   Zpad    : (optional) FFT zero-pad factor for the coarse search; default 4
%
% Outputs:
%   eta_hat              : 4P x 1 stacked estimate [aR;aI;tau;nu] per target,
%                          in the CLEAN-initialized (descending-strength) order
%   alpha_hat,tau_hat,nu_hat : P x 1 gains / delays / Dopplers

    if nargin < 7 || isempty(Zpad), Zpad = 4; end
    [N, M] = size(Y);

    s_tau = 1/(N*deltaf);                 % one delay   bin [s]  (variable scale)
    s_nu  = 1/(M*Ts);                     % one Doppler bin [Hz]

    % ===== Stage 1: CLEAN initialization (per-target periodogram + refine) =====
    tau0 = zeros(P,1);  nu0 = zeros(P,1);
    R = Y;                                % residual
    for p = 1:P
        % coarse periodogram peak via a zero-padded 2-D FFT of r = R Xhat*
        % metric ~ |sum_{n,m} r e^{j2pi n df tau} e^{-j2pi m Ts nu}|.
        r  = R .* conj(Xhat);
        Nz = Zpad*N;  Mz = Zpad*M;
        Gp = fft(ifft(r, Nz, 1), Mz, 2);          % +exp over n, -exp over m
        [~, idx] = max(abs(Gp(:)).^2);
        [ia, ib] = ind2sub([Nz, Mz], idx);
        a = ia - 1;  if a > Nz/2, a = a - Nz; end  % signed bin index
        b = ib - 1;  if b > Mz/2, b = b - Mz; end

        % single-target refine with the EXACT (ICI-aware) steering vector
        obj = @(th) -mf_metric(R, Xhat, deltaf, Ts, th(1)*s_tau, th(2)*s_nu, use_ici);
        opt = optimset('Display','off', 'TolX',1e-6, 'TolFun',1e-9, ...
                       'MaxFunEvals',400, 'MaxIter',400);
        th  = fminsearch(obj, [a/Zpad; b/Zpad], opt);   % init in bin units
        tau0(p) = th(1)*s_tau;  nu0(p) = th(2)*s_nu;

        % subtract this target's reconstruction from the residual
        c  = steer(Xhat, deltaf, Ts, tau0(p), nu0(p), use_ici);
        ap = (c(:)' * R(:)) / (c(:)' * c(:));
        R  = R - ap*c;
    end

    % ===== Stage 2: joint refinement of all (tau,nu) on the projection energy =
    if P > 1
        jobj = @(th) -proj_energy(Y, Xhat, deltaf, Ts, ...
                                  th(1:P).*s_tau, th(P+1:2*P).*s_nu, use_ici);
        jopt = optimset('Display','off', 'TolX',1e-6, 'TolFun',1e-9, ...
                        'MaxFunEvals',2000*P, 'MaxIter',2000*P);
        thJ  = fminsearch(jobj, [tau0./s_tau; nu0./s_nu], jopt);
        tau_hat = thJ(1:P).*s_tau;  nu_hat = thJ(P+1:2*P).*s_nu;
    else
        tau_hat = tau0;  nu_hat = nu0;
    end

    % ===== gains by least squares at the final delay-Doppler estimates =========
    C = zeros(N*M, P);
    for p = 1:P
        cp = steer(Xhat, deltaf, Ts, tau_hat(p), nu_hat(p), use_ici);
        C(:,p) = cp(:);
    end
    alpha_hat = (C' * C) \ (C' * Y(:));

    eta_hat = packEta(alpha_hat, tau_hat, nu_hat);
end

% ----- single-target matched-filter metric |c^H y|^2 / ||c||^2 -----------
function v = mf_metric(Y, Xhat, deltaf, Ts, tau, nu, use_ici)
    c  = steer(Xhat, deltaf, Ts, tau, nu, use_ici);
    cv = c(:);
    v  = abs(cv' * Y(:))^2 / (real(cv' * cv) + eps);
end

% ----- joint projection energy ||P_C y||^2 = y^H C (C^H C)^-1 C^H y -------
function v = proj_energy(Y, Xhat, deltaf, Ts, taus, nus, use_ici)
    [N, M] = size(Y);  P = numel(taus);
    C = zeros(N*M, P);
    for p = 1:P
        cp = steer(Xhat, deltaf, Ts, taus(p), nus(p), use_ici);
        C(:,p) = cp(:);
    end
    yv = Y(:);
    Cy = C' * yv;
    v  = real(Cy' * ((C' * C) \ Cy));
end

% ----- detected-symbol steering vector c_tilde(tau,nu) (ofdm_mmle.tex eq 6) --
function c = steer(Xhat, deltaf, Ts, tau, nu, use_ici)
    [N, M] = size(Xhat);
    ns = (0:N-1).';  ms = (0:M-1);
    phase = exp(-1j*2*pi*deltaf*ns*tau) * exp(1j*2*pi*Ts*ms*nu);   % N x M
    A = Xhat .* phase;
    if use_ici
        T = 1/deltaf;  [NP, NN] = ndgrid(0:N-1, 0:N-1);  diff_n = NP - NN;
        c = get_ici_coeffs(diff_n, nu, T).' * A;          % sum_{n'} A I(n'-n)
    else
        c = A;                                            % I = delta
    end
end

function eta = packEta(alpha, tau, nu)
    P = numel(alpha);  eta = zeros(4*P, 1);
    eta(1:4:end) = real(alpha);  eta(2:4:end) = imag(alpha);
    eta(3:4:end) = tau;          eta(4:4:end) = nu;
end

function [eta0, alpha0, tau0, nu0, resnorm, exitflag] = ...
        pseudoTrue_OFDM(X, Xhat, deltaf, Ts, alpha_bar, tau_bar, nu_bar, use_ici, eta_init, use_varpro)
% Pseudo-true parameter for the misspecified OFDM ISAC model (pseudo_true.tex).
% Solves the least-squares projection (last equation of the document)
%
%   eta0 = argmin_eta || mu(X, eta_bar) - mu_tilde(Xhat, eta) ||,
%
% i.e. the projection of the true noise-free mean (built from the transmitted
% symbols X and the TRUE parameters eta_bar) onto the manifold of assumed means
% (built from the detected symbols Xhat and the FREE parameter eta).  The means
% are the OFDM I/O means of ofdm_mcrb.tex eqs (3)-(4):
%
%   mu[n,m] = sum_p alpha_p sum_{n'} S[n',m] e^{-j2pi n' df tau_p}
%                                           e^{ j2pi m Ts nu_p} I_p(n'-n),
%
% with S = X (true mean) or S = Xhat (assumed mean).  All four cases are handled
% through (numel(alpha_bar), use_ici): single/multi target and no-ICI/ICI, where
% I_p(n) = delta[n] (no-ICI) or the unit-peak ICI kernel get_ici_coeffs (ICI).
%
% Parameter stacking: eta = [eta_1; ...; eta_P], eta_p = [aR_p; aI_p; tau_p; nu_p].
%
% Inputs:
%   X, Xhat   : N x M transmitted / detected symbol grids
%   deltaf,Ts : subcarrier spacing [Hz], OFDM symbol duration [s]
%   alpha_bar : P x 1 true complex path gains   (define the projection target)
%   tau_bar   : P x 1 true delays   [s]
%   nu_bar    : P x 1 true Dopplers [Hz]
%   use_ici   : true -> ICI kernel, false -> no-ICI (I = delta)
%   eta_init  : (optional) 4P x 1 initial guess; default = the true parameter
%   use_varpro: (optional) true -> variable-projection (gain-concentrated) solve
%               over the 2P delay/Doppler unknowns only, with the P complex gains
%               obtained in closed form by linear least squares at each (tau,nu)
%               (ofdm_mcrb.tex / cp_otfs_mcrb.tex Sec. "Pseudo-true parameter",
%               eqs. gain-concentration + cross-ambiguity, exact per-realization
%               form -- NO data averaging). false (default) -> full 4P solve.
%               Both reach the same optimum; VARPRO removes the ill-conditioned
%               gain block from the nonlinear problem and halves the unknowns.
%
% Outputs:
%   eta0            : 4P x 1 pseudo-true parameter (stacked [aR;aI;tau;nu])
%   alpha0,tau0,nu0 : P x 1 unpacked pseudo-true gains / delays / Dopplers
%   resnorm         : || mu - mu_tilde ||^2 at the optimum
%   exitflag        : lsqnonlin exit flag
%
% Solved with lsqnonlin (Gauss-Newton on the stacked real/imag residual).  The
% delay/Doppler unknowns are internally rescaled to ~one resolution bin so the
% optimizer operates on O(1) variables (tau ~ 1e-7 s, nu ~ 1e2 Hz, alpha ~ 1).

    alpha_bar = alpha_bar(:);  tau_bar = tau_bar(:);  nu_bar = nu_bar(:);
    P = numel(alpha_bar);

    % Projection target: true noise-free mean from X and the true parameters.
    mu_true = ofdm_mean(X, deltaf, Ts, alpha_bar, tau_bar, nu_bar, use_ici);
    mu_true = mu_true(:);

    % Variable scaling: delay/Doppler in resolution-bin units, gains O(1).
    [N, M] = size(X);
    s_tau = 1/(N*deltaf);                       % ~ one delay bin   [s]
    s_nu  = 1/(M*Ts);                           % ~ one Doppler bin [Hz]
    sc = repmat([1; 1; s_tau; s_nu], P, 1);     % 4P x 1

    if nargin < 9 || isempty(eta_init)
        eta_init = packEta(alpha_bar, tau_bar, nu_bar);
    end
    if nargin < 10 || isempty(use_varpro), use_varpro = false; end

    if use_varpro
        % ---- Variable projection: optimize (tau,nu) only, gains concentrated ----
        % For fixed (tau,nu) the objective is quadratic in the P complex gains, so
        % the gain-optimal alpha is the linear LS projection alpha0 = Chat \ mu_true
        % (ofdm_mcrb.tex eq. pt-gain / multi-target Gram). Substituting it back
        % leaves the reduced residual r(tau,nu) = (I - Chat Chat^+) mu_true over the
        % 2P delay/Doppler unknowns.  Exact per-realization form (no data averaging).
        [~, tinit, ninit] = unpackEta(eta_init);
        sc2 = [s_tau*ones(P,1); s_nu*ones(P,1)];     % 2P x 1 (tau block, nu block)
        th0 = [tinit(:); ninit(:)] ./ sc2;

        % Analytic Jacobian (Kaufman variable-projection form) supplied via the
        % two-output objective, so lsqnonlin does NOT finite-difference the
        % residual (which would re-evaluate the expensive ICI basis 2P+1 times
        % per iteration).  Tolerances relaxed from 1e-12/1e-15: the pseudo-true
        % feeds a trial-averaged bound, so sub-1e-8 precision on eta0 is wasted.
        resjac = @(th) varpro_resjac(Xhat, deltaf, Ts, th, sc2, P, use_ici, mu_true);
        opts = optimoptions('lsqnonlin', 'Display','off', ...
                            'SpecifyObjectiveGradient',true, ...
                            'FunctionTolerance',1e-10, 'StepTolerance',1e-9, ...
                            'OptimalityTolerance',1e-8, ...
                            'MaxFunctionEvaluations',200*(2*P), 'MaxIterations',200);
        [th_opt, resnorm, ~, exitflag] = lsqnonlin(resjac, th0, [], [], opts);

        thu   = th_opt .* sc2;
        tau0  = thu(1:P);  nu0 = thu(P+1:2*P);
        Chat0 = ofdm_basis(Xhat, deltaf, Ts, tau0, nu0, use_ici);   % NM x P
        alpha0 = Chat0 \ mu_true;                                   % concentrated gains
        eta0  = packEta(alpha0, tau0, nu0);
        [alpha0, tau0, nu0] = unpackEta(eta0);       % normalize orientation/order
    else
        % ---- Full 4P solve over [aR;aI;tau;nu] per target ----
        th0 = eta_init ./ sc;                        % scaled initial guess

        resfun = @(th) ri( mu_true - vecMean(Xhat, deltaf, Ts, th.*sc, use_ici) );
        opts = optimoptions('lsqnonlin', 'Display','off', ...
                            'FunctionTolerance',1e-15, 'StepTolerance',1e-12, ...
                            'OptimalityTolerance',1e-12, ...
                            'MaxFunctionEvaluations',1000*(4*P), 'MaxIterations',500);
        [th_opt, resnorm, ~, exitflag] = lsqnonlin(resfun, th0, [], [], opts);

        eta0 = th_opt .* sc;
        [alpha0, tau0, nu0] = unpackEta(eta0);
    end
end

% ----- variable-projection reduced residual + analytic Jacobian ----------
function [r, J] = varpro_resjac(Xhat, deltaf, Ts, th, sc2, P, use_ici, mu_true)
% Real-stacked residual r = ri(rho), rho = (I - C C^+) mu_true, and its Jacobian
% w.r.t. the scaled th = [tau; nu]./sc2.  Kaufman variable-projection form:
%   d rho / d theta_k = -a_p (I - C C^+) (dC[:,p]/d theta_k),   a = C \ mu_true,
% where only column p of C depends on theta_k, and
%   dPhi_p/dtau_p = -j2pi Psi_p,   dPhi_p/dnu_p = j2pi Omega_p
% with Psi_p, Omega_p the same kernel products as getAB_OFDM_core (valid for
% both no-ICI and ICI; the ICI dImat/dnu term enters Omega via D1).
    thu = th .* sc2;
    tau = thu(1:P);  nu = thu(P+1:2*P);
    [C, dCt, dCn] = ofdm_basis_jac(Xhat, deltaf, Ts, tau, nu, use_ici);  % NM x P each
    a   = C \ mu_true;                          % concentrated gains (P x 1)
    rho = mu_true - C*a;
    r   = [real(rho); imag(rho)];
    if nargout > 1
        Gd  = [dCt .* a.', dCn .* a.'];         % NM x 2P : (dC/dtheta_k) a
        Jc  = -(Gd - C*(C \ Gd));               % -(I - C C^+) Gd  (unscaled)
        Jc  = Jc .* sc2.';                       % chain rule to scaled th
        J   = [real(Jc); imag(Jc)];             % 2NM x 2P
    end
end

% ----- gain-stripped basis and its (tau,nu) column derivatives -----------
function [C, dCt, dCn] = ofdm_basis_jac(S, deltaf, Ts, tau, nu, use_ici)
% C(:,p)   = vec(Phi_p);
% dCt(:,p) = vec(dPhi_p/dtau_p) = -j2pi vec(Psi_p);
% dCn(:,p) = vec(dPhi_p/dnu_p)  =  j2pi vec(Omega_p).
    [N, M] = size(S);  P = numel(tau);
    ns = (0:N-1).';  ms = (0:M-1);
    C = zeros(N*M, P);  dCt = zeros(N*M, P);  dCn = zeros(N*M, P);
    if use_ici
        T = 1/deltaf;  [NP, NN] = ndgrid(0:N-1, 0:N-1);  diff_n = NP - NN;
    end
    for p = 1:P
        phase = exp(-1j*2*pi*deltaf*ns*tau(p)) * exp(1j*2*pi*Ts*ms*nu(p));  % N x M
        Ahat  = S .* phase;
        if use_ici
            Imat = get_ici_coeffs           (diff_n, nu(p), T);
            Ip   = get_ici_coeffs_derivative (diff_n, nu(p), T) / (1j*2*pi);   % I'
            Phi  = Imat.' * Ahat;
            Psi  = Imat.' * ((ns*deltaf) .* Ahat);
            D1   = Ip.'   * Ahat;
            Om   = (ms*Ts) .* Phi + D1;
        else
            Phi = Ahat;
            Psi = (ns*deltaf) .* Phi;
            Om  = (ms*Ts) .* Phi;
        end
        C(:,p)   = Phi(:);
        dCt(:,p) = -1j*2*pi * Psi(:);
        dCn(:,p) =  1j*2*pi * Om(:);
    end
end

% ----- mean builder (ofdm_mcrb.tex eqs (3)-(4)) --------------------------
function mu = ofdm_mean(S, deltaf, Ts, alpha, tau, nu, use_ici)
% Noise-free mean mu = sum_p alpha_p Phi_p = reshape(ofdm_basis * alpha).
    [N, M] = size(S);
    C  = ofdm_basis(S, deltaf, Ts, tau, nu, use_ici);   % NM x P
    mu = reshape(C * alpha(:), N, M);
end

% ----- gain-stripped basis: column p = vec(Phi_p(tau_p,nu_p)) ------------
function C = ofdm_basis(S, deltaf, Ts, tau, nu, use_ici)
    [N, M] = size(S);  P = numel(tau);
    ns = (0:N-1).';  ms = (0:M-1);  C = zeros(N*M, P);
    if use_ici
        T = 1/deltaf;  [NP, NN] = ndgrid(0:N-1, 0:N-1);  diff_n = NP - NN;
    end
    for p = 1:P
        phase = exp(-1j*2*pi*deltaf*ns*tau(p)) * exp(1j*2*pi*Ts*ms*nu(p));   % N x M
        A = S .* phase;
        if use_ici
            Imat = get_ici_coeffs(diff_n, nu(p), T);     % (n',n) = I(n'-n)
            Phi  = Imat.' * A;                           % sum_{n'} A I(n'-n)
        else
            Phi  = A;                                    % I = delta
        end
        C(:,p) = Phi(:);
    end
end

function v = vecMean(S, deltaf, Ts, eta, use_ici)
    [alpha, tau, nu] = unpackEta(eta);
    mu = ofdm_mean(S, deltaf, Ts, alpha, tau, nu, use_ici);
    v = mu(:);
end

function r = ri(z),  r = [real(z); imag(z)];  end

function eta = packEta(alpha, tau, nu)
    P = numel(alpha);  eta = zeros(4*P, 1);
    eta(1:4:end) = real(alpha);  eta(2:4:end) = imag(alpha);
    eta(3:4:end) = tau;          eta(4:4:end) = nu;
end

function [alpha, tau, nu] = unpackEta(eta)
    alpha = eta(1:4:end) + 1j*eta(2:4:end);
    tau   = eta(3:4:end);
    nu    = eta(4:4:end);
end

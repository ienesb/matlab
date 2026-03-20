function X_hat = mp_detector(Y_DD, delay_shifts, doppler_shifts, gains, ...
                              N0, constellation, n_iter, data_mask, ...
                              known_symbols, damping)
% MP_DETECTOR  Message-passing detector using DD channel triplets.
%
% Follows Raviteja et al., WCNC 2018 (Sec IV) and TVT 2019 (Sec II-D).
% The channel is represented as a sparse set of (delay, Doppler, gain)
% triplets, which can come from either perfect CSI or channel estimation.
%
% Inputs:
%   Y_DD           [N x M]  received DD grid
%   delay_shifts   [S x 1]  integer delay shifts (0-indexed, for circshift)
%   doppler_shifts [S x 1]  integer Doppler shifts (0-indexed, for circshift)
%   gains          [S x 1]  complex channel gains h_tilde
%   N0             scalar   noise variance per complex sample
%   constellation  [Q x 1]  constellation points
%   n_iter         scalar   number of MP iterations
%   data_mask      [N x M]  logical, true = data position
%   known_symbols  [N x M]  pilot value at pilot pos, 0 at guard/unknown
%   damping        scalar   damping factor in (0,1], default 0.5

if nargin < 10, damping = 0.5; end

[N, M] = size(Y_DD);
S = length(gains);
Q = numel(constellation);
E_x2 = mean(abs(constellation).^2);

%% Initialize beliefs
mu  = known_symbols;          % pilot: known; guard: 0
v   = zeros(N, M);            % pilot/guard: variance = 0
mu(data_mask)  = 0;
v(data_mask)   = E_x2;

prev_post = ones(N, M, Q) / Q;

%% MP iterations
log_post = zeros(N, M, Q);

for iter = 1:n_iter

    %% Step A — Aggregate expected received signal at each observation
    total_mean = zeros(N, M);
    total_var  = zeros(N, M);
    for s = 1:S
        mu_sh = circshift(mu, [delay_shifts(s), doppler_shifts(s)]);
        v_sh  = circshift(v,  [delay_shifts(s), doppler_shifts(s)]);
        total_mean = total_mean + gains(s) * mu_sh;
        total_var  = total_var  + abs(gains(s))^2 * v_sh;
    end

    %% Step B — Log-likelihood accumulation (observation→variable messages)
    log_post(:) = 0;
    for s = 1:S
        h = gains(s);
        if abs(h) < 1e-12, continue; end

        % Shift observation grid to variable-node position
        Y_sh  = circshift(Y_DD,       [-delay_shifts(s), -doppler_shifts(s)]);
        tm_sh = circshift(total_mean,  [-delay_shifts(s), -doppler_shifts(s)]);
        tv_sh = circshift(total_var,   [-delay_shifts(s), -doppler_shifts(s)]);

        % Extrinsic: remove current path's contribution
        interf_mean = tm_sh - h * mu;
        interf_var  = max(tv_sh - abs(h)^2 * v, 0);

        % Effective observation at each variable node
        y_eff      = (Y_sh - interf_mean) / h;
        sigma2_eff = (N0 + interf_var) / abs(h)^2;

        for q = 1:Q
            log_post(:,:,q) = log_post(:,:,q) ...
                - abs(y_eff - constellation(q)).^2 ./ sigma2_eff;
        end
    end

    %% Step C — Update posteriors with damping
    lp_max = max(log_post, [], 3);
    post   = exp(log_post - lp_max);
    post   = post ./ sum(post, 3);

    % Damping (Raviteja WCNC 2018, Eq. 21)
    post = damping * post + (1 - damping) * prev_post;
    post = post ./ sum(post, 3);
    prev_post = post;

    % Update mean and variance at data positions
    mu_new = zeros(N, M);
    v_new  = zeros(N, M);
    for q = 1:Q
        p_q    = post(:,:,q);
        mu_new = mu_new + p_q * constellation(q);
        v_new  = v_new  + p_q * abs(constellation(q))^2;
    end
    v_new = max(v_new - abs(mu_new).^2, 0);

    mu(data_mask) = mu_new(data_mask);
    v(data_mask)  = v_new(data_mask);
end

%% Hard decision (Step 5)
[~, idx] = max(log_post, [], 3);
X_hat = zeros(N, M);
for q = 1:Q
    X_hat(data_mask & idx == q) = constellation(q);
end

end

function [mse_tau, mse_nu, mse_alpha] = mse_est( ...
        gt_tau, gt_nu, gt_alpha, ...
        est_tau, est_nu, est_alpha, ...
        deltaf, N, M)
% MSE_EST  Match estimated paths to ground-truth by (tau, nu) proximity,
%          then compute per-parameter MSE over matched pairs.
%
%   Matching is done greedily in normalised delay-Doppler bin space:
%       delay bin = tau * N * deltaf
%     Doppler bin = nu  * M / deltaf
%   so both axes have comparable scale.  At each step the globally cheapest
%   unassigned (GT, estimate) pair is assigned first (nearest-neighbour).
%
%   Inputs
%     gt_tau / gt_nu / gt_alpha   : ground-truth delay (s), Doppler (Hz),
%                                   complex gain  — length-P vectors
%     est_tau / est_nu / est_alpha: estimated counterparts — length-K vectors
%     deltaf                      : subcarrier spacing (Hz)
%     N                           : number of subcarriers (delay axis size)
%     M                           : number of symbols  (Doppler axis size)
%
%   Outputs
%     mse_tau   : mean |Δτ|²   over matched pairs  (s²)
%     mse_nu    : mean |Δν|²   over matched pairs  (Hz²)
%     mse_alpha : mean |Δα|²   over matched pairs

    gt_tau  = gt_tau(:);   gt_nu  = gt_nu(:);   gt_alpha  = gt_alpha(:);
    est_tau = est_tau(:);  est_nu = est_nu(:);  est_alpha = est_alpha(:);

    n_gt  = numel(gt_tau);
    n_est = numel(est_tau);

    if n_est == 0 || n_gt == 0
        mse_tau = NaN;  mse_nu = NaN;  mse_alpha = NaN;
        return;
    end

    % Cost matrix (n_gt × n_est): Euclidean distance in bin space
    d_tau = (gt_tau - est_tau') * N * deltaf;   % n_gt × n_est  [delay bins]
    d_nu  = (gt_nu  - est_nu')  * M / deltaf;   % n_gt × n_est  [Doppler bins]
    C = sqrt(d_tau.^2 + d_nu.^2);

    % Greedy nearest-neighbour matching
    matched_gt  = false(n_gt,  1);
    matched_est = false(n_est, 1);
    n_pairs = 0;
    pairs   = zeros(min(n_gt, n_est), 2);

    for k = 1:min(n_gt, n_est)
        Cm = C;
        Cm(matched_gt,  :) = Inf;
        Cm(:, matched_est) = Inf;
        [cost, idx] = min(Cm(:));
        if isinf(cost), break; end
        [r, c] = ind2sub([n_gt, n_est], idx);
        matched_gt(r)  = true;
        matched_est(c) = true;
        n_pairs = n_pairs + 1;
        pairs(n_pairs, :) = [r, c];
    end

    pairs = pairs(1:n_pairs, :);

    if n_pairs == 0
        mse_tau = NaN;  mse_nu = NaN;  mse_alpha = NaN;
        return;
    end

    gi = pairs(:, 1);
    ei = pairs(:, 2);

    mse_tau   = mean(abs(gt_tau(gi)   - est_tau(ei)).^2);
    mse_nu    = mean(abs(gt_nu(gi)    - est_nu(ei)).^2);
    mse_alpha = mean(abs(gt_alpha(gi) - est_alpha(ei)).^2);
end

function [Xhat, info] = mp_detect_dd(Y, taps, alphabet, N0, known_mask, known_val, n_iter, damp, damp_mode, tol)
% MP_DETECT_DD  Message-passing OTFS detector on the sparse delay-Doppler
% factor graph (Raviteja et al., WCNC 2018, Sec. IV, eqs. 18-21), delay = rows.
%
%   [Xhat, info] = mp_detect_dd(Y, taps, alphabet, N0, known_mask, known_val, ...
%                               n_iter, damp, damp_mode, tol)
%   Y          N x M received DD grid (delay rows, Doppler cols)
%   taps       struct dl, dk (0-based shifts), gain (S x 1): y[l+dl,k+dk] += gain*x[l,k]
%   alphabet   Q x 1 constellation (already scaled to the data amplitude)
%   N0         noise variance per cell
%   known_mask / known_val   pilot+guard cells (deterministic nodes)
%   n_iter     max iterations (default 20);  damp  Delta in (0,1] (default 0.7)
%   damp_mode  'edge' (default): per-edge messages, damping on the
%              variable->observation messages, WCNC eq. (21) -- port of
%              impl/otfs_ber/mp_detect_otfs.m transposed to delay rows;
%              'posterior': damping on the node posteriors -- the legacy
%              old/eva/utils/mp_detector.m (0.3, 30 iterations) for regression.
%   tol        early stop when max |Delta p| < tol (default 1e-4; 0 = never)
% info: iterations run, final max message change, S, memory estimate.
% Memory ('edge'): two N x M x S x Q single arrays -> 8*N*M*S*Q bytes.
    if nargin < 7 || isempty(n_iter), n_iter = 20; end
    if nargin < 8 || isempty(damp), damp = 0.7; end
    if nargin < 9 || isempty(damp_mode), damp_mode = 'edge'; end
    if nargin < 10 || isempty(tol), tol = 1e-4; end
    [N, M] = size(Y);
    S = numel(taps.gain);  Q = numel(alphabet);
    alphabet = alphabet(:);
    info.S = S;  info.mem_bytes = 8*N*M*S*Q;
    if S == 0
        Xhat = zeros(N, M);  Xhat(known_mask) = known_val(known_mask);  info.iters = 0;  return;
    end
    if strcmpi(damp_mode, 'posterior')
        [Xhat, info.iters] = mp_posterior(Y, taps, alphabet, N0, known_mask, known_val, n_iter, damp);
        return;
    end
    dl = taps.dl(:);  dk = taps.dk(:);  gS = single(taps.gain(:));  absg2 = single(abs(taps.gain(:)).^2);
    ds  = single(reshape(alphabet, 1, 1, 1, Q));
    ds2 = single(reshape(abs(alphabet).^2, 1, 1, 1, Q));
    Gsq = single(taps.gain(:) * alphabet.');            % S x Q
    floorv = single(N0*1e-6);
    hasKnown = any(known_mask(:));
    if hasKnown, kidx = find(known_mask);  kval = single(known_val(kidx)); end

    % Observation seen from variable (l,k) through tap s: y[l+dl, k+dk].
    Yshift = zeros(N, M, S, 'single');
    for s = 1:S, Yshift(:,:,s) = single(circshift(Y, [-dl(s), -dk(s)])); end

    P = ones(N, M, S, Q, 'single') / Q;                  % variable -> observation messages
    logL = zeros(N, M, S, Q, 'single');
    maxdiff = NaN;  it = 0;
    for it = 1:n_iter
        muX = sum(P .* ds, 4);                            % N x M x S
        vX  = sum(P .* ds2, 4) - abs(muX).^2;
        if hasKnown
            muX = reshape(muX, N*M, S);  muX(kidx, :) = repmat(kval, 1, S);  muX = reshape(muX, N, M, S);
            vX  = reshape(vX,  N*M, S);  vX(kidx, :)  = 0;                    vX  = reshape(vX,  N, M, S);
        end
        Mtot = zeros(N, M, 'single');  Vtot = zeros(N, M, 'single');
        for s = 1:S
            Mtot = Mtot + gS(s)    * circshift(muX(:,:,s), [dl(s), dk(s)]);
            Vtot = Vtot + absg2(s) * circshift(vX(:,:,s),  [dl(s), dk(s)]);
        end
        Vtot = Vtot + single(N0);
        sumLogL = zeros(N, M, Q, 'single');
        for s = 1:S
            mu_excl  = circshift(Mtot, [-dl(s), -dk(s)]) - gS(s)*muX(:,:,s);
            var_excl = max(circshift(Vtot, [-dl(s), -dk(s)]) - absg2(s)*vX(:,:,s), floorv);
            res = Yshift(:,:,s) - mu_excl;
            for q = 1:Q
                d  = res - Gsq(s, q);
                ll = -(real(d).^2 + imag(d).^2) ./ var_excl;
                logL(:,:,s,q) = ll;
                sumLogL(:,:,q) = sumLogL(:,:,q) + ll;
            end
        end
        maxdiff = 0;
        for s = 1:S
            lm = sumLogL - reshape(logL(:,:,s,:), N, M, Q);   % extrinsic (exclude edge s)
            lm = lm - max(lm, [], 3);
            w  = exp(lm);  w = w ./ sum(w, 3);
            Pnew = (1-damp)*P(:,:,s,:) + damp*reshape(w, N, M, 1, Q);
            maxdiff = max(maxdiff, max(abs(Pnew - P(:,:,s,:)), [], 'all'));
            P(:,:,s,:) = Pnew;
        end
        if tol > 0 && maxdiff < tol, break; end
    end
    [~, qhat] = max(sumLogL, [], 3);
    Xhat = reshape(alphabet(qhat), N, M);
    Xhat(known_mask) = known_val(known_mask);
    info.iters = it;  info.maxdiff = double(maxdiff);
end

function [Xhat, it] = mp_posterior(Y, taps, alphabet, N0, known_mask, known_val, n_iter, damp)
% Legacy posterior-damped MP (old/eva/utils/mp_detector.m), delay rows.
    [N, M] = size(Y);  S = numel(taps.gain);  Q = numel(alphabet);
    dl = taps.dl(:);  dk = taps.dk(:);  gains = taps.gain(:);
    data_mask = ~known_mask;
    mu = known_val;  v = zeros(N, M);
    mu(data_mask) = 0;  v(data_mask) = mean(abs(alphabet).^2);
    prev_post = ones(N, M, Q)/Q;  log_post = zeros(N, M, Q);
    for it = 1:n_iter
        total_mean = zeros(N, M);  total_var = zeros(N, M);
        for s = 1:S
            total_mean = total_mean + gains(s) * circshift(mu, [dl(s), dk(s)]);
            total_var  = total_var  + abs(gains(s))^2 * circshift(v, [dl(s), dk(s)]);
        end
        log_post(:) = 0;
        for s = 1:S
            h = gains(s);  if abs(h) < 1e-12, continue; end
            Y_sh  = circshift(Y, [-dl(s), -dk(s)]);
            tm_sh = circshift(total_mean, [-dl(s), -dk(s)]);
            tv_sh = circshift(total_var,  [-dl(s), -dk(s)]);
            interf_mean = tm_sh - h*mu;
            interf_var  = max(tv_sh - abs(h)^2*v, 0);
            y_eff = (Y_sh - interf_mean)/h;
            s2eff = (N0 + interf_var)/abs(h)^2;
            for q = 1:Q
                log_post(:,:,q) = log_post(:,:,q) - abs(y_eff - alphabet(q)).^2 ./ s2eff;
            end
        end
        post = exp(log_post - max(log_post, [], 3));  post = post ./ sum(post, 3);
        post = damp*post + (1-damp)*prev_post;  post = post ./ sum(post, 3);
        prev_post = post;
        mu_new = zeros(N, M);  v_new = zeros(N, M);
        for q = 1:Q
            mu_new = mu_new + post(:,:,q)*alphabet(q);
            v_new  = v_new  + post(:,:,q)*abs(alphabet(q))^2;
        end
        v_new = max(v_new - abs(mu_new).^2, 0);
        mu(data_mask) = mu_new(data_mask);  v(data_mask) = v_new(data_mask);
    end
    [~, idx] = max(log_post, [], 3);
    Xhat = known_val;
    Xd = reshape(alphabet(idx), N, M);
    Xhat(data_mask) = Xd(data_mask);
end

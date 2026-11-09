function Xhat = mp_detect_otfs(Y, taps, dataSyms, sigma2, knownMask, knownVal, n_iter, damp)
% MP_DETECT_OTFS  Vectorised message-passing detector for OTFS, following
% [MP] Sec. IV (eqs. 18-21).  Interference + noise seen by each symbol is
% approximated as Gaussian; messages are passed on the sparse DD factor
% graph defined by the channel taps.
%
% Inputs:
%   Y         : N x M received DD grid
%   taps      : struct with dk,dl (0-based shifts) and complex gain (S taps)
%   dataSyms  : Q x 1 data constellation (already scaled to sqrt(SNRd))
%   sigma2    : noise variance
%   knownMask : N x M logical, true where the symbol is a known pilot/guard
%   knownVal  : N x M, the known symbol values (0 on guards, x_p at pilot)
%   n_iter    : MP iterations
%   damp      : damping factor Delta (weight on the fresh message)
%
% Output:
%   Xhat      : N x M hard-decided data symbols (knownVal kept at known pos.)
%
% Big arrays are kept in single precision (memory-bandwidth bound).

[N, M] = size(Y);
S = numel(taps.gain);
Q = numel(dataSyms);
dk = taps.dk; dl = taps.dl; g = taps.gain;

ds   = single(reshape(dataSyms, 1, 1, 1, Q));     % 1x1x1xQ
ds2  = single(reshape(abs(dataSyms).^2, 1, 1, 1, Q));
Gsq  = g(:) * dataSyms.';                          % S x Q : gain(s)*dataSyms(q)
gS   = single(g(:)); absg2 = single(abs(g).^2);
floorv = single(sigma2*1e-6);

hasKnown = any(knownMask(:));
if hasKnown
    kidx = find(knownMask);
    kval = single(knownVal(kidx));                 % column of known values
end

% Precompute the (iteration-invariant) shifted observations y_e for each tap.
Yshift = zeros(N, M, S, 'single');
for s = 1:S
    Yshift(:,:,s) = single(circshift(Y, [-dk(s), -dl(s)]));
end

% Messages P(k,l,s,q): from variable (k,l) toward obs reached via tap s.
P    = ones(N, M, S, Q, 'single') / Q;
logL = zeros(N, M, S, Q, 'single');

for it = 1:n_iter
    % ---- (1) per-message symbol mean / variance --------------------
    muX = sum(P .* ds, 4);                          % NxMxS single complex
    vX  = sum(P .* ds2, 4) - abs(muX).^2;           % NxMxS single real
    if hasKnown                                     % deterministic known nodes
        muX = reshape(muX, N*M, S); muX(kidx,:) = repmat(kval,1,S); muX = reshape(muX,N,M,S);
        vX  = reshape(vX,  N*M, S); vX(kidx,:)  = 0;             vX  = reshape(vX, N,M,S);
    end

    % ---- (2) total mean / variance at each observation -------------
    Mtot = zeros(N, M, 'single'); Vtot = zeros(N, M, 'single');
    for s = 1:S
        Mtot = Mtot + gS(s)    * circshift(muX(:,:,s), [dk(s), dl(s)]);
        Vtot = Vtot + absg2(s) * circshift(vX(:,:,s),  [dk(s), dl(s)]);
    end
    Vtot = Vtot + single(sigma2);

    % ---- (3,4) edge log-likelihoods and their sum ------------------
    sumLogL = zeros(N, M, Q, 'single');
    for s = 1:S
        mu_excl  = circshift(Mtot, [-dk(s), -dl(s)]) - gS(s)*muX(:,:,s);
        var_excl = circshift(Vtot, [-dk(s), -dl(s)]) - absg2(s)*vX(:,:,s);
        var_excl = max(var_excl, floorv);
        res      = Yshift(:,:,s) - mu_excl;
        for q = 1:Q
            d  = res - single(Gsq(s,q));
            ll = -(real(d).^2 + imag(d).^2) ./ var_excl;
            logL(:,:,s,q)  = ll;
            sumLogL(:,:,q) = sumLogL(:,:,q) + ll;
        end
    end

    % ---- (5) variable -> observation messages (extrinsic) ----------
    for s = 1:S
        lm = sumLogL - reshape(logL(:,:,s,:), N, M, Q);
        lm = lm - max(lm, [], 3);                          % stabilise
        w  = exp(lm);
        w  = w ./ sum(w, 3);
        P(:,:,s,:) = (1-damp)*P(:,:,s,:) + damp*reshape(w, N, M, 1, Q);
    end
end

% ---- final decision ---------------------------------------------------
[~, qhat] = max(sumLogL, [], 3);          % NxM index of best symbol
Xhat = reshape(dataSyms(qhat), N, M);
Xhat(knownMask) = knownVal(knownMask);
end

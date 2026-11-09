function G = build_effective_dd(ch, p, Emax)
% BUILD_EFFECTIVE_DD  2D effective delay-Doppler channel response G.
%
% Implements the input-output relation [MP] eq. (13) (= Embedded eq. 3):
%   y[k,l] = sum_{k',l'} G[k',l'] * x[[k-k']_N, [l-l']_M] + v[k,l]
% where, for a path i with gain h_i, integer delay tap alpha_i and
% (fractional) Doppler tap kreal_i = beta_i + gamma_i,
%   G[k', alpha_i] += h'_i * alpha(q),   q = (beta_i - k') mod N,
%   h'_i  = h_i * exp(-j2*pi*kreal_i*alpha_i/(N*M)),
%   alpha(q) = (e^{j2pi(-q-gamma)} - 1) / (N(e^{j(2pi/N)(-q-gamma)} - 1)).  ([MP] eq. 13)
%
% With Emax = Inf (default) every Doppler index k' is filled, giving the
% EXACT effective response (used to apply the true channel).  With a finite
% Emax the inter-Doppler-interference spread is truncated to |q| <= Emax
% around each path's Doppler tap ([MP] truncates the IDI to E neighbours),
% giving the sparse response the detector works with.
%
% G is N x M; only delay columns occupied by a path are nonzero.

if nargin < 3, Emax = Inf; end
N = p.N; M = p.M;
G = zeros(N, M);
kp = (0:N-1).';                      % output Doppler indices k'

for i = 1:numel(ch.gain)
    li     = ch.ltap(i);             % integer delay tap alpha_i
    kreal  = ch.kreal(i);
    beta   = round(kreal);
    gamma  = kreal - beta;           % fractional Doppler in (-0.5, 0.5]
    hp     = ch.gain(i) * exp(-1i*2*pi*kreal*li/(N*M));

    if abs(gamma) < 1e-9
        % Integer Doppler -> single tap at k' = beta (no Doppler spread).
        kidx = mod(beta, N);
        G(kidx+1, li+1) = G(kidx+1, li+1) + hp;
    else
        q = mod(beta - kp, N);                 % q for each output index k'
        q = mod(q + N/2, N) - N/2;             % map to [-N/2, N/2)
        a = -q - gamma;
        col = hp * (exp(1i*2*pi*a) - 1) ./ (N*(exp(1i*(2*pi/N)*a) - 1));
        if isfinite(Emax)
            col(abs(q) > Emax) = 0;            % truncate IDI spread
        end
        G(:, li+1) = G(:, li+1) + col;
    end
end
end

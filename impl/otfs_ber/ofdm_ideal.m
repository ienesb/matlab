function [nbits, nerr] = ofdm_ideal(ch, p, SNRd_lin_vec)
% OFDM_IDEAL  BER of OFDM with IDEAL (known) channel over the same EVA
% delay-Doppler realization, used as the reference curve in Fig. 9.
%
% The OTFS frame is viewed as N OFDM symbols of M subcarriers ([MP] Sec.
% IV-B).  Within OFDM symbol n the received frequency-domain vector is
%   y_n = H_n x_n + v,        H_n = W H_t,n W^H   ([MP] eqs. 22-23)
% with H_n built exactly from the known path parameters.  Doppler makes
% H_n a (circulant) banded matrix with inter-carrier interference; the
% known-channel receiver removes it via MMSE equalisation.
%
% SNRd_lin_vec : vector of linear data SNRs (H_n is built once and reused).
% Returns nbits, nerr : row vectors (one entry per SNR_d).

N = p.N; M = p.M;
nS = numel(SNRd_lin_vec);
W = 2*p.k_nu;                  % ICI half-bandwidth kept (subcarriers)

hi = ch.gain; li = ch.ltap; kreal = ch.kreal;
phi = kreal/(N*M);             % per-sample Doppler phase [cycles]
P = numel(hi);

% Dirichlet ICI kernel Dir_i(Delta), Delta in [-W,W].
Deltas = (-W:W).';
DirK = zeros(numel(Deltas), P);
for i = 1:P
    if abs(kreal(i)) < 1e-9
        DirK(:,i) = (Deltas==0);
    else
        num = exp(1i*2*pi*kreal(i)/N) - 1;
        den = exp(1i*2*pi*(Deltas/M + phi(i))) - 1;
        DirK(:,i) = num ./ (M*den);
    end
end

a = (0:M-1).';
nd = numel(Deltas);
nbits = zeros(1,nS); nerr = zeros(1,nS);

for n = 0:N-1
    n_off = n*(M+p.cp_len) + p.cp_len;
    cn = hi .* exp(1i*2*pi*phi*n_off);     % h_i e^{j2pi phi_i n_off}

    % Build sparse banded circulant H_n (once per OFDM symbol).
    ii = zeros(M*nd,1); jj = zeros(M*nd,1); vv = zeros(M*nd,1);
    for d = 1:nd
        Delta = Deltas(d);
        diagv = zeros(M,1);
        for i = 1:P
            diagv = diagv + cn(i)*DirK(d,i)*exp(-1i*2*pi*(a+Delta)*li(i)/M);
        end
        rng = (d-1)*M + (1:M);
        ii(rng) = a+1;  jj(rng) = mod(a+Delta,M)+1;  vv(rng) = diagv;
    end
    Hn = sparse(ii, jj, vv, M, M);
    HhH = Hn'*Hn;  Hh = Hn';

    % One random data vector, reused across SNR_d (only noise/scale change).
    bits = randi([0 1], M, p.Qbits);
    for si = 1:nS
        Es  = SNRd_lin_vec(si);
        sym = qam4_mod(bits, p.alphabet*sqrt(Es), p.bitmap);
        y   = Hn*sym + sqrt(1/2)*(randn(M,1)+1i*randn(M,1));
        xhat = (HhH + (1/Es)*speye(M)) \ (Hh*y);
        bhat = qam4_demod(xhat, p.alphabet*sqrt(Es), p.bitmap);
        nbits(si) = nbits(si) + numel(bits);
        nerr(si)  = nerr(si)  + sum(bits(:) ~= bhat(:));
    end
end
end

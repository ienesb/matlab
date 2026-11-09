function Qh = liu_precoders(N, N1, N2, c1, c2)
% LIU_PRECODERS  Q^H = F_N U for the unitary bases of Liu et al. eqs. (8)-(10):
% the subcarrier symbols are X = Q^H s and the FIM depends on s only through
% |Q^H s|^2.
%   Qh = liu_precoders(N, N1, N2, c1, c2)   (defaults N1 = 32, N2 = 4, c1 = 1/16, c2 = 1/8)
%   Qh.OFDM = I            (U = F_N^H)       alpha_OFDM = 0
%   Qh.SC   = F_N          (U = I)           alpha_SC   = sqrt((N^2-1)/(12N^2))
%   Qh.OTFS = F_N (F_N1^H kron I_N2)         alpha_OTFS = sqrt((N2^2-1)/(12N2^2))
%   Qh.AFDM = F_N Lam_c1^H F_N^H Lam_c2^H,  Lam_c = diag(e^{-j2pi c n^2})  (2 N c1 integer)
% F_N = fft(eye(N))/sqrt(N) ([F]_{m,n} = e^{-j2pi(m-1)(n-1)/N}/sqrt(N), eq. 9).
% Same construction as impl/2605.14451v1/run_sim.m.  Qh.names lists the order.
    if nargin < 2 || isempty(N1), N1 = 32; end
    if nargin < 3 || isempty(N2), N2 = N/N1; end
    if nargin < 4 || isempty(c1), c1 = 1/16; end
    if nargin < 5 || isempty(c2), c2 = 1/8; end
    assert(N1*N2 == N, 'liu_precoders:N', 'N1*N2 must equal N');
    assert(abs(2*N*c1 - round(2*N*c1)) < 1e-12, 'liu_precoders:c1', '2 N c1 must be an integer');
    F_N  = fft(eye(N))/sqrt(N);
    F_N1 = fft(eye(N1))/sqrt(N1);
    n = (0:N-1).';
    Lam_c1 = diag(exp(-1j*2*pi*c1*n.^2));
    Lam_c2 = diag(exp(-1j*2*pi*c2*n.^2));
    Qh.OFDM = eye(N);
    Qh.SC   = F_N;
    Qh.OTFS = F_N * kron(F_N1', eye(N2));
    Qh.AFDM = F_N * Lam_c1' * F_N' * Lam_c2';
    Qh.names = {'SC', 'AFDM', 'OTFS', 'OFDM'};
    Qh.alpha_spread = struct('OFDM', 0, 'SC', sqrt((N^2-1)/(12*N^2)), ...
        'OTFS', sqrt((N2^2-1)/(12*N2^2)), 'AFDM', sqrt((1-4*c1^2)/12));
end

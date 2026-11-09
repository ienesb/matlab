function J = getFIM_OFDM_noICI(X, deltaf, Ts, sigma2, alpha, tau, nu)
% Single-target FIM for OFDM, no-ICI channel (Section 3.1.1 of ofdm_fim.tex).
% Parameter order: eta = [alphaR, alphaI, tau, nu].
%
% Inputs:
%   X       : N x M transmitted symbol grid (subcarrier x OFDM symbol)
%   deltaf  : subcarrier spacing [Hz]
%   Ts      : OFDM symbol duration (including CP) [s]
%   sigma2  : noise variance
%   alpha   : complex path gain (scalar)
%   tau     : path delay [s]
%   nu      : path Doppler [Hz]
%
% Output:
%   J       : 4 x 4 Fisher information matrix.

    [N, M] = size(X);

    ns = (0:N-1).';
    ms = (0:M-1);

    % Common phase factor: e^{-j2pi n Df tau} * e^{+j2pi m Ts nu}
    phase = exp(-1j*2*pi*deltaf*ns*tau) * exp(1j*2*pi*Ts*ms*nu);   % N x M

    % d mu / d alphaR
    g_aR = X .* phase;                                              % N x M
    % d mu / d alphaI = j * g_aR
    g_aI = 1j * g_aR;
    % d mu / d tau   = -j2pi alpha n Df * g_aR
    g_tau = (-1j*2*pi*deltaf) * (ns .* g_aR) * alpha;
    % d mu / d nu    = +j2pi alpha m Ts * g_aR
    g_nu  = (1j*2*pi*Ts) * (g_aR .* ms) * alpha;

    G = [g_aR(:), g_aI(:), g_tau(:), g_nu(:)];                      % NM x 4

    J = (2/sigma2) * real(G' * G);
end

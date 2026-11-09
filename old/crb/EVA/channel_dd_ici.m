function y_DD = channel_dd_ici(x_DD, alpha_gts, tau_gts, nu_gts, deltaf, Ts, sigma2)
% CHANNEL_DD_ICI  DD-domain OTFS channel under ICI (no ISI), rectangular pulses.
%
% Implements the last equation of otfs_io_equations.tex:
%
%   y[k,l] = (1/MN) sum_p alpha_p sum_{m,n,n',k',l'} x[k',l']
%            * exp( j2pi*( m*l'/M - n'*k'/N ) )      % ISFFT kernel
%            * exp(-j2pi*n'*deltaf*tau_p )           % delay phase
%            * exp( j2pi*m*Ts*nu_p )                 % Doppler phase
%            * I_p( n - n' ; nu_p )                  % ICI coefficient
%            * exp(-j2pi*( m*l/M - n*k/N ) )         % SFFT kernel
%            + z[k,l].
%
% Computed efficiently as ISFFT -> TF ICI channel (channel.m) -> SFFT.
% The dummy-index sign on the delay phase / ICI coefficient follows the
% codebase convention (otfs_io.tex Eq. 18, channel.m); the two forms are
% statistically equivalent.

    [N, M] = size(x_DD);

    % ISFFT: x[k,l] -> X[n,m]
    X_TF = isfft(x_DD);

    % TF-domain ICI channel (no AWGN here; injected in DD below)
    Y_TF = channel(X_TF, alpha_gts, tau_gts, nu_gts, deltaf, Ts, 0);

    % SFFT: Y[n,m] -> y[k,l]
    y_DD = sfft(Y_TF);

    % AWGN z[k,l]. SFFT is unitary so this matches TF-domain AWGN in distribution.
    if sigma2 > 0
        y_DD = y_DD + (randn(N, M) + 1j*randn(N, M)) * sqrt(sigma2/2);
    end
end

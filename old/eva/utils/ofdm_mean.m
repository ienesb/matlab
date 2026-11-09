function Y = ofdm_mean(X, deltaf, Ts, alpha, tau, nu, use_ici)
% OFDM_MEAN  Noise-free OFDM time-frequency input/output mean (ofdm_mcrb.tex
% eqs. (3)-(4)) -- the SAME parametric ICI channel used by the MCRB machinery
% (getAB_OFDM_core), so the BER pipeline and the (later) MCRB share one model.
%
%   Y[n,m] = sum_p alpha_p sum_{n'} X[n',m] e^{-j2pi n' deltaf tau_p}
%                                            e^{ j2pi m  Ts    nu_p } I_p(n'-n)
%
% with the subcarrier index n = 0..N-1 (delay axis, rows), the OFDM-symbol
% index m = 0..M-1 (Doppler axis, cols), and the ICI kernel
%   I_p(Delta) = get_ici_coeffs(Delta, nu_p, T),   T = 1/deltaf
% (unit peak at Delta = 0; use_ici = false collapses it to I = delta).
%
% Inputs:
%   X       : N x M transmitted TF symbol grid
%   deltaf  : subcarrier spacing [Hz]
%   Ts      : OFDM symbol duration (incl. CP) [s]
%   alpha   : P x 1 complex path gains
%   tau     : P x 1 path delays   [s]
%   nu      : P x 1 path Dopplers [Hz]
%   use_ici : true -> ICI kernel, false -> no-ICI (I = delta)
%
% Output:
%   Y       : N x M noise-free received TF grid.

    [N, M] = size(X);  P = numel(alpha);
    ns = (0:N-1).';  ms = (0:M-1);  Y = zeros(N, M);
    if use_ici
        T = 1/deltaf;  [NP, NN] = ndgrid(0:N-1, 0:N-1);  diff_n = NP - NN;
    end
    for p = 1:P
        phase = exp(-1j*2*pi*deltaf*ns*tau(p)) * exp(1j*2*pi*Ts*ms*nu(p));  % N x M
        A = X .* phase;
        if use_ici
            Y = Y + alpha(p) * (get_ici_coeffs(diff_n, nu(p), T).' * A);
        else
            Y = Y + alpha(p) * A;
        end
    end
end

function ch = gen_eva_channel(p)
% GEN_EVA_CHANNEL  Draw one random realization of the EVA delay-Doppler
% channel used in the Embedded paper (Sec. IV).
%
% Each of the 9 EVA taps gets:
%   - integer delay tap   l_i  = round(tau_i * M * deltaf)
%   - complex path gain   h_i  ~ CN(0, P_i)  with EVA power profile P_i
%   - a single Doppler    nu_i = nu_max*cos(theta_i), theta_i ~ U[-pi,pi]
%     (Jakes' model, Sec. IV), expressed as a (generally fractional) tap
%     k_real_i = nu_i * N * T.
%
% Returns struct ch with column vectors: gain, ltap (int delay), kreal
% (real-valued Doppler tap).

P = numel(p.eva_delays_ns);

% Integer delay taps (paper rounds delays to nearest tap; 1/(M*deltaf)
% resolution is sufficient, Sec. II-C).
ltap = round(p.eva_delays_ns*1e-9 * p.M * p.deltaf);   % 1xP
ltap = min(ltap, p.l_tau);                              % clip to l_tau

% Path gains: Rayleigh fading with the EVA average power profile,
% normalized to unit total average power.
lin_pow = 10.^(p.eva_pow_dB/10);
lin_pow = lin_pow/sum(lin_pow);
g = sqrt(lin_pow/2) .* (randn(1,P) + 1i*randn(1,P));    % CN(0,P_i)

% Doppler via Jakes' model, uniform angle of arrival over [-pi,pi].
theta = -pi + 2*pi*rand(1,P);
nu    = p.nu_max*cos(theta);                            % [Hz]
kreal = nu * p.N * p.T;                                 % real Doppler tap

ch.gain  = g(:);
ch.ltap  = ltap(:);
ch.kreal = kreal(:);
end

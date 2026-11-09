% diag_ofdm_floor.m
% Why does the OFDM BER floor at high SNR?  Decompose the floor into its two
% signal-proportional candidates by swapping one impairment at a time:
%
%   A  real CE  + ICI channel   : the curve seen in ber_ce_channels (baseline)
%   B  GENIE CE + ICI channel   : removes ALL channel-estimation error
%                                 -> whatever floor remains is pure ICI (the
%                                    1-tap equalizer cannot cancel it)
%   C  GENIE CE + no-ICI channel: removes ICI too -> should show NO floor
%                                 (pure AWGN reference)
%
% Genie H is the n'=n coefficient of the OFDM-ICI mean, i.e. the single complex
% gain the 1-tap MMSE equalizer is *trying* to estimate:
%   H[n,m] = sum_p alpha_p e^{-j2pi n df tau_p} e^{j2pi m Ts nu_p} I_p(0)
clear; clc;
thisdir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisdir, 'utils'));

SNR_dbs = [10 20 30 35];
SNRp_dB = 50;   nTrials = 20;   M_order = 4;   nbits = log2(M_order);
seed_base = 9100;

% ---- EVA grid/channel (same as ber_ce_channels 'eva') ----
N = 512;  M = 128;  deltaf = 15e3;  Ts = 1/deltaf;  T = 1/deltaf;
fc = 4e9;  c = 3e8;  v = 120/3.6;  nu_max = v*fc/c;
eva_delays_ns = [0 150 370];  eva_pow_db = [0 -1.4 -0.6];
delay_taps = round(eva_delays_ns*1e-9 * N * deltaf);
tau = (delay_taps/(N*deltaf)).';
amp = (10.^(eva_pow_db/20)).';  amp = amp/norm(amp);
P = numel(tau);  l_tau = max(delay_taps);  k_nu = ceil(nu_max*M/deltaf);

otfs_nondata = (2*l_tau+1)*M;
[dN_p, dM_p, n_pilot_cells] = match_pilot_spacing(N, M, otfs_nondata, l_tau, k_nu);
n_pilots = 1:dN_p:N;  m_pilots = 1:dM_p:M;
pilot_mask = false(N,M);  pilot_mask(n_pilots,m_pilots) = true;
data_mask  = ~pilot_mask;

fprintf('EVA: pilot spacing (%d,%d), %d pilots.  nu_max*T = %.4f -> ICI ~ %.1f dB\n', ...
        dN_p, dM_p, n_pilot_cells, nu_max*T, 10*log10((pi*nu_max*T)^2/3));
fprintf('  phase rotation between freq pilots (longest path): %.0f deg\n\n', ...
        360*dN_p*deltaf*max(tau));

variants = {'A real CE + ICI','B GENIE CE + ICI','C GENIE CE + no ICI'};
ber = zeros(numel(SNR_dbs), 3);

for si = 1:numel(SNR_dbs)
    N0  = 1/db2pow(SNR_dbs(si));
    a_p = sqrt(db2pow(SNRp_dB) * N0 / n_pilot_cells);
    for vi = 1:3
        use_ici_ch = (vi <= 2);          % C uses a no-ICI channel
        genie      = (vi >= 2);          % B,C use genie CE
        be = zeros(1,nTrials);
        parfor t = 1:nTrials
            rs = RandStream('Threefry','Seed',seed_base);  rs.Substream = 1000*si+t;
            RandStream.setGlobalStream(rs);
            nu    = (nu_max * cos(2*pi*rand(1,P) - pi)).';
            alpha = amp .* exp(1j*2*pi*rand(P,1));
            sym = randi(M_order, N, M) - 1;
            Xd  = qammod(sym, M_order, UnitAveragePower=true);
            Z   = sqrt(N0/2)*(randn(N,M) + 1j*randn(N,M));

            d = Xd;  d(pilot_mask) = a_p;
            Y = ofdm_mean(d, deltaf, Ts, alpha, tau, nu, use_ici_ch) + Z;

            if genie
                % n'=n coefficient: what the 1-tap equalizer is trying to estimate
                ns = (0:N-1).';  ms = (0:M-1);
                H_hat = zeros(N,M);
                for p = 1:P
                    ph = exp(-1j*2*pi*deltaf*ns*tau(p)) * exp(1j*2*pi*Ts*ms*nu(p));
                    if use_ici_ch
                        I0 = get_ici_coeffs(0, nu(p), T);   % I_p(0)
                    else
                        I0 = 1;
                    end
                    H_hat = H_hat + alpha(p) * ph * I0;
                end
            else
                H_ls = Y(n_pilots,m_pilots) ./ d(n_pilots,m_pilots);
                F_H  = griddedInterpolant({n_pilots,m_pilots}, H_ls, 'linear','nearest');
                [Nq,Mq] = ndgrid(1:N,1:M);
                H_hat = F_H(Nq,Mq);
            end
            Xsoft = conj(H_hat).*Y ./ (abs(H_hat).^2 + N0);
            sh = qamdemod(Xsoft(data_mask), M_order, UnitAveragePower=true);
            be(t) = mean(int2bit(sym(data_mask),nbits) ~= int2bit(sh,nbits), 'all');
        end
        ber(si,vi) = mean(be);
    end
    fprintf('SNR=%2d dB | A(real CE,ICI)=%.3e | B(genie,ICI)=%.3e | C(genie,noICI)=%.3e\n', ...
            SNR_dbs(si), ber(si,1), ber(si,2), ber(si,3));
end

fprintf('\nIf B stays high  -> the floor is ICI (1-tap equalizer cannot cancel it).\n');
fprintf('If B drops a lot -> the floor is CHANNEL-ESTIMATION (LS + linear interp).\n');

clear;
close all;
clc;

rcs2_dB = 1.5;
nIter = 1;
pilot_ratio = 0.05;
monteCarlo = 500;
is_genie = 0;
is_data_only = 0;

params = init_simulation_params(rcs2_dB, pilot_ratio, nIter, monteCarlo, is_genie, is_data_only);
% params.nus = [375.8984, 377.3253];
v2struct(params);

[X, pilot_mask] = generate_ofdm_symbols(params);
Xp = zeros(size(X));
Xp(pilot_mask) = X(pilot_mask);

H = generate_channel(params);
Y = generate_received_signal(X, H, params);


[Hhat, noiseVar, Hpilots] = estimateChannelFromPilots(Y, Xp, pilot_mask);

X_hat = data_demodulation(Y, Hhat, X, pilot_mask, params);
X_hat(pilot_mask) = X(pilot_mask);

getSer(X, X_hat, pilot_mask);






% x = ifft(X);
% 
% h = ifft(H);


% x = ifft(X(:, 1));
% 
% h = ifft(H(:, 1));
% 
% y1 = ifft(Y1(:, 1));
% y2 = cconv(x.', h.', N);
% y2 = y2.';
% 
% isequal(y2, y1)


% y = 0;
% for col_idx = 1:M
%     temp = cconv(x(:, 1).', h(:, 1).', N);
%     y = [y; temp.'];
% end
% y = y(2:end);
% Y2 = fft(reshape(y, [N, M]));
% 
% diff = sum(abs(Y2(:) - Y1(:)))
% sum(abs(Y2(:)) - abs(Y1(:)))





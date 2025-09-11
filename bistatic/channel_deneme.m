clear;
close all;
clc;

params = init_params();

v2struct(params);

X = generate_ofdm_symbols(params);

H = generate_channel(params);
Y1 = X .* H;

x = ifft(X);

x = x(:);




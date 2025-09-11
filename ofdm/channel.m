x = randn(10000, 1);

channel_delays = [0, 3, 5];

channel_gains = [1 0.6 0.4];

y1 = zeros(length(x)+max(channel_delays), 1);

for idx = 1:length(channel_delays)
    delay = channel_delays(idx);
    gain = channel_gains(idx);
    y1(1:length(x)+delay) = y1(1:length(x)+delay) + [zeros(delay, 1); gain*x];
end

h = [1 0 0 0.6 0 0.4];

y2 = conv(x, h);

% figure; hold on;
% stem(x);
% stem(y1);
% grid on;

y3 = cconv(x, h, 10005);

X = fft(x, 10005);
H = fft(h, 10005);
H = H.';
Y = X .* H;
y4 = ifft(Y);

isequal(y3, y4)

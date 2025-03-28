m = 1:10000;

r = floor(log2(3.^m)) ./ m;

b = ones(size(m)) * log2(3);

figure; hold on;
plot(m, r);
plot(m, b);
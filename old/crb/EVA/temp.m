clear;
close all;


v = [(1 + 1j), (-1 + 1j), (1 - 1j), (-1 -1j)] ./ sqrt(2);
v = v.';

t = conj(v) .* (v.');
t = t(:);

rhos = 0:0.001:1;
results = zeros(size(rhos));
for rho_idx = 1:length(rhos)
    % rho = 1e-2;
    rho = rhos(rho_idx);
    
    e0 = (1-rho)^2;
    e1 = rho*(1-rho);
    e2 = rho^2;
    p = [e0; e1; e1; e2; e1; e0; e2; e1; e1; e2; e0; e1; e2; e1; e1; e0];
    
    results(rho_idx) = sum(t .* p * 0.25);
end

figure;
plot(rhos, results);

% OTFS

N = 10;
M = 10;
results = zeros(N, M);


for n = 0:(N-1)
    for m = 0:(M-1)
        s = 0;
        for k = 0:(N-1)
            for l = 0:(M-1)
                for kp = 0:(N-1)
                    for lp = 0:(M-1)
                        s = s + exp(1j*2*pi*(m*(lp - l)/M - n*(kp - k)/N));
                    end
                end
            end
        end
        results(n+1, m+1) = abs(s);
    end
end



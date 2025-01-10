function y = twconv(h, x)

if any(size(h)~=size(x))
  error('Input matrices must be same size.');
end

M = size(x, 1);
N = size(x, 2);

y = zeros(M, N);

for k=0:M-1
    for l=0:N-1
        for kp=0:M-1
            for lp=0:N-1
                if k-kp >= 0 
                    y(k+1,l+1) = y(k+1,l+1) + h(kp+1,lp+1)*x(mod(k-kp,M)+1, mod(l-lp,N)+1)*exp(1j*2*pi*(k-kp)*lp/M/N);
                else 
                    y(k+1,l+1) = y(k+1,l+1) + h(kp+1,lp+1)*x(mod(k-kp,M)+1, mod(l-lp,N)+1)*exp(1j*2*pi*(k-kp)*lp/M/N)*exp(-1j*2*pi*mod(l-lp,N)/N);
                end
            end
        end
    end
end
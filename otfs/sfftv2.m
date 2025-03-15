function Xdd = sfftv2(Xtf, M, N)
    Xdd = zeros(M,N);
    for k=0:N-1
        for l=0:M-1
            for n=0:N-1
                for m=0:M-1
                    Xdd(l+1,k+1) = Xdd(l+1,k+1) + Xtf(m+1, n+1)*exp(1j*2*pi*(m*l/M - n*k/N));
                end
            end
        end
    end
    % Xdd = Xdd ./ (M*N);
end
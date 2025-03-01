function Xtf = isfftv2(Xdd, M, N)
    Xtf = zeros(M,N);
    for n=0:N-1
        for m=0:M-1
            for k=0:N-1
                for l=0:M-1
                    Xtf(m+1,n+1) = Xtf(m+1,n+1) + Xdd(l+1, k+1)*exp(1j*2*pi*(n*k/N - m*l/M));
                end
            end
        end
    end
    Xtf = Xtf ./ (M*N);
end
function c = dd_crosstalk_coefficients(k, kp, l, lp, nu, tau, T, N, M)    
    
    df = 1/T;
    c = 0;
    for n = 0:N-1
        for np = 0:N-1
            for m = 0:M-1
                for mp = 0:M-1
                    c = c + cross_ambg((n-np)*T-tau, (m-mp)*df-nu, T)/(N*M) * exp(1j*2*pi*(np * T * nu - mp * df * tau + ( (np * kp) / N - (mp * lp) / M ) - ( (n * k) / N - (m * l) / M )));
                end
            end
        end
    end

end
function C = cross_ambg(tau, nu, T)
    if tau < -T || tau > T
        C = 0;
    elseif tau < 0
        if nu == 0
            C = T + tau;
        else
            C = 1j*(exp(-1j*2*pi*nu*(T+tau))-1)/(2*pi*nu);
        end
    else
        if nu == 0
            C = T - tau;
        else
            C = 1j*(exp(-1j*2*pi*nu*T)-exp(-1j*2*pi*nu*tau))/(2*pi*nu);
        end
    end
    C = C / T;
end


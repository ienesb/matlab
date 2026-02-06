function CAF = getCAF(tau, nu, T) % T
    if abs(tau) < T
        if tau >= 0
            op1 = T;
            op2 = tau;
        else
            op1 = T + tau;
            op2 = 0;
        end
        if nu ~= 0
            CAF = 1j * (exp(-1j*2*pi*nu*op1) - exp(-1j*2*pi*nu*op2)) / (2*pi*nu);
        else
            CAF = op1 - op2;
        end
        CAF = CAF / T;
    else
        CAF = 0;
    end

end
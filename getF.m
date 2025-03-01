function F = getF(Na, Ns)
    G = 90;
    
    theta_min = -pi/4;
    theta_max = pi/4;
    
    n = (1:Na)';
    phis = linspace(theta_min, theta_max, G);
    
    A = exp(1j.*((n-1)*phis));
    
    b = ones(G, 1);
    
    f_unc = (A * A') \ (A * b);
    f_opt = f_unc / sqrt(f_unc' * (A * A') * f_unc);
    
    F = f_opt;
end
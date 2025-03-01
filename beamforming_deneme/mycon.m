function [c,ceq] = mycon(x)
    Na = 64;
    G = 36;
    n = (1:Na)';
    phis = linspace(-pi/2, pi/2, G);

    A = exp(1j.*((n-1)*phis));

    ceq = norm(A'*x)^2-1;
    c = -1;
end
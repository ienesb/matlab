function J = cost_tau_nu(eta, deltaf, T, X_use, y)
% eta = [nu; tau]
    nu  = eta(1);
    tau = eta(2);

    N = size(X_use,1);
    M = size(X_use,2);
    ns = (0:(N-1)).';
    ms = 0:(M-1);

    C = (exp(-1j*2*pi*deltaf*ns*tau) * exp(1j*2*pi*T*ms*nu)) .* X_use;
    c = C(:);

    cy = c' * y;
    cc = c' * c;

    % concentrated NLL up to constant: ||y||^2 - |c^H y|^2/(c^H c)
    J = real( (y' * y) - (abs(cy)^2 / cc) );
end
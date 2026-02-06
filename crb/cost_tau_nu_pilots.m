function J = cost_tau_nu_pilots(eta, deltaf, T, X, y, pilot_indices)
    nu  = eta(1);
    tau = eta(2);

    [N,M] = size(X);
    ns = (0:(N-1)).';
    ms = 0:(M-1);

    P = pilot_indices;               % logical N×M
    Y = reshape(y, N, M);

    C = (exp(-1j*2*pi*deltaf*ns*tau) * exp(1j*2*pi*T*ms*nu)) .* X;

    % only pilots
    c = C(P);
    yP = Y(P);

    cy = c' * yP;
    cc = c' * c;

    J = real( (yP' * yP) - (abs(cy)^2 / cc) );
end
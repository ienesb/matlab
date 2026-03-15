function Y_TF = channel(X_TF, alpha_gts, tau_gts, nu_gts, deltaf, Ts, sigma2)
    T = 1 / deltaf;
    [N, M] = size(X_TF);
    ns = (0:(N-1)).';
    ms = 0:(M-1);    

    %% ICI case
    Y_TF = zeros(N, M);
    P = length(tau_gts);
    for p = 1:P
        tau_gt = tau_gts(p);
        nu_gt = nu_gts(p);
        alpha_gt = alpha_gts(p);
        for n = 0:(N-1)
            I = get_ici_coeffs(n-ns, nu_gt, T);
    
            temp = (exp(-1j*2*pi*deltaf*ns*tau_gt) * exp(1j*2*pi*Ts*ms*nu_gt));
            temp = X_TF .* temp .* I;
            temp = sum(temp, 1);
            Y_TF(n+1, :) = Y_TF(n+1, :) + alpha_gt * temp;
        end
    end
    %% AWGN
    Z = (randn(N, M) + 1j*randn(N, M)) * sqrt(sigma2/2);
    Y_TF = Y_TF + Z;
end
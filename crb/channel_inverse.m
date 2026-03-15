function Y_TF = channel_inverse(X_TF, alpha_gts, tau_gts, nu_gts, deltaf, Ts, N0)
    T = 1 / deltaf;
    [N, M] = size(X_TF);
    ns = (0:(N-1)).';
    ms = 0:(M-1);    
    
    P = length(alpha_gts);

    %% ICI case
    Y_TF = zeros(N, M);
    for m = 0:(M-1)
        Hb = zeros(N, N);
        for p = 1:P
            tau_gt = tau_gts(p);
            nu_gt = nu_gts(p);
            alpha_gt = alpha_gts(p);
            Hb = Hb + getChannelSubMatrix(N, m ,tau_gt, nu_gt, alpha_gt, deltaf, Ts);
        end
        % Y_TF(:, (m+1)) = pinv(Hb) * X_TF(:, (m+1)); 
        
        % MMSE Equalization instead of pinv (Zero-Forcing)
        Hb_H = Hb';
        Y_TF(:, (m+1)) = (Hb_H * Hb + N0 * eye(N)) \ (Hb_H * X_TF(:, (m+1)));
    end
end
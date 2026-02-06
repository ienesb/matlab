function [detectionReport, X_hat] = gradient_update(X_hat, Y, detectionReport, c, lambda, params_OFDM) % X_hat, H_hat, Y, detectionReport (tau_hats, nu_hats, alpha_hats), params_OFDM

N = 400;
M = 60;

% detectionReport(end+1, 1) = 59.9298 + 1.85;
% detectionReport(end+1, 2) = 2.0214 + 1.4;
% detectionReport(end+1, 3) = 0;
% detectionReport(end+1, 4) = 5e-8;

detectionReport(:, 1) = detectionReport(:, 1) * 2 / c;
detectionReport(:, 2) = detectionReport(:, 2) * 2 / lambda;
% detectionReport(end+1, 1) = 0.3995e-6 + 1.9e-8;
% detectionReport(end+1, 2) = 377.3253;
% detectionReport(end+1, 3) = 0;
% detectionReport(end+1, 4) = 5e-8;

nDet = size(detectionReport, 1);

H_hat = zeros(N, M);
for ii = 1:nDet
    tau_hat = detectionReport(ii, 1);
    nu_hat = detectionReport(ii, 2);
    alpha_hat = detectionReport(ii, 4);

    b_tau = freqSteeringOFDM(tau_hat, params_OFDM);
    c_nu = conj(timeSteeringOFDM(nu_hat, params_OFDM));

    H_hat = H_hat + alpha_hat * (b_tau * c_nu.');
end

detectionReport_new = detectionReport;

for e = 1:100
    for ii = 1:nDet
        tau_hat = detectionReport(ii, 1);
        nu_hat = detectionReport(ii, 2);
        alpha_hat = detectionReport(ii, 4);
    
        [tau_hat_new, nu_hat_new, X_hat_new] = gradient_step(tau_hat, nu_hat, alpha_hat, X_hat, Y);
        detectionReport_new(ii, 1) = tau_hat_new;
        detectionReport_new(ii, 2) = nu_hat_new;
    end
    detectionReport = detectionReport_new; 
    X_hat = X_hat_new;
end

detectionReport(:, 1) = detectionReport(:, 1) * c / 2;
detectionReport(:, 2) = detectionReport(:, 2) * lambda / 2;


end
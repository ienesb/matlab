function H_hat = remove_target(H_hat, data_mask, tau_hat, nu_hat, alpha_hat, params)
    v2struct(params);

    b = getb(tau_hat, params);
    c = getc(nu_hat, params);

    temp = alpha_hat * b * c.';
    temp(data_mask) = 0;
    
    H_hat = H_hat - temp;
end
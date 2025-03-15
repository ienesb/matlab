Na = 64;
Nrf = 4;
B = 6;
G = Na;

theta_min = -pi/4;
theta_max = pi/4;
fov = theta_max - theta_min;

ns = (1:Na)';
phis = linspace(theta_min, theta_max, G);


% b = 6;
% n = 2;

for b = 1:B

    % figure; hold on;
    
    for n = 1:Nrf
    
        t_min = theta_min + (b-1)*fov/B + (n-1)*fov/B/Nrf;
        t_max = theta_min + (b-1)*fov/B + n*fov/B/Nrf;
        
        
        [~, idx1] = min(abs(phis - t_min));
        [~, idx2] = min(abs(phis - t_max));
        
        idx1 = max(idx1 - 1, 1);
        idx2 = min(idx2 + 1, G);
        
        A = exp(1j.*((ns-1)*phis));
        
        bv = zeros(G, 1);
        bv(idx1:idx2) = 1;
        
        % Given: A (matrix), b (vector)
        u_unc = pinv(A) * bv; % Unconstrained solution
        
        % Normalize to satisfy f^H A A^H f = 1
        u_opt = u_unc ./ sqrt((u_unc' * (A * A') * u_unc));
        
        disp('****************')
        abs(u_opt' * (A * A') * u_opt)
        % (u_unc' * (A * A') * u_unc)
        disp('****************')
        
        norm(A'*u_opt-bv)^2;

        
        % plot(phis*180/pi, abs(A'*u_opt))
    
    end
end
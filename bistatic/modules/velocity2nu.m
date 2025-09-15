function nu = velocity2nu(vel_target, pos_target, pos_tx, pos_rx, lambda)
    u_tx = (pos_target - pos_tx) / norm(pos_target - pos_tx);
    u_rx = (pos_rx - pos_target) / norm(pos_rx - pos_target);
    
    nu = vel_target.' * (u_tx - u_rx) / lambda;

    % u_tx = (pos_target - pos_tx); 
    % u_rx = -(pos_target - pos_rx);
    % 
    % u_tx = u_tx / norm(u_tx);
    % u_rx = u_rx / norm(u_rx);
    % 
    % nu = 2 * dot(vel_target, (u_tx + u_rx)) / lambda;
end
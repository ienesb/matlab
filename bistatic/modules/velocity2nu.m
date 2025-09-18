function nu = velocity2nu(vel_target, pos_target, pos_tx, pos_rx, lambda)
    u_tx = (pos_tx - pos_target) / norm(pos_tx - pos_target);
    u_rx = (pos_rx - pos_target) / norm(pos_rx - pos_target);
    
    % nu = vel_target.' * (u_tx + u_rx) / lambda;
    nu = vel_target.' * (-u_tx + u_rx)*2 / lambda;
end
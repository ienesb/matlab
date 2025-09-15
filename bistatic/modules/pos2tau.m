function tau = pos2tau(pos_target, pos_tx, pos_rx)
    d1 = norm(pos_target - pos_tx);
    d2 = norm(pos_target - pos_rx);
    
    % c = physconst('LightSpeed');
    c = 3e8;
    d = d1 + d2;
    tau = d./c;
end
function tau = pos2tau(pos_target, pos_tx, pos_rx)
    d1 = norm(pos_target - pos_tx);
    d2 = norm(pos_target - pos_rx);
    
    c = physconst('LightSpeed');

    tau = (d1 + d2)./c;
end
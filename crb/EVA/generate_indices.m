if type == "OFDM"
    % pilot_periods = [10, 2, 2, 1; 5, 5, 1, 1];
    np = 2;
    mp = 5;
    
    pilot_indices = zeros(N, M);
    pilot_indices(1:np:end, 1:mp:end) = 1;
    % pilot_indices(:, 1:10:41) = 1;
    pilot_indices = logical(pilot_indices);
    
    data_indices = ~pilot_indices;
    
else % OTFS
    pilot_indices = zeros(N, M);
    pilot_indices(35, 25) = 1;
    pilot_indices = logical(pilot_indices);
    
    guard_indices = zeros(N, M);
    guard_indices(25:45, 15:35) = 1;
    guard_indices(35, 25) = 0;
    guard_indices = logical(guard_indices);
    
    Xp = sqrt(sum(guard_indices(:)) + 1);
    
    data_indices = guard_indices + pilot_indices;
    data_indices = ~data_indices;
end
function C = my_kron_nd_fast(A, B)
    % Get sizes of A and B
    [m, n, l] = size(A);
    [p, q] = size(B);
    
    % Reshape A to (m, 1, n, 1, l) for broadcasting
    A_reshaped = reshape(A, m, 1, n, 1, l);
    
    % Reshape B to (1, p, 1, q) for broadcasting
    B_reshaped = reshape(B, 1, p, 1, q);
    
    % Perform element-wise multiplication using implicit broadcasting
    C = A_reshaped .* B_reshaped;
    
    % Reshape to final size (mp, nq, l)
    C = reshape(C, m * p, n * q, l);
end

function C = my_kron_nd(A, B)
    % Get the sizes of A and B
    [m, n, l, L] = size(A);
    [p, q, P, r1, r2] = size(B);
    
    % Ensure A and B have the expected third dimension (singleton)
    if L ~= 1 || P ~= 1
        error('A and B must have singleton third dimensions for broadcasting.');
    end
    
    % Expand A for broadcasting
    A_expanded = repelem(A, p, q, 1, 1, 1); % Expand A in first two dims

    % Expand B to match (mp, nq, l, r1, r2)
    B_tiled = repmat(B, m, n, l, 1, 1); % Tile B in first 3 dims
    
    % Compute Kronecker product
    C = A_expanded .* B_tiled;
end

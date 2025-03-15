function C = my_kron(A, B)
    % Get the sizes of matrices A and B
    [m, n] = size(A);
    [p, q] = size(B);
    
    % Expand A to match the Kronecker pattern
    A_expanded = repelem(A, p, q);
    
    % Create a large B matrix repeated across dimensions
    B_tiled = repmat(B, m, n);
    
    % Compute Kronecker product
    C = A_expanded .* B_tiled;
end
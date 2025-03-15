function C = my_kron_nd(A, B)
    % Get the sizes of A and B
    [m, n, l] = size(A);
    [p, q] = size(B);
    
    A_expanded = repelem(A, p, q, 1);
    
    B_tiled = repmat(B, m, n, l);
    
    C = A_expanded .* B_tiled;
end

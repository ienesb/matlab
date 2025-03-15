function result = kron_product_mult(A, B, x)
    % m = size(A, 1);
    n = size(A, 2);
    % p = size(B, 1);
    q = size(B, 2);

    X = reshape(x, [q, n]);

    Y = B * X * A.';
    result = Y(:);
end
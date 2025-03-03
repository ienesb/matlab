A = rand(2,2,3,1);       % A has size (2,2,3,1)
B = rand(2,2,1,4,5);     % B has size (2,2,1,4,5)

C2 = my_kron_nd(A, B);    % Output should have size (4,4,3,4,5)
disp(size(C2));           % Should print: [4 4 3 4 5]

% Get the sizes of A and B
[m, n, l, L] = size(A);
[p, q, P, r1, r2] = size(B);

% Ensure A and B have singleton third dimensions
if L ~= 1 || P ~= 1
    error('A and B must have singleton third dimensions for broadcasting.');
end

% Initialize the output matrix
C1 = zeros(m * p, n * q, l, r1, r2);

% Compute Kronecker product for each (l, r1, r2) slice
for i = 1:l
    for j = 1:r1
        for k = 1:r2
            for a = 1:m
                for b = 1:n
                    % Extract scalar from A
                    A_elem = A(a, b, i, 1);

                    % Multiply with full B(:,:,1,j,k)
                    C1((a-1)*p+1:a*p, (b-1)*q+1:b*q, i, j, k) = A_elem * B(:,:,1,j,k);
                end
            end
        end
    end
end

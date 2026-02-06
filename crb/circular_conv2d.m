function C = circular_conv2d(A, B)
    % % Get dimensions (assuming both are N x M)
    % [N, M] = size(A);
    % 
    % % Initialize output matrix
    % C = zeros(N, M);
    % 
    % % Loop only over the elements of A (the kernel)
    % % We shift matrix B for every element in A
    % for i = 1:N
    %     for j = 1:M
    %         if A(i, j) ~= 0  % Optimization: skip zeros
    %             % circshift(B, [row_shift, col_shift])
    %             % Shift B so that B(1,1) aligns with A(i,j) at the target output
    %             C = C + A(i, j) * circshift(B, [i-1, j-1]);
    %         end
    %     end
    % end
    % C = C / sqrt(N*M);

    a = isfft(A);
    b = isfft(B);
    c = a .* b;
    C = sfft(c);
end
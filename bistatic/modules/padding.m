function x_padded = padding(x, p)
    x_padded = [x(end-p+1:end, :) ; x ; x(1:p, :)];
    x_padded = [x_padded(:, end-p+1:end) x_padded x_padded(:, 1:p)];
end
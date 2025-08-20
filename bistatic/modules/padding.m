function x_padded = padding(x, p1, p2)
    x_padded = [x(end-p1+1:end, :) ; x ; x(1:p1, :)];
    x_padded = [x_padded(:, end-p2+1:end) x_padded x_padded(:, 1:p2)];
end
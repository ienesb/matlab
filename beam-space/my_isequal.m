function e = my_isequal(x, y)
    th = 10^-5;
    diff = sum(sum(sum(sum(sum(sum(abs(x - y)))))))
    e = diff < th;
end
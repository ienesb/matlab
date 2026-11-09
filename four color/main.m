c = combinations([0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1]);

c = table2array(c);

num_of_planars5 = 0;
for idx5 = 1:size(c, 1)
    A = zeros(5, 5);
    A(1, 2:5) = c(idx5, 1:4);
    A(2, 3:5) = c(idx5, 5:7);
    A(3, 4:5) = c(idx5, 8:9);
    A(4, 5) = c(idx5, 10);
    A = A + A.';
    if isPlanar(A)
        num_of_planars5 = num_of_planars5 + 1;
    end
end




c = combinations([0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1]);

c = table2array(c);

num_of_planars6 = 0;
for idx6 = 1:size(c, 1)
    A = zeros(6, 6);
    A(1, 2:6) = c(idx6, 1:5);
    A(2, 3:6) = c(idx6, 6:9);
    A(3, 4:6) = c(idx6, 10:12);
    A(4, 5:6) = c(idx6, 13:14);
    A(4, 6) = c(idx6, 15);
    A = A + A.';
    if isPlanar(A)
        num_of_planars6 = num_of_planars6 + 1;
    end
end



c = combinations([0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1]);

c = table2array(c);

num_of_planars7 = 0;
for idx7 = 1:size(c, 1)
    A = zeros(7, 7);
    A(1, 2:7) = c(idx7, 1:6);
    A(2, 3:7) = c(idx7, 7:11);
    A(3, 4:7) = c(idx7, 12:15);
    A(4, 5:7) = c(idx7, 16:18);
    A(4, 6:7) = c(idx7, 19:20);
    A(4, 7) = c(idx7, 21);
    A = A + A.';
    if isPlanar(A)
        num_of_planars7 = num_of_planars7 + 1;
    end
end
function piEst = computePi(m,n)
pointsInCircle = 0;

for i = 1:n
    % Generate random points.
    x = rand(m,1);
    y = rand(m,1);

    % Determine whether the points lie inside the unit circle.
    r = x.^2 + y.^2;
    pointsInCircle = pointsInCircle + sum(r<=1);
end

piEst = 4/(m*n) * pointsInCircle;
end
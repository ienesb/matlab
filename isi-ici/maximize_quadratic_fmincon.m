function b_opt = maximize_quadratic_fmincon(D)
    n = size(D, 1);
    
    % Initial guess (random unit vector)
    b0 = randn(n,1);
    b0 = b0 / norm(b0);  

    % Objective function: maximize b^T D b (converted to minimization)
    objFun = @(b) -b' * D * b;

    % Constraint: norm(b)^2 <= 1
    nonlcon = @(b) deal([], b' * b - 1);

    % Solve using fmincon
    options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter');
    b_opt = fmincon(objFun, b0, [], [], [], [], [], [], nonlcon, options);
end
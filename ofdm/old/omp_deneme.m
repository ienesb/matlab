A = [1 1/2 -1/sqrt(2); 0 sqrt(3)/2 -1/sqrt(2)];

v = [1; 1/2];

r = [0; 0];
errors = zeros(1, 10);
for k = 1:10
    inner_product = (v-r).'*A;
    [~, idx] = max(abs(inner_product));
    r = r + A(:, idx)*inner_product(idx);
    % figure; hold on;
    % quiver(0, 0, v(1), v(2))
    % quiver(0, 0, r(1), r(2))
    % legend("v", "r")
    errors(k) = norm(v - r)^2;
end

figure
plot(errors);
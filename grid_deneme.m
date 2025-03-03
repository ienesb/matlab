x1 = 1:3;
x2 = 1:5;
x3 = 1:4;

result = 0;
for i1 = x1
    for i2 = x2
        for i3 = x3
            result = result + i1^2 + i2^2 + i3^2;
        end
    end
end
result;

[X1v1, X2v1, X3v1] = meshgrid(x1, x2, x3);
[X1v2, X2v2, X3v2] = ndgrid(x1, x2, x3);

result2 = sum(sum(sum(X1v1.^2 + X2v1.^2 + X3v1.^2)));

U = eye(4, 64);

N = eye(64,1);

phis = linspace(-pi/4, pi/4, 50);



av2 = array_responsev2(phis, 64);

av1 = zeros(64, 50);

for idx = 1:length(phis)
    phi = phis(idx);
    av1(:, idx) = array_response(phi, 64);
end

av1 = av1.';

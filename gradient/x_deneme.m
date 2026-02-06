clear;
close all;
clc;

sigma2 = 0.1;
z = sqrt(sigma2/2) * (randn + 1j * randn);
h = randn;

x = (1 + 1j) / sqrt(2);
y = x*h + z;

xr_hat = -3:0.01:5;
xc_hat = -3:0.01:5;

J = zeros(length(xr_hat), length(xc_hat));

for xr_idx = 1:length(xr_hat)
    for xc_idx = 1:length(xc_hat)
        xr = xr_hat(xr_idx);
        xc = xc_hat(xc_idx);
        x_hat = xr + 1j*xc;
        y_hat = x_hat*h;
        J(xr_idx, xc_idx) = abs(y - y_hat).^2;
    end
end

[xc_grid, xr_grid] = meshgrid(xc_hat, xr_hat);

figure;
mesh(xc_grid, xr_grid, -J);
xlabel('X_c');
ylabel('X_r');
zlabel('J');
grid on;

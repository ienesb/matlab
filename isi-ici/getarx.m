function arx = getarx(theta, NR)
    ns = 1:NR;
    [theta_grid, ns_grid] = ndgrid(theta, ns);
    arx = exp(1j*pi*NR*sin(theta_grid).*(ns_grid-1));
    arx = arx.';
end
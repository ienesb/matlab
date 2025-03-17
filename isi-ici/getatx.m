function atx = getatx(theta, NT)
    ns = 1:NT;
    [theta_grid, ns_grid] = ndgrid(theta, ns);
    atx = exp(1j*pi*sin(theta_grid).*(ns_grid-1));
    atx = atx.';
end
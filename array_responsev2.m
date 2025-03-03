function a = array_responsev2(phi, Na)
    ns = (1:Na)';
    [phi_grid, ns_grid] = ndgrid(phi, ns);
    a = exp(1j*2*pi*sin(phi_grid).*(ns_grid-1));
    a = a.';
end
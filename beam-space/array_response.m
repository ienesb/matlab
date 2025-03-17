function a = array_response(phi, Na)
    ns = 1:Na;
    [phi_grid, ns_grid] = ndgrid(phi, ns);
    a = exp(1j*pi*sin(phi_grid).*(ns_grid-1));
    a = a.';
end
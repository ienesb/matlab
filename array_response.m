function a = array_response(phi, Na)
    ns = (1:Na)';
    a = exp(1j*2*pi*sin(phi)*(ns-1));
end
function dd = dd_crosstalk_coefficientsv2(nu, tau, T, N, M)    
    
    df = 1/T;

    ns = (0:(N-1))';
    nps = 0:(N-1);
    
    ms = (0:(M-1))';
    mps = 0:(M-1);
    
    arg1 = (ns - nps)*T-tau;
    arg1 = arg1(:);
    
    arg2 = (ms - mps)*df-nu;
    arg2 = arg2(:);
    arg2 = arg2.';
    
    C = cross_ambgv2(arg1, arg2, T);
    
    C = reshape(C, [N, N, M, M]);
    
    exp1 = exp(1j*2*pi*nps*T*nu);
    exp1 = reshape(exp1, [1, N, 1, 1]);
    C = C.*exp1;
    
    exp2 = exp(-1j*2*pi*mps*df*tau);
    exp2 = reshape(exp2, [1, 1, 1, M]);
    C = C.*exp2;
    
    dd = C;
    dd = fft(dd, N, 1);
    dd = ifft(dd, N, 2);
    dd = ifft(dd, M, 3);
    dd = fft(dd, M, 4);

end
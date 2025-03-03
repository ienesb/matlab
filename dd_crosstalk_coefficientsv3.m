function ddm = dd_crosstalk_coefficientsv3(nu, tau, T, N, M)    
    df = 1/T;
    Nnu = length(nu);
    Ntau = length(tau);

    ns = (0:(N-1))';
    nps = 0:(N-1);
    
    ms = (0:(M-1))';
    mps = 0:(M-1);
    
    arg1 = (ns - nps)*T;
    arg1 = arg1(:);
    arg1 = arg1 - tau;
    arg1 = arg1(:);

    arg2 = (ms - mps)*df;
    arg2 = arg2(:);
    arg2 = arg2 - nu;
    arg2 = arg2(:);
    arg2 = arg2.';

    C = cross_ambgv2(arg1, arg2, T);
    
    C = reshape(C, [N, N, Ntau, M, M, Nnu]);
    C = permute(C, [1,2,4,5,3,6]);
    
    nps = nps.';
    exp1 = exp(1j*2*pi*T*nps*nu);
    exp1 = permute(exp1, [3,1,4,5,6,2]);
    C = C.*exp1;

    mps = mps.';
    exp2 = exp(-1j*2*pi*df*mps*tau);
    exp2 = permute(exp2, [3,4,5,1,2,6]);
    C = C.*exp2;

    dd = C;
    dd = fft(dd, N, 1);
    dd = ifft(dd, N, 2);
    dd = ifft(dd, M, 3);
    dd = fft(dd, M, 4);


    ddm = permute(dd, [3, 1, 4, 2, 5, 6]);
    ddm = reshape(ddm, [M*N, M*N, Ntau, Nnu]);

end
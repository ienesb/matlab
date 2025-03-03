N = 64;
M = 64;
df = 15e3;
T = 1/df;

ks = 0:3;
ls = 0:3;

nus = ks ./ (N*T);
taus = ls ./ (M*df);

Nnu = 4;
Ntau = 4;
tic
ddm = zeros(N*M, N*M, Ntau, Nnu);
for idxnu = 1:Nnu
    for idxtau = 1:Ntau
        tau = taus(idxtau);
        nu = nus(idxnu);
        dd = dd_crosstalk_coefficientsv2(nu, tau, T, N, M);
        ddm(:,:,idxtau, idxnu) = dd2m(dd,M,N);  
    end
end

toc
tic
ddmv2 = dd_crosstalk_coefficientsv3(nus, taus, T, N, M);
toc















% ns = (0:(N-1))';
% nps = 0:(N-1);
% 
% arg1 = (ns - nps)*T;
% arg1 = arg1(:);
% arg1 = arg1 - tau;
% arg1 = arg1(:);
% 
% [nsg, npsg, taug] = ndgrid(ns, nps, tau);
% 
% arg1v2 = (nsg - npsg) .*T - taug;
% arg1v2 = arg1v2(:);

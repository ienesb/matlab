M = 64;          % number of subcarriers
N = 30;          % number of subsymbols/frame

Xdd = zeros(M,N);
Xdata = randi([0,1],2*M,N);
Xdd(1:M,:) = pskmod(Xdata,4,pi/4,InputType="bit");


Xtfv1 = isfft(Xdd, M);
Xtfv2 = isfftv2(Xdd, M, N);

difftf = Xtfv1 - Xtfv2;

figure;
surf(abs(Xtfv1))

figure;
surf(abs(Xtfv2))

figure;
surf(abs(difftf)) % sonuclar ayni


Xddv1 = sfft(Xtfv1, M);
Xddv2 = sfftv2(Xtfv1, M, N);

diffdd = Xddv1 - Xddv2;

figure;
surf(abs(Xdd))

figure;
surf(abs(Xddv1))

figure;
surf(abs(Xddv2))

figure;
surf(abs(diffdd)) % sonuclar ayni





clear;
close all;

Na = 64;
phi = 0.1;
phis = (-45:45) .* (pi/180);

U = getU(Na, 4, 6, 1);  
F = getF(Na, 1);

aaH1 = zeros(Na, Na, length(phis));
C1 = zeros(4, 1, length(phis));
for idx = 1:length(phis)
    phi = phis(idx);
    a1 = array_response(phi, Na);
    aaH1(:, :, idx) = a1 * a1';
    C1(:, :, idx) = U' * a1 * a1' * F;
end

a2 = array_responsev2(phis, Na);
a2 = permute(a2, [1,3,2]);
a2H = conj(permute(a2, [2,1,3]));

aaH2 = pagemtimes(a2, a2H);

diffaaH = aaH1 - aaH2;

C2 = pagemtimes(U', pagemtimes(aaH2, F));

diffC = C1 - C2;
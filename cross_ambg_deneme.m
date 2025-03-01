T = 1;
taus = (-10:0.1:10)';
nus = -10:0.1:10;

C = zeros(length(taus), length(nus));

for idx1 = 1:length(taus)
    for idx2 = 1:length(nus)
        C(idx1,idx2) = cross_ambg(taus(idx1), nus(idx2), T);
    end
end

figure
surf(abs(C))

xlabel('x')
ylabel('y')
C1 = C;

C = cross_ambgv2(taus, nus, T);

figure
surf(abs(C))

xlabel('x')
ylabel('y')

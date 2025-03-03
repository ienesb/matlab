function S = getS(nu, tau, phi, ys, x, B)
    num = 0;
    denum = 0;
    for b = 1:B
        G = getGb(nu, tau, phi, b, B);
        num = num + ys(:, b)'*G*x;
        denum = denum + norm(G*x)^2;
    end
    num = abs(num)^2;
    S = num / denum;
end
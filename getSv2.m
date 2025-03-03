function S = getSv2(nu, tau, phi, ys, x)
    B = 6;

    num = 0;
    denum = 0;
    for b = 1:B
        Gbs = getGbv2(nu, tau, phi, b);
        temp = pagemtimes(Gbs,x);
        num = num + pagemtimes(ys(:, b)', temp);
        denum = denum + pagenorm(temp).^2;
    end
    num = abs(num).^2;
    num = num(:);
    denum = denum(:);
    S = num ./ denum;
end
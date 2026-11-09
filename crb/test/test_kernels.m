function r = test_kernels()
% ICI / ISI kernels: quadrature identity, singular limits, split identity,
% and the reduced-derivative relation to the legacy full derivatives (Phase 0.2).
    r = [];
    df = 15e3;  T = 1/df;  N = 16;
    [NP, NN] = ndgrid(0:N-1, 0:N-1);  Q = NP - NN;

    % (a) quadrature identity I^(k)(q,nu) = (1/T) int_0^T t^k e^{j2pi(q df + nu)t} dt
    %     including singular points q + nu T = 0.
    nus = [444.4, 0, -df, 0.37*df, 1e-7, -2*df + 3e-3];
    qs  = [0 1 -1 2 -3 5];
    errs = zeros(3,1);
    for i = 1:numel(nus)
        nu = nus(i);
        [I, Ip, Ipp] = ici_kernel(qs, nu, T);
        for k = 1:numel(qs)
            q = qs(k);  w = 2*pi*(q*df + nu);
            f0 = integral(@(t) exp(1j*w*t),        0, T, 'AbsTol', 1e-16, 'RelTol', 1e-12)/T;
            f1 = integral(@(t) t.*exp(1j*w*t),     0, T, 'AbsTol', 1e-20, 'RelTol', 1e-12)/T;
            f2 = integral(@(t) t.^2.*exp(1j*w*t),  0, T, 'AbsTol', 1e-24, 'RelTol', 1e-12)/T;
            % absolute errors on the natural scales (peak values 1, T/2, T^2/3);
            % the kernels vanish exactly at integer q with nu = 0, so a
            % relative error is undefined there.
            errs(1) = max(errs(1), abs(I(k)   - f0));
            errs(2) = max(errs(2), abs(Ip(k)  - f1)/T);
            errs(3) = max(errs(3), abs(Ipp(k) - f2)/T^2);
        end
    end
    r = tc(r, 'quadrature identity I  vs integral', errs(1), 1e-9);
    r = tc(r, 'quadrature identity I'' vs integral', errs(2), 1e-9);
    r = tc(r, 'quadrature identity I'''' vs integral', errs(3), 1e-9);

    % singular limits exactly
    [I0, Ip0, Ipp0] = ici_kernel(0, 0, T);
    r = tc(r, 'limit I(0,0) = 1',       abs(I0 - 1), 1e-15);
    r = tc(r, 'limit I''(0,0) = T/2',   abs(Ip0 - T/2)/T, 1e-15);
    r = tc(r, 'limit I''''(0,0) = +T^2/3', abs(Ipp0 - T^2/3)/T^2, 1e-15);
    [Im1, ~, ~] = ici_kernel(-1, df, T);   % q + nu T = 0 away from the origin
    r = tc(r, 'limit I(-1, df) = 1',     abs(Im1 - 1), 1e-15);

    % no NaN/Inf on the full lag matrix at nu = 0 and at tiny nu (series branch)
    [Iz, Ipz, Ippz] = ici_kernel(Q, 0, T);
    [It, Ipt, Ippt] = ici_kernel(Q, 1e-9, T);
    bad = ~all(isfinite([Iz(:); Ipz(:); Ippz(:); It(:); Ipt(:); Ippt(:)]));
    r = tc(r, 'no NaN/Inf over lags at nu = 0 / 1e-9', double(bad), 0);
    r = tc(r, 'nu -> 0 continuity of I''''', relerr(Ippt, Ippz), 1e-6);

    % (b)-(c) ISI kernel: quadrature, tau = T reduces to the full kernel, limits
    tau = 0.37*T;  nu = 444.4;
    [Ii, Iip, Eb] = isi_kernel(qs, nu, T, tau);
    e_i = 0; e_ip = 0; e_b = 0;
    for k = 1:numel(qs)
        w = 2*pi*(qs(k)*df + nu);
        f0 = integral(@(t) exp(1j*w*t),    0, tau, 'AbsTol', 1e-16, 'RelTol', 1e-12)/T;
        f1 = integral(@(t) t.*exp(1j*w*t), 0, tau, 'AbsTol', 1e-20, 'RelTol', 1e-12)/T;
        e_i  = max(e_i,  abs(Ii(k)  - f0)/abs(f0));
        e_ip = max(e_ip, abs(Iip(k) - f1)/abs(f1));
        e_b  = max(e_b,  abs(Eb(k) - exp(1j*w*tau)/T)*T);
    end
    r = tc(r, 'ISI kernel vs integral',           e_i,  1e-9);
    r = tc(r, 'ISI reduced nu-derivative vs integral', e_ip, 1e-9);
    r = tc(r, 'ISI boundary term Ebnd',           e_b,  1e-12);
    [IT, IpT] = isi_kernel(Q, nu, T, T);
    [IF, IpF] = ici_kernel(Q, nu, T);
    r = tc(r, 'I^ISI(tau = T) == I',              relerr(IT, IF), 1e-12);
    r = tc(r, 'I^ISI''(tau = T) == I''',          relerr(IpT, IpF), 1e-12);
    [Il, Ilp] = isi_kernel(0, 0, T, tau);
    r = tc(r, 'limit I^ISI(0,0) = tau/T',         abs(Il - tau/T)/(tau/T), 1e-14);
    r = tc(r, 'limit I^ISI''(0,0) = tau^2/(2T)',  abs(Ilp - tau^2/(2*T))/(tau^2/(2*T)), 1e-14);

    % (e) relation to the legacy full derivatives
    oracle_paths('add', 'old');
    c = onCleanup(@() oracle_paths('remove', 'old'));
    nu = 444.4;
    [I, Ip, Ipp] = ici_kernel(Q, nu, T);
    r = tc(r, 'I == old get_ici_coeffs',                  relerr(I,   get_ici_coeffs(Q, nu, T)), 1e-13);
    r = tc(r, 'I''  == old derivative / (j2pi)',          relerr(Ip,  get_ici_coeffs_derivative(Q, nu, T)/(1j*2*pi)), 1e-12);
    r = tc(r, 'I'''' == old derivative2 / (j2pi)^2',      relerr(Ipp, get_ici_coeffs_derivative2(Q, nu, T)/(1j*2*pi)^2), 1e-11);
    [Ii, Iip] = isi_kernel(Q, nu, T, tau);
    r = tc(r, 'I^ISI == old get_isi_coeffs',              relerr(Ii,  get_isi_coeffs(Q, nu, T, tau)), 1e-13);
    r = tc(r, 'I^ISI'' == old isi derivative / (j2pi)',   relerr(Iip, get_isi_coeffs_derivative(Q, nu, T, tau)/(1j*2*pi)), 1e-12);
end

function r = test_fim_mt40_liu()
% MT40 bridge to Liu's formulation (port of impl/2605.14451v1/verify_my_derivation.m):
% our nu-masked 4P no-ICI FIM at M = 1 equals Liu's J = (2/s2) Re{H^H Diag(|Q^H s|^2) H}
% with H = [beta adot, a, j a] (eqs. 29-32) after the [tau; R; I] -> [R; I; tau]
% permutation and the tau rescale by Ts; OFDM + PSK gives J == Jbar exactly.
    r = [];
    N = 128;  B = 160e6;  g = make_grid(N, 1, B/N, 0, 'cp');
    ch = chan_mt40(g);  cr = chan_realize(ch, g, []);
    Qh = liu_precoders(N, 32, 4, 1/16, 1/8);
    sigma2 = 1;
    idx = eta_index(cr.P);  active = true(4*cr.P, 1);  active(idx.nu) = false;
    n = (0:N-1).';
    for mt = {'16qam', '16psk'}
        mdm = modem(mt{1});
        s = mdm.random(N, crn_stream(7000, 0, 1)).';
        for wn = {'OFDM', 'SC', 'OTFS', 'AFDM'}
            X = Qh.(wn{1}) * s;
            J = fim(X, 'cp-ofdm', g, cr, false, sigma2);
            Ja = J(active, active);
            % Liu's H in OUR per-target ordering [alpha_R, alpha_I, tau(s)]
            H = zeros(N, 3*cr.P);
            for l = 1:cr.P
                a  = exp(-1j*2*pi*cr.tau(l)/g.Ts*n/N);           % a(tau) with tau in samples
                ad = (-1j*2*pi*n/N) .* a / g.Ts;                  % d a / d tau[s]
                H(:, 3*(l-1)+(1:3)) = [a, 1j*a, cr.alpha(l)*ad];
            end
            Jl = (2/sigma2) * real(H' * (abs(X).^2 .* H));
            r = tc(r, sprintf('4P-masked FIM == Liu eq. 32 (%s, %s)', wn{1}, mdm.name), relerr(Ja, Jl), 1e-10);
        end
        % Jensen: |X|^2 = 1
        Jb = fim(ones(N,1), 'cp-ofdm', g, cr, false, sigma2);
        Xo = Qh.OFDM * s;
        Jo = fim(Xo, 'cp-ofdm', g, cr, false, sigma2);
        if strcmp(mt{1}, '16psk')
            r = tc(r, 'OFDM + 16-PSK realization FIM == Jbar', relerr(Jo, Jb), 1e-12);
        else
            fprintf('    (info) OFDM + 16-QAM realization FIM vs Jbar: rel %.3e (nonzero expected)\n', relerr(Jo, Jb));
        end
    end
    % spread constants
    r = tc(r, 'alpha_OTFS(N2=4) = 0.2795', abs(Qh.alpha_spread.OTFS - 0.27951)/0.27951, 1e-4);
    r = tc(r, 'alpha_AFDM(c1=1/16) = 0.2864', abs(Qh.alpha_spread.AFDM - 0.28641)/0.28641, 1e-4);
    for wn = {'SC', 'OTFS', 'AFDM'}
        r = tc(r, sprintf('Q^H %s unitary', wn{1}), relerr(Qh.(wn{1})'*Qh.(wn{1}), eye(N)), 1e-12);
    end
end

function r = test_rng()
% CRN streams and the in-house modem (Phase 0.7).
    r = [];
    a = rand(crn_stream(7, 2, 5), 5, 1);  b = rand(crn_stream(7, 2, 5), 5, 1);  c = rand(crn_stream(7, 2, 6), 5, 1);
    r = tc(r, 'crn_stream reproducible for equal (seed,sweep,trial)', norm(a - b), 0);
    r = tc(r, 'crn_stream differs across trials', double(norm(a - c) == 0), 0);

    for mt = {'4qam', '16qam', '16psk'}
        mdm = modem(mt{1});
        s = mdm.random(1e5, crn_stream(3, 0, 1));
        r = tc(r, sprintf('%s unit power', mdm.name),      abs(mean(abs(s).^2) - 1), 1e-2);
        r = tc(r, sprintf('%s E[s] = 0', mdm.name),        abs(mean(s)), 2e-2);
        r = tc(r, sprintf('%s E[s^2] = 0', mdm.name),      abs(mean(s.^2)), 2e-2);
        bits = randi([0 1], mdm.k, 1000);
        r = tc(r, sprintf('%s map/demap round trip', mdm.name), double(~isequal(mdm.demap(mdm.map(bits)), bits)), 0);
        % Gray property: nearest neighbours differ in exactly one bit
        A = mdm.alphabet;  Q = numel(A);  worst = 0;
        for i = 1:Q
            d = abs(A - A(i));  d(i) = Inf;  nn = find(abs(d - min(d)) < 1e-9);
            for j = nn.'
                nb = sum(bitget(bitxor(i-1, j-1), 1:mdm.k));
                worst = max(worst, nb);
            end
        end
        r = tc(r, sprintf('%s Gray: nearest neighbours differ by 1 bit', mdm.name), worst - 1, 0);
    end
    r = tc(r, '16-QAM kurtosis 1.32',  abs(modem('16qam').kurtosis - 1.32), 1e-3);
    r = tc(r, '16-PSK kurtosis 1.00',  abs(modem('16psk').kurtosis - 1.00), 1e-12);
    r = tc(r, '4-QAM  kurtosis 1.00',  abs(modem('4qam').kurtosis - 1.00), 1e-12);
    [S, bits] = gen_symbols(8, 4, '16qam', crn_stream(3, 0, 2));
    r = tc(r, 'gen_symbols bits consistent', double(~isequal(reshape(modem('16qam').map(bits), 8, 4), S)), 0);
end

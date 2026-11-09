function [dN, dM, cnt] = pilot_spacing_match(N, M, target, l_tau, k_nu)
% PILOT_SPACING_MATCH  OFDM pilot-lattice spacing (dN, dM) whose pilot count
% |1:dN:N| * |1:dM:M| is closest to `target` cells, subject to the channel-
% estimation oversampling limits dN < N/(2 l_tau), dM < M/(2 k_nu) and dN,dM >= 2.
% Used to match the OFDM pilot overhead to the OTFS pilot+guard cell count.
% Exact equality is generally unattainable (full EVA: target 5248, nearest
% 5472); callers assert a tolerance.  Verbatim port of old/eva/utils/match_pilot_spacing.m.
    dN_max = max(2, floor(N/(2*max(l_tau, 1))));
    dM_max = max(2, floor(M/(2*max(k_nu, 1))));
    best = inf;  dN = 2;  dM = 2;  cnt = 0;
    for a = 2:dN_max
        cn = floor((N-1)/a) + 1;
        for b = 2:dM_max
            cm = floor((M-1)/b) + 1;
            c = cn*cm;
            if abs(c - target) < best
                best = abs(c - target);  dN = a;  dM = b;  cnt = c;
            end
        end
    end
end

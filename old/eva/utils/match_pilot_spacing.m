function [dNp, dMp, cnt] = match_pilot_spacing(N, M, target, l_tau, k_nu)
% MATCH_PILOT_SPACING  Pick the OFDM pilot sublattice spacing (dNp,dMp) whose
% pilot count |1:dNp:N|*|1:dMp:M| is closest to `target` cells, subject to the
% channel-estimation oversampling limits
%   dNp < N/(2 l_tau)  (resolve delay),  dMp < M/(2 k_nu)  (resolve Doppler),
% and dNp,dMp >= 2 (avoid degenerate full-row/column pilot blocks).
% Used to overhead-match OFDM pilots to the OTFS guard+pilot cell count.
    % max(.,1) guards the flat/AWGN case (l_tau=0 or k_nu=0): no delay/Doppler
    % resolution constraint -> spacing bounded only by the grid, not by 1/0=Inf.
    dN_max = max(2, floor(N/(2*max(l_tau,1))));
    dM_max = max(2, floor(M/(2*max(k_nu,1))));
    best = inf;  dNp = 2;  dMp = 2;  cnt = 0;
    for dN = 2:dN_max
        cn = floor((N-1)/dN) + 1;
        for dM = 2:dM_max
            cm = floor((M-1)/dM) + 1;
            c  = cn*cm;
            if abs(c - target) < best
                best = abs(c - target);  dNp = dN;  dMp = dM;  cnt = c;
            end
        end
    end
end

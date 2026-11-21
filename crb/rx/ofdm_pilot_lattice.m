function [d, masks, a_p, info] = ofdm_pilot_lattice(Xd, g, ch, n_target, E_pil, tol, lattice)
% OFDM_PILOT_LATTICE  Scattered pilot lattice for OFDM, overhead-matched to the
% OTFS pilot+guard budget.
%   [d, masks, a_p, info] = ofdm_pilot_lattice(Xd, g, ch, n_target, E_pil, tol, lattice)
% n_target: OTFS non-data cell count; E_pil: total pilot energy (data-referenced,
% shared by both waveforms; OTFS puts it in one cell x_p = sqrt(E_pil), OFDM
% spreads it as a_p = sqrt(E_pil / n_pilot)).  Fairness: the cell count is
% matched within `tol` (default 5 %; printed as info.mismatch) and the pilot
% ENERGY exactly.  masks: pilot, data, n_pilots, m_pilots, dN, dM, n_pilot.
% lattice (optional): [dN dM] pins the spacing instead of the search in
% pilot_spacing_match, whose answer depends on ch.k_nu and so on the speed
% (EVA-reduced, 512 x 128: (37, 2) at 120-360 km/h but (4, 19) at 60 km/h, both
% 896 pilots).  A pinned lattice must still meet the count tolerance and the
% oversampling limits dN < N/(2 l_tau), dM < M/(2 k_nu).
    if nargin < 6 || isempty(tol), tol = 0.05; end
    [N, M] = size(Xd);
    if nargin < 7 || isempty(lattice)
        [dN, dM, cnt] = pilot_spacing_match(N, M, n_target, ch.l_tau, ch.k_nu);
    else
        dN = lattice(1);  dM = lattice(2);
        assert(numel(lattice) == 2 && all(lattice == round(lattice)) && all(lattice >= 2), 'ofdm_pilot_lattice:lattice', 'lattice must be [dN dM], integers >= 2');
        assert(dN <= max(2, floor(N/(2*max(ch.l_tau, 1)))) && dM <= max(2, floor(M/(2*max(ch.k_nu, 1)))), 'ofdm_pilot_lattice:limits', ...
            'pinned lattice (%d, %d) is outside the oversampling limits of this channel', dN, dM);
        cnt = (floor((N-1)/dN) + 1)*(floor((M-1)/dM) + 1);
    end
    n_pilots = 1:dN:N;  m_pilots = 1:dM:M;
    mask = false(N, M);  mask(n_pilots, m_pilots) = true;
    assert(nnz(mask) == cnt, 'ofdm_pilot_lattice:count');
    info.mismatch = abs(cnt - n_target)/n_target;
    assert(info.mismatch <= tol, 'ofdm_pilot_lattice:fairness', ...
        'pilot count %d vs target %d (%.1f%% > %.0f%% tolerance)', cnt, n_target, 100*info.mismatch, 100*tol);
    a_p = sqrt(E_pil/cnt);
    d = Xd;  d(mask) = a_p;
    masks.pilot = mask;  masks.data = ~mask;
    masks.n_pilots = n_pilots;  masks.m_pilots = m_pilots;
    masks.dN = dN;  masks.dM = dM;  masks.n_pilot = cnt;  masks.a_p = a_p;
    info.n_pilot = cnt;  info.n_target = n_target;  info.E_pil = cnt*a_p^2;
end

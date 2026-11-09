function H_hat = ofdm_ce(Y, d, masks, l_tau, N0, mdm, mode, n_dd)
% OFDM_CE  OFDM channel estimation from the scattered pilot lattice
% (old/eva/ofdm_ce_demod.tex; port of the ce_ofdm local of old/eva/mcrb_ce_ici.m).
%   'interp' : LS at the pilots -> 2-D linear interpolation.  Chords the phase
%              rotation between frequency pilots and floors the BER (~1e-2).
%   'dft'    : LS fit of the (l_tau+1) delay taps across the pilot comb per
%              pilot symbol, then linear interpolation in time.  Needs
%              l_tau < N/dN (delay Nyquist on the comb).
%   'dftdd'  : 'dft' + n_dd decision-directed passes re-fitting the taps from
%              ALL N subcarriers of every symbol (pilots known, data
%              hard-decided).  Reaches the genie-CE bound.
% mdm = modem(...) supplies the hard decisions.
    if nargin < 8 || isempty(n_dd), n_dd = 2; end
    [N, M] = size(Y);
    n_pilots = masks.n_pilots;  m_pilots = masks.m_pilots;
    H_ls = Y(n_pilots, m_pilots) ./ d(n_pilots, m_pilots);
    [Nq, Mq] = ndgrid(1:N, 1:M);
    Afull = exp(-1j*2*pi*(0:N-1).'*(0:l_tau)/N);
    switch lower(mode)
        case 'interp'
            F_H = griddedInterpolant({n_pilots, m_pilots}, H_ls, 'linear', 'nearest');
            H_hat = F_H(Nq, Mq);
            return
        case {'dft', 'dftdd'}
            assert(l_tau < N/masks.dN, 'ofdm_ce:nyquist', 'dft CE needs l_tau < N/dN (%d vs %g)', l_tau, N/masks.dN);
            Abas = exp(-1j*2*pi*(n_pilots(:)-1)*(0:l_tau)/N);
            Hp   = Afull * (Abas \ H_ls);
            F_H  = griddedInterpolant({(1:N).', m_pilots}, Hp, 'linear', 'nearest');
            H_hat = F_H(Nq, Mq);
        otherwise
            error('ofdm_ce:mode', 'unknown ce mode %s', mode);
    end
    if strcmpi(mode, 'dftdd')
        for it = 1:n_dd
            Xs = conj(H_hat).*Y ./ (abs(H_hat).^2 + N0);
            dhat = reshape(mdm.map(mdm.demap(Xs(:))), N, M);
            dhat(masks.pilot) = d(masks.pilot);
            H_hat = Afull * (Afull \ (Y ./ dhat));
        end
    end
end

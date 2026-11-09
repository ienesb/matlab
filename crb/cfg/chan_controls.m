function ch = chan_controls(g, kind)
% CHAN_CONTROLS  Idealized control channels (channels.md Sec. 4), tests only.
%   ch = chan_controls(g, 'awgn' | 'flat_rayleigh')
%   awgn          : one unit-magnitude tap, random phase, tau = nu = 0
%   flat_rayleigh : one tap alpha ~ CN(0,1), tau = nu = 0
    ch.family = 'control';
    ch.kind   = lower(kind);
    ch.name   = upper(kind);
    ch.P = 1;  ch.tau = 0;  ch.tau_phys = 0;  ch.l_tap = 0;  ch.l_tau = 0;
    ch.nu_max = 0;  ch.k_nu = 0;  ch.P_lin = 1;  ch.E = 0;  ch.snap = true;
    ch.doppler_scale = 0;
    switch ch.kind
        case 'awgn',          ch.gain_model = 'constmag';
        case 'flat_rayleigh', ch.gain_model = 'rayleigh';
        otherwise, error('chan_controls:kind', 'unknown control %s', kind);
    end
    ch.meta.collisions = {};
    ch.pilot.l_p = 0;  ch.pilot.k_p = floor(g.M/2);
end

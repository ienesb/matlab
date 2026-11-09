function ch = chan_mt40(g, varargin)
% CHAN_MT40  MT40 channel specification (channels.md Sec. 3): L static
% delay-only targets, continuous distinct delays in [0, N-1] samples.
%   ch = chan_mt40(g, 'name', value, ...)
% Options (defaults):
%   'L'           40          number of targets (10 for MT40-10)
%   'gain_model'  'unit'      | 'randphase' | 'cn01'
%   'doppler'     'none'      | 'mt40d' (nu_l ~ U[-nu_max, nu_max], nu_max = 0.5/(M TS))
%   'seed'        20260514    geometry seed (MATLAB default generator)
%   'seed_dopp'   20260707    Doppler seed for 'mt40d'
%   'min_sep'     1           minimum delay separation [samples]
%   'integer'     false       round delays to samples (MT40-int)
% Geometry is drawn ONCE here and is fixed for every trial; chan_realize only
% draws gains when gain_model ~= 'unit'.  Bit-identical to
% impl/2605.14451v1/run_sim.m and the old 'paper' scenario for the defaults.
    o = struct('L', 40, 'gain_model', 'unit', 'doppler', 'none', 'seed', 20260514, ...
        'seed_dopp', 20260707, 'min_sep', 1, 'integer', false);
    for i = 1:2:numel(varargin)
        assert(isfield(o, varargin{i}), 'chan_mt40:opt', 'unknown option %s', varargin{i});
        o.(varargin{i}) = varargin{i+1};
    end
    L = o.L;

    rng(o.seed);                                   % default (twister) generator
    usable   = (g.N - 1) - (L - 1)*o.min_sep;
    u        = sort(usable * rand(L, 1));
    tau_samp = u + o.min_sep*(0:L-1).';
    if o.integer, tau_samp = round(tau_samp); end

    ch.family   = 'mt40';
    ch.name     = sprintf('MT40%s', ternary(L ~= 40, sprintf('-%d', L), ''));
    if o.integer, ch.name = [ch.name '-int']; end
    ch.P        = L;
    ch.tau_samp = tau_samp;
    ch.tau      = tau_samp * g.Ts;
    ch.tau_phys = ch.tau;
    ch.l_tap    = round(tau_samp);
    ch.l_tau    = g.N - 1;                          % delays may span the block
    ch.gain_model = lower(o.gain_model);
    ch.P_lin    = ones(L, 1);
    ch.beta     = ones(L, 1);

    switch lower(o.doppler)
        case 'none'
            ch.nu = zeros(L, 1);  ch.nu_max = 0;
        case 'mt40d'
            ch.name = [ch.name '-D'];
            ch.nu_max = 0.5/(g.M*g.TS);
            rng(o.seed_dopp);
            ch.nu = ch.nu_max * (2*rand(L, 1) - 1);
        otherwise
            error('chan_mt40:doppler', 'unknown doppler option %s', o.doppler);
    end
    ch.k_nu = ceil(ch.nu_max * g.M * g.TS);
    ch.range_unit_m = 3e8/(2*g.B);                  % 0.9375 m at 160 MHz (monostatic)
    ch.E = 0;
    ch.opts = o;
    ch.meta.collisions = {};
    ch.snap = o.integer;
end

function out = ternary(c, a, b)
    if c, out = a; else, out = b; end
end

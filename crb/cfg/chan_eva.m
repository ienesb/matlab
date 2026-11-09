function ch = chan_eva(g, varargin)
% CHAN_EVA  EVA channel specification (channels.md Sec. 2).
%   ch = chan_eva(g, 'name', value, ...)
% Options (defaults):
%   'reduced'         false   3-tap delay-sorted EVA-reduced [0 150 370] ns
%   'gain_model'      'rayleigh' | 'constmag'
%   'snap'            true    delays snapped to the tap grid l_p = round(tau N df)
%   'doppler_scale'   1       multiplies nu_max (0 = no Doppler)
%   'fc'              4e9     carrier [Hz]
%   'v_kmh'           120     speed [km/h]
%   'merge_colliding' false   sum gains of taps that share a delay bin
%   'legacy_preset'   false   old/eva/verify_crb.m deterministic 3-tap preset
%                             (power-sorted, fc = 2 GHz, fixed angles/phases,
%                             physical delays, UN-normalized gains)
%   'E'               10      DD Doppler-spread truncation (receiver constant)
%   'normalize'       true    sum P_l = 1
%   'taus_ns'/'pows_db'  []   custom delay [ns] / power [dB] profile with the EVA
%                             draw rules; 'name' labels it (default 'EVA-custom').
%                             Overrides 'reduced'.
%   'variant'         ''      named preset: 'full' | 'reduced' | 'reduced-fast'
%                             ('reduced-fast' = EVA-reduced at 240 km/h)
% The per-trial draw is done by chan_realize(ch, rs).
    o = struct('reduced', false, 'gain_model', 'rayleigh', 'snap', true, ...
        'doppler_scale', 1, 'fc', 4e9, 'v_kmh', 120, 'merge_colliding', false, ...
        'legacy_preset', false, 'E', 10, 'normalize', true, ...
        'taus_ns', [], 'pows_db', [], 'name', 'EVA-custom', 'variant', '');
    o = set_opts(o, varargin{:});
    % named presets (channels.md Sec. 2.5); they only set the options below
    switch lower(o.variant)
        case '',        % nothing: use the explicit options
        case 'full',    o.reduced = false;
        case 'reduced', o.reduced = true;
        case {'reduced-fast', 'reduced_fast'}
            % EVA-reduced at 240 km/h: nu_max = 888.9 Hz, k_nu = 8, nu_max*T = 0.059.
            % The Doppler spread per path doubles while the delay profile is
            % unchanged, which is the regime where CP-OFDM's lower bound beats
            % CP-OTFS's in BOTH the delay and the Doppler block (see README).
            o.reduced = true;  o.v_kmh = 240;
        otherwise, error('chan_eva:variant', 'unknown variant %s', o.variant);
    end

    ch.family = 'eva';
    ch.legacy = o.legacy_preset;
    if o.legacy_preset
        % old/eva/verify_crb.m: power-sorted pairing assigned BY INDEX.
        ch.name      = 'EVA-legacy';
        ch.tau_phys  = [0; 370; 150]*1e-9;
        ch.pow_db    = [0; -0.6; -1.4];
        ch.theta     = [0.5; 2.0; 3.5];
        ch.phi       = [0.7; 1.1; -0.5];
        o.fc = 2e9;  o.snap = false;  o.normalize = false;  o.gain_model = 'constmag';
    elseif ~isempty(o.taus_ns)
        % custom EVA-derived profile (channels.md Sec. 2.5): arbitrary delays and
        % powers with the EVA draw rules (Rayleigh gains, Jakes Doppler, snapping)
        assert(numel(o.taus_ns) == numel(o.pows_db), 'chan_eva:custom', 'taus_ns and pows_db must have the same length');
        ch.name      = o.name;
        ch.tau_phys  = o.taus_ns(:)*1e-9;
        ch.pow_db    = o.pows_db(:);
    elseif o.reduced
        ch.name      = 'EVA-reduced';
        if any(strcmpi(o.variant, {'reduced-fast', 'reduced_fast'})), ch.name = 'EVA-reduced-fast'; end
        ch.tau_phys  = [0; 150; 370]*1e-9;
        ch.pow_db    = [0; -1.4; -0.6];
    else
        ch.name      = 'EVA';
        ch.tau_phys  = [0; 30; 150; 310; 370; 710; 1090; 1730; 2510]*1e-9;
        ch.pow_db    = [0; -1.5; -1.4; -3.6; -0.6; -9.1; -7.0; -12.0; -16.9];
    end
    P_lin = 10.^(ch.pow_db/10);
    if o.normalize, P_lin = P_lin / sum(P_lin); end
    ch.P_lin = P_lin;                      % E|alpha_l|^2
    ch.P     = numel(P_lin);

    c0 = 3e8;
    ch.fc = o.fc;  ch.v = o.v_kmh/3.6;
    ch.nu_max = ch.v * ch.fc / c0;         % 444.44 Hz at 120 km/h, 4 GHz
    ch.doppler_scale = o.doppler_scale;
    ch.gain_model = lower(o.gain_model);
    ch.snap = o.snap;
    ch.E = o.E;

    ch.l_tap = round(ch.tau_phys * g.N * g.deltaf);
    if o.snap
        ch.tau = ch.l_tap * g.Ts;
    else
        ch.tau = ch.tau_phys;
    end
    ch.l_tau = ceil(max(ch.tau_phys) * g.N * g.deltaf);     % guard bound (20 at 512/15 kHz)
    ch.k_nu  = ceil(ch.nu_max * g.M * g.TS);                 % Doppler guard bound
    ch.pilot.l_p = ch.l_tau;
    ch.pilot.k_p = floor(g.M/2);

    % Tap collisions (paths sharing a delay bin after snapping).
    [ut, ~, ic] = unique(ch.l_tap);
    coll = {};
    for i = 1:numel(ut)
        members = find(ic == i);
        if numel(members) > 1, coll{end+1} = members.'; end %#ok<AGROW>
    end
    ch.meta.collisions = coll;
    ch.merge_colliding = o.merge_colliding;
    if o.merge_colliding && ~isempty(coll)
        keep = true(ch.P, 1);
        ch.merge_groups = coll;
        for i = 1:numel(coll)
            keep(coll{i}(2:end)) = false;
        end
        ch.merge_keep = keep;
    end
    ch.opts = o;
end

function o = set_opts(o, varargin)
    for i = 1:2:numel(varargin)
        key = varargin{i};
        assert(isfield(o, key), 'chan_eva:opt', 'unknown option %s', key);
        o.(key) = varargin{i+1};
    end
end

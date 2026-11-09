function g = make_grid(N, M, deltaf, Tcp, cp_mode)
% MAKE_GRID  Time-frequency / delay-Doppler grid constants (latex/io README).
%
%   g = make_grid(N, M, deltaf, Tcp, cp_mode)
%
% Inputs
%   N       subcarriers  = delay bins   (rows of every N x M grid)
%   M       OFDM symbols = Doppler bins (columns)
%   deltaf  subcarrier spacing [Hz]
%   Tcp     cyclic-prefix length [s].  'cp' mode: one CP per OFDM symbol;
%           'rcp' mode: one CP per frame (does not enter the symbol period).
%   cp_mode 'cp' (default) or 'rcp'
%
% Output struct g
%   T       = 1/deltaf          useful symbol duration
%   Ts      = 1/(N*deltaf)      sample period            (delay resolution)
%   TS      = T + Tcp ('cp') | T ('rcp')   symbol period used in every
%             Doppler phase  e^{j2pi m TS nu}  and in the Doppler kernel B_p
%   B       = N*deltaf          bandwidth
%   bin_tau = Ts,  bin_nu = 1/(M*TS)
%   Tframe  = M*TS ('cp') | M*T + Tcp ('rcp')
%
% Three distinct time constants: never alias Ts (sample) with TS (symbol).
    if nargin < 4 || isempty(Tcp),     Tcp = 0;       end
    if nargin < 5 || isempty(cp_mode), cp_mode = 'cp'; end
    cp_mode = lower(cp_mode);
    assert(any(strcmp(cp_mode, {'cp','rcp'})), 'make_grid:mode', 'cp_mode must be ''cp'' or ''rcp''.');
    assert(Tcp >= 0, 'make_grid:Tcp', 'Tcp must be >= 0.');

    g.N = N;  g.M = M;  g.deltaf = deltaf;  g.Tcp = Tcp;  g.cp_mode = cp_mode;
    g.T  = 1/deltaf;
    g.Ts = 1/(N*deltaf);
    g.B  = N*deltaf;
    switch cp_mode
        case 'cp'
            g.TS     = g.T + Tcp;
            g.Tframe = M*g.TS;
        case 'rcp'
            g.TS     = g.T;
            g.Tframe = M*g.T + Tcp;
    end
    g.bin_tau = g.Ts;
    g.bin_nu  = 1/(M*g.TS);
    g.Lcp     = round(Tcp/g.Ts);      % CP length in samples (informational)
end

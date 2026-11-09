function sys = sys_params(M, N, df, maxdelay)
% SYS_PARAMS  System parameter struct for the WPMC-2023 CRLB study.
%   M        : number of subcarriers (delay bins)               [paper: M]
%   N        : number of symbols      (Doppler bins)            [paper: N]
%   df       : subcarrier spacing (SCS) in Hz                   [paper: Delta f]
%   maxdelay : maximum path delay (s), used to size the CP      (CP-OFDM only)
%
% Physical resource: bandwidth B = M*df, symbol duration T = 1/df,
% sample spacing T/M (we take M samples per OFDM/OTFS symbol).  See Sec. II.
T   = 1/df;                 % symbol duration  T = 1/Delta f
dt  = T/M;                  % sample spacing   T/M
Lcp = max(1, ceil(maxdelay/dt));   % CP length in samples (>= max delay)
Tcp = Lcp*dt;               % CP duration
Tp  = Tcp + T;              % CP-OFDM symbol duration  T' = Tcp + T  (Sec. II.B)

sys = struct('M',M,'N',N,'df',df,'T',T,'dt',dt, ...
             'Tcp',Tcp,'Tp',Tp,'NM',N*M);
end

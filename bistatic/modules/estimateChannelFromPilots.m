function [Hhat, noiseVar, Hpilots] = estimateChannelFromPilots(Y, Xp, pilotMask, varargin)
%ESTIMATECHANNELFROMPILOTS  Pilot-aided OFDM channel estimation (LS + 2D interp)
%
% Inputs:
%   Y         : Nsc x Nsym complex received grid (FFT output per OFDM symbol)
%   Xp        : Nsc x Nsym grid with pilot symbols at pilot REs (zeros elsewhere ok)
%   pilotMask : logical Nsc x Nsym, true where pilots are present
% Optional name-values:
%   'Interp'  : 'linear' (default) | 'natural' | 'nearest' | 'cubic'
%   'AveragingWindow' : [F T] neighborhood (in REs) to locally average pilot LS estimates (default [0 0] = none)
%
% Outputs:
%   Hhat      : Nsc x Nsym estimated channel over all REs
%   noiseVar  : scalar noise variance estimate (per RE) from pilot residuals
%   Hpilots   : Nsc x Nsym grid with LS pilot-only estimates (NaN elsewhere)

opts = inputParser;
opts.addParameter('Interp','linear');
opts.addParameter('AveragingWindow',[0 0]);
opts.parse(varargin{:});
method   = opts.Results.Interp;
avgWin   = opts.Results.AveragingWindow;

[Nsc,Nsym] = size(Y);
assert(isequal(size(Xp),[Nsc Nsym]) && isequal(size(pilotMask),[Nsc Nsym]), ...
  'Y, Xp, pilotMask must be Nsc x Nsym');

% --- 1) LS estimate at pilot positions: H_LS(p) = Yp ./ Xp
Yp = Y(pilotMask);
Xp_il = Xp(pilotMask);
assert(all(Xp_il~=0),'Xp must be nonzero at pilot positions.');
Hls_p = Yp ./ Xp_il;

% Place into a sparse grid for clarity
Hpilots          = nan(Nsc,Nsym);
Hpilots(pilotMask)= Hls_p;

% --- 2) (Optional) local pilot averaging like LTE docs do before 2-D interp
if any(avgWin>0)
    [Fgrid,Tgrid] = ndgrid(1:Nsc,1:Nsym);
    Hp = Hpilots;
    for f = 1:Nsc
        f0 = max(1, f-avgWin(1)); f1 = min(Nsc, f+avgWin(1));
        for t = 1:Nsym
            t0 = max(1, t-avgWin(2)); t1 = min(Nsym, t+avgWin(2));
            maskWin = pilotMask(f0:f1,t0:t1);
            if maskWin(:)'*1>0 && pilotMask(f,t)
                vals = Hp(f0:f1,t0:t1);
                Hpilots(f,t) = mean(vals(maskWin),'omitnan');
            end
        end
    end
end

% --- 3) 2-D interpolation over time-frequency grid (as in LTE examples)
% Gather pilot coordinates and values
[fp,tp] = find(pilotMask);
vp      = Hpilots(pilotMask);

% Use scatteredInterpolant (robust to irregular pilot patterns)
F = scatteredInterpolant(double(fp), double(tp), vp, method, 'nearest');
[FF,TT] = ndgrid(1:Nsc, 1:Nsym);
Hhat = F(double(FF), double(TT));

% --- 4) Noise variance estimate from pilot residuals (after smoothing)
% Residual on pilots: r_p = Yp - Xp .* Hhat(pilots)
Hhat_p = Hhat(pilotMask);
rp     = Yp - Xp_il .* Hhat_p;
noiseVar = mean(abs(rp).^2);

end
function out = ct_crlb_floor(g, cr, snr_lin)
% CT_CRLB_FLOOR  Continuous-time CRLB floor and flat-spectrum reference
% (latex/fim/ct_crlb_benchmark.tex eqs. floor, flat, and Sec. "Where OFDM/OTFS land").
%   out = ct_crlb_floor(g, cr, snr_lin)     snr_lin = 1/sigma2 (unit-power data)
% Per target p, with frame SNR  SNR_p = |alpha_p|^2 N M / sigma2 :
%   floor:  CRLB(tau) >= 1/(2 pi^2 W^2 SNR_p),  CRLB(nu) >= 1/(2 pi^2 Tfr^2 SNR_p)
%   flat :  3 x floor  (rectangular envelope / flat spectrum; OFDM and OTFS
%           attain it up to N^2/(N^2-1), M^2/(M^2-1))
% W = N df, Tfr = M TS.  Fields: tau_floor, nu_floor, tau_flat, nu_flat (P x 1)
% and their sums tr_*.  NOTE: the tex prose (line 312) says "6x above"; its own
% eq. (flat) and the EFI entry give the factor 3 -- the 3 is used here.
    W   = g.N*g.deltaf;
    Tfr = g.M*g.TS;
    SNRp = abs(cr.alpha(:)).^2 * g.N * g.M * snr_lin;
    out.tau_floor = 1 ./ (2*pi^2*W^2*SNRp);
    out.nu_floor  = 1 ./ (2*pi^2*Tfr^2*SNRp);
    out.tau_flat  = 3*out.tau_floor;
    out.nu_flat   = 3*out.nu_floor;
    out.tr_tau_floor = sum(out.tau_floor);  out.tr_nu_floor = sum(out.nu_floor);
    out.tr_tau_flat  = sum(out.tau_flat);   out.tr_nu_flat  = sum(out.nu_flat);
end

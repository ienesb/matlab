function [m0, m1, m2] = exp_moments(D)
% EXP_MOMENTS  Normalized exponential moments  m_k(D) = int_0^1 u^k e^{D u} du,
% k = 0, 1, 2, evaluated elementwise for complex D with the removable
% singularity at D = 0 handled by a power series (m_k(0) = 1/(k+1)).
%
% Every kernel of latex/io (cp_ofdm_io.tex, rcp_ofdm_io.tex) is one of these
% moments up to a scale:  I(q) = m0(D), I' = T m1(D), I'' = T^2 m2(D) with
% D = j2pi(q + nu T), and the split kernel I^ISI = (tau/T) m0(D tau/T).
%
% Closed forms (used for |D| >= 1):
%   m0 = (e^D - 1)/D
%   m1 = (D e^D - e^D + 1)/D^2
%   m2 = (D^2 e^D - 2 D e^D + 2 e^D - 2)/D^3
% Series (used for |D| < 1, avoids catastrophic cancellation):
%   m_k = sum_{n>=0} D^n / (n! (n+k+1)).
    m0 = zeros(size(D));  m1 = m0;  m2 = m0;
    small = abs(D) < 1;

    % ---- series branch ---------------------------------------------------
    if any(small(:))
        Ds = D(small);
        s0 = zeros(size(Ds));  s1 = s0;  s2 = s0;
        term = ones(size(Ds));            % D^n / n!
        for n = 0:30
            if n > 0, term = term .* Ds / n; end
            s0 = s0 + term / (n+1);
            s1 = s1 + term / (n+2);
            s2 = s2 + term / (n+3);
        end
        m0(small) = s0;  m1(small) = s1;  m2(small) = s2;
    end

    % ---- closed-form branch ------------------------------------------------
    if any(~small(:))
        Dl = D(~small);
        E  = exp(Dl);
        m0(~small) = (E - 1) ./ Dl;
        m1(~small) = (Dl.*E - E + 1) ./ Dl.^2;
        m2(~small) = (Dl.^2 .* E - 2*Dl.*E + 2*E - 2) ./ Dl.^3;
    end
end

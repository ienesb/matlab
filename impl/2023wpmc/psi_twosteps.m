function Psi = psi_twosteps(h, tau, nu, sys)
% PSI_TWOSTEPS  Channel coefficient matrix Psi^p for Two-Steps OTFS (eq.(25),(26)).
%   DD symbols x[k,l] (col=k*M+l+1) -> received DD samples Y[k',l'] (row=k'*M+l'+1)
%     Psi^p_{k',k}[l',l] = (1/(MN)) h e^{j2pi nu (tau + l dt)}
%                          * Dir(nu N T - k' + k, N) * Dir(l'-l-tau M df, M) * fac
%   fac = 1 for l in I1 (l <= M-1-ceil(M tau/T)),
%         e^{-j2pi(nu T + k/N)} for l in I2 (l >= M-ceil(M tau/T)).
%   The e^{j2pi nu tau} factor (eq.(26)) is kept explicitly: it couples tau and nu.
%
% NOTE on the global constant.  Eq.(26) prints "1/MN", but we use 1/sqrt(MN):
%   * Zak-OTFS (eq.(40)) and Two-Steps OTFS share the SAME transmit signal and
%     differ only by an invertible receiver, so their CRLBs must be ~equal.
%   * 1/(MN) makes Two-Steps 144x (=MN) worse than Zak -- physically impossible
%     and opposite to the paper's own conclusion that Two-Steps is best for
%     Doppler.  1/sqrt(MN) (matching eq.(40)) restores Zak ~ Two-Steps.
% The printed "1/MN" is therefore interpreted as 1/sqrt(MN).
M=sys.M; N=sys.N; NM=sys.NM; df=sys.df; T=sys.T; dt=T/M;
Lt = ceil(M*tau/T);
DM = zeros(M,M); zM=(0:M-1);
for lp=0:M-1, for l=0:M-1
    DM(lp+1,l+1) = sum(exp(1j*2*pi*(lp-l-tau*M*df)*zM/M));
end, end
DN = zeros(N,N); zN=(0:N-1);
for kp=0:N-1, for k=0:N-1
    DN(kp+1,k+1) = sum(exp(1j*2*pi*(nu*N*T - kp + k)*zN/N));
end, end
Psi = zeros(NM,NM);
for kp=0:N-1
    for lp=0:M-1
        row = kp*M + lp + 1;
        for k=0:N-1
            dn = DN(kp+1,k+1);
            for l=0:M-1
                if l >= M-Lt, fac = exp(-1j*2*pi*(nu*T + k/N)); else, fac = 1; end
                pref = (1/sqrt(M*N))*h*exp(1j*2*pi*nu*(tau + l*dt));
                col = k*M + l + 1;
                Psi(row,col) = pref*DM(lp+1,l+1)*dn*fac;
            end
        end
    end
end
end

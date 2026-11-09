function Psi = psi_zak(h, tau, nu, sys)
% PSI_ZAK  Channel coefficient matrix Psi^p for Zak-OTFS (eq.(38),(40)).
%   DD symbols x[k,l] (col=k*M+l+1) -> received DD samples Y[k',l'] (row=k'*M+l'+1)
%     Psi^p_{k',k}[l',l] = (1/sqrt(MN)) h e^{j2pi nu l' dt}
%                          * Dir(l'-l-tau df M, M) * Dir(k-k'+nu T N, N) * fac
%   fac = 1 for l' in I1' (l' >= ceil(M tau/T)), e^{-j2pi k'/N} for l' in I2'.
% Global constant 1/sqrt(MN) exactly as printed in eq.(40).
M=sys.M; N=sys.N; NM=sys.NM; df=sys.df; T=sys.T; dt=T/M;
Lt = ceil(M*tau/T);
% Dirichlet tables
DM = zeros(M,M); zM=(0:M-1);
for lp=0:M-1, for l=0:M-1
    DM(lp+1,l+1) = sum(exp(1j*2*pi*(lp-l-tau*df*M)*zM/M));
end, end
DN = zeros(N,N); zN=(0:N-1);
for k=0:N-1, for kp=0:N-1
    DN(k+1,kp+1) = sum(exp(1j*2*pi*(k-kp+nu*T*N)*zN/N));
end, end
Psi = zeros(NM,NM);
for kp=0:N-1
    for lp=0:M-1
        row = kp*M + lp + 1;
        if lp >= Lt, fac = 1; else, fac = exp(-1j*2*pi*kp/N); end   % I1'/I2'
        pref = (1/sqrt(M*N))*h*exp(1j*2*pi*nu*lp*dt)*fac;
        for k=0:N-1
            dn = DN(k+1,kp+1);
            for l=0:M-1
                col = k*M + l + 1;
                Psi(row,col) = pref*DM(lp+1,l+1)*dn;
            end
        end
    end
end
end

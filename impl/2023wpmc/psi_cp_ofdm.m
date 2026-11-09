function Psi = psi_cp_ofdm(h, tau, nu, sys)
% PSI_CP_OFDM  Channel coefficient matrix Psi^p for CP-OFDM (eq.(10),(11)).
%   x[n,m] (col=n*M+m+1) -> y[n',k] (row=n'*M+k+1).  CP (>= max delay) removes
%   ISI, so only the current symbol n=n' contributes:
%     Psi^p_{n',k}[m] = (1/sqrt(M)) h e^{-j2pi m df tau} e^{j2pi nu (n'T'+Tcp)}
%                       * Dir(m-k+nu T, M)
% Global constant 1/sqrt(M) exactly as printed in eq.(11).
M=sys.M; N=sys.N; NM=sys.NM; df=sys.df; T=sys.T; Tp=sys.Tp; Tcp=sys.Tcp;
zz = (0:M-1);
Psi = zeros(NM,NM);
for np = 0:N-1
    pref0 = (1/sqrt(M))*h*exp(1j*2*pi*nu*(np*Tp+Tcp));
    for k = 0:M-1
        row = np*M + k + 1;
        for m = 0:M-1
            D = sum(exp(1j*2*pi*(m-k+nu*T)*zz/M));      % Dir(m-k+nu T, M)
            col = np*M + m + 1;
            Psi(row,col) = pref0*exp(-1j*2*pi*m*df*tau)*D;
        end
    end
end
end

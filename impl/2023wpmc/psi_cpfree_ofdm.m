function Psi = psi_cpfree_ofdm(h, tau, nu, sys)
% PSI_CPFREE_OFDM  Channel coefficient matrix Psi^p for CP-Free OFDM (eq.(6),(7)).
%   Maps TF information symbols x[n,m] (col = n*M+m+1) to OFDM-demodulated
%   samples y[n',k] (row = n'*M+k+1).  Two contributions per path:
%     Psi^{p,2}  multiplies x[n', m]   (current symbol, eq.(7))
%     Psi^{p,1}  multiplies x[n'-1,m]  (ISI from previous symbol)
%   split at i = ceil(M*tau/T):  i in [ceil,M-1] current, i in [0,ceil-1] ISI.
% Global constant 1/sqrt(M) exactly as printed in eq.(7).
M=sys.M; N=sys.N; NM=sys.NM; df=sys.df; T=sys.T;
Lt = ceil(M*tau/T);                 % ceil(M tau / T)
ii = (0:M-1);
Psi = zeros(NM,NM);
for np = 0:N-1
    pref0 = (1/sqrt(M))*h*exp(1j*2*pi*nu*np*T);
    for k = 0:M-1                   % output subcarrier
        row = np*M + k + 1;
        for m = 0:M-1               % input subcarrier
            a = m - k + nu*T;       % combined modulation/Doppler/demod coeff
            E = exp(1j*2*pi*a*ii/M);
            S2 = sum(E(Lt+1:end));  % i = Lt..M-1  -> current symbol
            S1 = sum(E(1:Lt));      % i = 0..Lt-1  -> previous symbol (ISI)
            pref = pref0*exp(-1j*2*pi*m*df*tau);
            cc = np*M + m + 1;                 % x[n',m]
            Psi(row,cc) = Psi(row,cc) + pref*S2;
            if np >= 1                         % x[n'-1,m] (no symbol -1)
                cp = (np-1)*M + m + 1;
                Psi(row,cp) = Psi(row,cp) + pref*S1;
            end
        end
    end
end
end

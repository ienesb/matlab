function J = getFIM_OFDM_ICI_expected(N, M, deltaf, Ts, sigma2, alpha, nu)
% Single-target EXPECTED OFDM Fisher information matrix, ICI channel.
% Closed form for i.i.d. unit-energy, zero-mean QPSK symbols, i.e. the
% data-averaged FIM of Section 3.1.2 / eqs. (139)-(148) of ofdm_fim.tex.
% Parameter ordering: eta = [alphaR, alphaI, tau, nu].
%
% The cross-bin data products X[n',m] X*[n'',m] vanish on average except on
% the diagonal n'=n'', collapsing the random FIM onto the data-independent
% kernel sums below.  The expected FIM depends on the true Doppler nu (through
% the ICI kernel I and its nu-derivative I') but NOT on the true delay tau.
%
% Uses the normalized (unit-peak) ICI kernel of get_ici_coeffs, consistent
% with the realization function getFIM_OFDM_ICI:
%   I(l)  = (e^{j2pi nu T} - 1) / (j2pi (l + nu T)),   I(0)|_{nu=0} = 1.
%
% Inputs:
%   N, M    : grid size (subcarriers x OFDM symbols)
%   deltaf  : subcarrier spacing [Hz]
%   Ts      : OFDM symbol duration including CP [s]
%   sigma2  : noise variance
%   alpha   : complex path gain (scalar)
%   nu      : true path Doppler [Hz] (operating point of the kernel)
%
% Output:
%   J       : 4 x 4 Fisher information matrix.

    aR = real(alpha);  aI = imag(alpha);  a2 = abs(alpha)^2;
    T  = 1/deltaf;

    % Kernel and its nu-derivative on the full N x N lag grid.
    % Convention: row index = n' (0..N-1), column index = n (0..N-1).
    [NP, NN] = ndgrid(0:N-1, 0:N-1);
    diff_n   = NP - NN;                                   % l = n' - n
    I_mat  = get_ici_coeffs(diff_n, nu, T);              % N x N
    % PDF (eq. 8-9) defines I' := (1/j2pi) dI/dnu; get_ici_coeffs_derivative
    % returns the full dI/dnu = j2pi*I', so scale by 1/(j2pi) to obtain I'.
    Ip_mat = get_ici_coeffs_derivative(diff_n, nu, T) / (1j*2*pi);   % N x N

    npcol = (0:N-1).';                                    % n' values (rows)

    % Data-independent kernel sums (eqs. 133-136).  Sum over n (columns) first,
    % then weight by powers of n' (rows).
    %   S_k = sum_{n,n'} (n')^k |I(n'-n)|^2 ,            k = 0,1,2  (real)
    %   T_k = sum_{n,n'} (n')^k I*(n'-n) I'(n'-n) ,      k = 0,1    (complex)
    %   U   = sum_{n,n'} |I'(n'-n)|^2                               (real)
    rowAbsI = sum(abs(I_mat).^2, 2);                      % per n'
    S0 = sum(rowAbsI);
    S1 = sum(npcol    .* rowAbsI);
    S2 = sum(npcol.^2 .* rowAbsI);

    rowTI = sum(conj(I_mat) .* Ip_mat, 2);                % per n'
    Tk0 = sum(rowTI);
    Tk1 = sum(npcol .* rowTI);

    U = sum(abs(Ip_mat(:)).^2);

    % Symbol-time moments that survive the expectation.
    M1 = M*(M-1)/2;          % sum_m m
    M2 = M*(M-1)*(2*M-1)/6;  % sum_m m^2

    J = zeros(4,4);

    % alpha-alpha block (eqs. 139-140)
    J(1,1) = (2*M/sigma2) * S0;
    J(2,2) = J(1,1);
    J(1,2) = 0;  J(2,1) = 0;

    % tau-tau (eq. 141)
    J(3,3) = (8*pi^2*a2*M/sigma2) * deltaf^2 * S2;

    % nu-nu (eq. 142)
    J(4,4) = (8*pi^2*a2/sigma2) * ( Ts^2 * M2 * S0 ...
                                  + Ts   * 2*M1 * real(Tk0) ...
                                  + M    * U );

    % tau-nu (eq. 143)
    J(3,4) = -(8*pi^2*a2*deltaf/sigma2) * ( Ts * M1 * S1 + M * real(Tk1) );
    J(4,3) =  J(3,4);

    % alpha-tau (eqs. 144, 146)
    J(1,3) =  (4*pi*aI*M*deltaf/sigma2) * S1;
    J(2,3) = -(4*pi*aR*M*deltaf/sigma2) * S1;

    % alpha-nu (eqs. 145, 147) with W := Ts*M1*S0 + M*Tk0  (complex)
    W = Ts*M1*S0 + M*Tk0;
    J(1,4) = -(4*pi/sigma2) * imag(alpha * W);
    J(2,4) =  (4*pi/sigma2) * real(alpha * W);

    % Symmetric completion.
    J(3,1) = J(1,3);  J(4,1) = J(1,4);
    J(3,2) = J(2,3);  J(4,2) = J(2,4);
end

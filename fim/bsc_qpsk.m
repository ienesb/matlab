function rx = bsc_qpsk(tx, rho)
% BSC_QPSK  Pass Gray-coded QPSK symbols through a binary symmetric channel.
%   tx  : array of QPSK symbols in {+/-1 +/-j}/sqrt(2) (any shape)
%   rho : BSC bit error probability (per bit, independent on I and Q)
%
% Gray map:  b1 = (real<0),  b2 = (imag<0)
%   00 -> ( 1+j)/sqrt(2)
%   01 -> ( 1-j)/sqrt(2)
%   10 -> (-1+j)/sqrt(2)
%   11 -> (-1-j)/sqrt(2)

    sz = size(tx);
    tx = tx(:);
    n  = numel(tx);

    % Symbols -> bits
    b1 = real(tx) < 0;
    b2 = imag(tx) < 0;

    % Independent BSC on I and Q bits
    flips = rand(n, 2) < rho;
    b1 = xor(b1, flips(:,1));
    b2 = xor(b2, flips(:,2));

    % Bits -> symbols
    re = 1 - 2*b1;
    im = 1 - 2*b2;
    rx = reshape((re + 1j*im) / sqrt(2), sz);
end

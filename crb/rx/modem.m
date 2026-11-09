function mdm = modem(modtype)
% MODEM  Unit-power Gray-labelled constellations without toolbox dependencies.
%   mdm = modem('4qam' | '16qam' | '16psk')
% Fields: name, Q (order), k (bits/symbol), alphabet (Q x 1, alphabet(i+1) is
% the point carrying Gray label i), kurtosis E|s|^4, and handles
%   s    = mdm.map(bits)      bits: k x n logical/0-1 matrix -> symbols 1 x n
%   bits = mdm.demap(s)       nearest-point hard decision -> k x n bits
%   s    = mdm.random(n, rs)  i.i.d. uniform symbols (rs = RandStream or [])
% 4-QAM: (+-1 +-j)/sqrt(2), Gray by construction.  16-QAM: per-axis Gray
% 00,01,11,10 -> -3,-1,+1,+3, /sqrt(10) (kurtosis 1.32).  16-PSK: Gray phase
% index, e^{j2pi k/16} (kurtosis 1).  Liu Assumptions 1-3 hold for all three.
    modtype = lower(strrep(modtype, '-', ''));
    switch modtype
        case {'4qam','qpsk'}
            k = 2;  lev = [1 -1];                       % bit 0 -> +1, bit 1 -> -1
            labels = 0:3;
            b1 = bitshift(labels, -1);  b0 = bitand(labels, 1);
            alphabet = (lev(b1+1) + 1j*lev(b0+1)) / sqrt(2);
            name = '4-QAM';
        case '16qam'
            k = 4;
            % Per-axis Gray labels: 00 -> -3, 01 -> -1, 11 -> +1, 10 -> +3
            gray2 = zeros(1,4);
            gray2(0+1) = -3;  gray2(1+1) = -1;  gray2(3+1) = 1;  gray2(2+1) = 3;
            labels = 0:15;
            hi = bitshift(labels, -2);  lo = bitand(labels, 3);
            alphabet = (gray2(hi+1) + 1j*gray2(lo+1)) / sqrt(10);
            name = '16-QAM';
        case '16psk'
            k = 4;
            labels = 0:15;
            idx = bitxor(labels, bitshift(labels, -1));  % binary -> Gray code value
            % alphabet(label+1) must be the point whose Gray code is `label`:
            % invert: point index m has Gray label gray(m); place e^{j2pi m/16} at gray(m).
            alphabet = zeros(1,16);
            for m = 0:15
                alphabet(idx(m+1)+1) = exp(1j*2*pi*m/16);
            end
            name = '16-PSK';
        otherwise
            error('modem:type', 'unknown modulation %s', modtype);
    end
    alphabet = alphabet(:);
    alphabet = alphabet / sqrt(mean(abs(alphabet).^2));   % unit average power

    mdm.name = name;  mdm.Q = numel(alphabet);  mdm.k = k;
    mdm.alphabet = alphabet;
    mdm.kurtosis = mean(abs(alphabet).^4);
    w = 2.^(k-1:-1:0).';                                   % MSB first
    mdm.map = @(bits) reshape(alphabet(w.' * double(reshape(bits, k, [])) + 1), 1, []);
    mdm.demap = @(s) demap_local(s, alphabet, k);
    mdm.random = @(n, rs) random_local(n, rs, alphabet);
end

function bits = demap_local(s, alphabet, k)
    s = s(:).';
    [~, i] = min(abs(s - alphabet), [], 1);               % Q x n distances
    lab = i - 1;
    bits = zeros(k, numel(s));
    for b = 1:k
        bits(b, :) = bitand(bitshift(lab, -(k-b)), 1);
    end
end

function s = random_local(n, rs, alphabet)
    Q = numel(alphabet);
    if isempty(rs), i = randi(Q, 1, n); else, i = randi(rs, Q, 1, n); end
    s = reshape(alphabet(i), 1, n);
end

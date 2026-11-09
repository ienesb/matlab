function bits = qam4_demod(xhat, syms, bitmap)
% QAM4_DEMOD  Nearest-constellation hard decision of XHAT (K x 1) against
% SYMS (4x1), de-mapped to bits via BITMAP (4x2).  Returns K x 2 bits.
[~, idx] = min(abs(xhat(:) - syms.'), [], 2);
bits = bitmap(idx, :);
end

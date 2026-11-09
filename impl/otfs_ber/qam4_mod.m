function sym = qam4_mod(bits, syms, bitmap)
% QAM4_MOD  Map rows of BITS (K x 2) to constellation points SYMS (4x1)
% according to BITMAP (4x2 Gray map).  Returns K x 1 symbols.
idx = zeros(size(bits,1),1);
for r = 1:size(bitmap,1)
    idx(all(bits == bitmap(r,:), 2)) = r;
end
sym = syms(idx);
end

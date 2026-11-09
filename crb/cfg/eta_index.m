function idx = eta_index(P)
% ETA_INDEX  Index sets into the stacked 4P parameter vector.
%   idx.aR, idx.aI, idx.tau, idx.nu : 1 x P (one entry per target)
%   idx.gain = [aR aI] interleaved per target (1 x 2P), idx.kin = [tau nu] (1 x 2P)
%   idx.target(p) = 4*(p-1) + (1:4)
    base    = 4*(0:P-1);
    idx.aR  = base + 1;
    idx.aI  = base + 2;
    idx.tau = base + 3;
    idx.nu  = base + 4;
    idx.gain = reshape([idx.aR; idx.aI], 1, []);
    idx.kin  = reshape([idx.tau; idx.nu], 1, []);
    idx.target = @(p) 4*(p-1) + (1:4);
end

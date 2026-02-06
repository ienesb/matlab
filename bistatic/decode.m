function [uhat_all, info] = decode(X_hat, params, noiseVar, maxIter)
%DECODE  Soft-decode an N×M QPSK grid with QC-LDPC using ldpcDecode.
%
% [uhat_all, info] = decode(X_hat, params, noiseVar, maxIter)
%
% Inputs:
%   X_hat    : N×M equalized QPSK symbols (complex). Same packing as generator.
%   params   : struct from generate_coded_ofdm_symbols (must contain H,n,k,numCW)
%   noiseVar : per-complex-symbol noise variance AFTER equalization (Es-plane).
%              For AWGN: noiseVar = E{|noise|^2}. If unknown, estimate from pilots or
%              decision errors. (Critical for good LLR scaling.)
%   maxIter  : (optional) LDPC iterations, default 100.
%
% Outputs:
%   uhat_all : k×numCW decoded information bits (logical)
%   info     : struct with fields (ber_hard_pre, meanLLRmag, itersUsed if available)

if nargin < 4 || isempty(maxIter), maxIter = 100; end

% unpack
H      = params.H;
n      = params.n;
k      = params.k;
numCW  = params.numCW;
N      = params.N;
M      = params.M;
mapTag = params.mapping;

% sanity
assert(all(size(X_hat) == [N M]), 'X_hat size mismatch.');
assert(noiseVar > 0, 'noiseVar must be positive.');

% ---------- QPSK LLRs (Gray; bit 0 => +1)
switch mapTag
    case 'QPSK_Gray_bit0_to_plus1'
        y = X_hat(:);
        % For AWGN, LLR for bit=0 vs 1 on I and Q:
        %   LLR_I = (2/σ^2) * Re(y),  LLR_Q = (2/σ^2) * Im(y)
        LLR_I =  (2/ noiseVar) * real(y);
        LLR_Q =  (2/ noiseVar) * imag(y);
        llr_frame = reshape([LLR_I.'; LLR_Q.'], [], 1);  % interleave I,Q per symbol
    otherwise
        error('Unknown mapping tag: %s', mapTag);
end

% ---------- split into codewords and decode
assert(mod(numel(llr_frame), n) == 0 && numel(llr_frame)/n == numCW, ...
    'LLR stream/frame alignment mismatch. Check N,M and LDPC shape.');

llr_cw = reshape(llr_frame, n, numCW);
uhat_all = false(k, numCW);

for cw = 1:numCW
    uhat_all(:,cw) = ldpcDecode(llr_cw(:,cw), H, maxIter);
end

% ---------- simple diagnostics (optional)
% Pre-decode hard BER proxy (on coded bits) if you want:
hardBits = llr_frame < 0; % bit=1 if LLR<0 (since positive supports bit=0)
info = struct();
info.ber_hard_pre = [];  % (not computed vs ground-truth here)
info.meanLLRmag   = mean(abs(llr_frame));
info.itersUsed    = [];  % (ldpcDecode does not return iterations; left blank)
end

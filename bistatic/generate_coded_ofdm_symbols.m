function [X, params, u_all] = generate_coded_ofdm_symbols(N, M, varargin)
    %GENERATE_CODED_OFDM_SYMBOLS  Build rate~0.1 QC-LDPC and QPSK symbols on an N×M grid.
    %
    % [X, params, u_all] = generate_coded_ofdm_symbols(N, M, 'Rate', 0.1, 'Z', 240, 'Seed', 42)
    %
    % Outputs:
    %   X       : N×M complex QPSK symbols (no pilots; insert yours if desired)
    %   params  : struct with LDPC objects & shapes (pass to decode.m)
    %   u_all   : k×numCW matrix of info bits per codeword (logical)
    %
    % Notes:
    %   - Default design is 9×10 prototype -> nominal rate ≈ 1 - 9/10 = 0.1.
    %   - QPSK Gray mapping with bit 0 -> +1, bit 1 -> -1 on each axis.
    %   - If your N,M change, ensure (2*N*M) is a multiple of n = nb*z.
    
    % ---- options
    p = inputParser;
    p.addParameter('Rate', 0.1, @(x)isnumeric(x)&&x>0&&x<1);
    p.addParameter('Z', 240, @(x)isnumeric(x)&&x>=8&&mod(x,1)==0);
    p.addParameter('Seed', 42, @(x)isnumeric(x)&&isscalar(x));
    p.parse(varargin{:});
    RATE = p.Results.Rate;   %#ok<NASGU> % (kept for future variants)
    z    = p.Results.Z;
    seed = p.Results.Seed;
    
    % ---- grid & modulation
    Mmod = 4;                  % QPSK
    bitsPerSym = log2(Mmod);   % = 2
    numSyms = N*M;
    numBits = numSyms*bitsPerSym;
    
    % ---- choose prototype for ~0.1 rate
    mb = 9; nb = 10;                   % 9x10 -> ~0.1
    n  = nb*z;                         % codeword length
    m  = mb*z;                         % parity-check rows
    k  = n - m;                        % info bits per CW
    
    % sanity: make frame length an integer number of codewords
    assert(mod(numBits, n) == 0, ...
        'Frame length (2*N*M=%d) is not a multiple of codeword length n=%d. Adjust N,M or Z.', numBits, n);
    numCW = numBits / n;
    
    % ---- build a prototype with moderate row/col weights
    rng(seed);
    B = -1*ones(mb, nb);
    % simple mask (1 => place a circulant; 0 => zero block)
    mask = [ ...
        1 1 1 0 1 0 1 0 1 0;
        0 1 0 1 1 1 0 1 0 1;
        1 0 1 0 1 0 1 1 0 1;
        1 1 0 1 0 1 0 1 1 0;
        0 1 1 0 1 1 0 0 1 1;
        1 0 1 1 0 1 1 0 0 1;
        0 1 0 1 1 0 1 1 0 1;
        1 0 1 0 1 1 0 1 1 0;
        1 1 0 1 0 0 1 0 1 1];
    for i = 1:mb
        for j = 1:nb
            if mask(i,j) == 1
                B(i,j) = randi([0 z-1]); % random circulant shift
            end
        end
    end
    
    % ---- parity-check matrix (sparse)
    H = ldpcQuasiCyclicMatrix(z, B);
    cfgLDPCEnc = ldpcEncoderConfig(H);
    
    % ---- encode numCW codewords
    u_all = false(k,  numCW);
    c_all = false(n,  numCW);
    for cw = 1:numCW
        u = logical(randi([0 1], k, 1));     % info bits
        c = ldpcEncode(u, cfgLDPCEnc);                % coded bits
        u_all(:,cw) = u;
        c_all(:,cw) = c;
    end
    
    % ---- serialize bits -> QPSK symbols (Gray; bit 0 => +1, bit 1 => -1)
    codedBits = reshape(c_all, [], 1);         % length = numBits
    b = reshape(codedBits, 2, []).';
    sym = ((1 - 2*double(b(:,1)))/sqrt(2)) + 1j*((1 - 2*double(b(:,2)))/sqrt(2));
    
    % ---- map to N×M grid
    X = reshape(sym, N, M);
    
    % ---- pack params for the decoder
    params = struct();
    params.N = N; params.M = M;
    params.mb = mb; params.nb = nb; params.z = z;
    params.n  = n;  params.m = m;  params.k = k;
    params.numCW = numCW;
    params.H = H;
    params.mapping = 'QPSK_Gray_bit0_to_plus1'; % for consistent LLRs
    
    % display
    fprintf('LDPC H: %d x %d (rate≈%.3f), z=%d, n=%d, k=%d, numCW=%d\n', ...
        m, n, 1 - m/n, z, n, k, numCW);
end

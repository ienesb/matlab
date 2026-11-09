function Pa = power_sums(K, a)
% POWER_SUMS  P_a(K) = sum_{k=0}^{K-1} k^a  (cp_ofdm_fim.tex eq. 23).
%   P_0 = K, P_1 = K(K-1)/2, P_2 = K(K-1)(2K-1)/6; higher a by direct sum.
    switch a
        case 0, Pa = K;
        case 1, Pa = K*(K-1)/2;
        case 2, Pa = K*(K-1)*(2*K-1)/6;
        otherwise, Pa = sum((0:K-1).^a);
    end
end

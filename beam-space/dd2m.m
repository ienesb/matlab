function ddm = dd2m(dd, M, N)
    ddm = permute(dd, [3, 1, 4, 2]);
    ddm = reshape(ddm, [M*N, M*N]);
end

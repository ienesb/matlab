function [peak_rows, peak_cols] = find_cluster_peaks(Y, threshold)
% FIND_CLUSTER_PEAKS  Cluster above-threshold entries of Y and return peak indices.
%
%   [peak_rows, peak_cols] = find_cluster_peaks(Y)
%   [peak_rows, peak_cols] = find_cluster_peaks(Y, threshold)
%
%   Y          : N x M matrix (real or complex)
%   threshold  : minimum |Y| to consider non-zero  (default: 0)
%
%   peak_rows  : 1 x K row    indices of the magnitude peak of each cluster
%   peak_cols  : 1 x K column indices of the magnitude peak of each cluster
%   K          : number of clusters found

    if nargin < 2
        threshold = 0;
    end

    [N, M] = size(Y);
    mag     = abs(Y);
    mask    = mag > threshold;
    visited = false(N, M);

    peak_rows = [];
    peak_cols = [];

    % 8-connectivity offsets
    dr = [-1, -1, -1,  0,  0,  1,  1,  1];
    dc = [-1,  0,  1, -1,  1, -1,  0,  1];

    [nz_r, nz_c] = find(mask);

    for i = 1:numel(nz_r)
        r0 = nz_r(i);
        c0 = nz_c(i);
        if visited(r0, c0), continue; end

        % BFS: collect all connected pixels in this cluster
        Q           = [r0, c0];
        visited(r0, c0) = true;
        head        = 1;

        while head <= size(Q, 1)
            r    = Q(head, 1);
            c    = Q(head, 2);
            head = head + 1;

            for d = 1:8
                nr = r + dr(d);
                nc = c + dc(d);
                if nr >= 1 && nr <= N && nc >= 1 && nc <= M ...
                        && mask(nr, nc) && ~visited(nr, nc)
                    visited(nr, nc) = true;
                    Q(end+1, :) = [nr, nc];  %#ok<AGROW>
                end
            end
        end

        % Peak = index of maximum magnitude within the cluster
        cluster_lin = sub2ind([N, M], Q(:, 1), Q(:, 2));
        [~, pk]     = max(mag(cluster_lin));
        peak_rows(end+1) = Q(pk, 1);  %#ok<AGROW>
        peak_cols(end+1) = Q(pk, 2);  %#ok<AGROW>
    end

    peak_rows = peak_rows(:).';
    peak_cols = peak_cols(:).';
end

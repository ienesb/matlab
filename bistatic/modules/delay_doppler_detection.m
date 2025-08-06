function [DD_map, peaks] = delay_doppler_detection(H, params)
    N = params.N;
    M = params.M;

    Fn = dftmtx(N)./sqrt(N);
    Fm = dftmtx(M)./sqrt(M);
    
    % Delay-Doppler transformation: FFT over frequency, IFFT over time
    DD_map_complex = Fn' * H * Fm;
    DD_map = abs(DD_map_complex).^2;
    
    % Simple thresholding for detection
    threshold = mean(DD_map(:)) + 3 * std(DD_map(:));
    [rows, cols] = find(DD_map > threshold);
    peaks = [rows, cols];

    % Estimate delays and dopplers from bin indices
    delay_bins = (0:N-1) * params.delay_res;
    doppler_bins = (-M/2:M/2-1) * params.doppler_res;
end

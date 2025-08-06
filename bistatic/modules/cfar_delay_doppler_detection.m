function [DD_map, peaks, est_delays, est_dopplers] = cfar_delay_doppler_detection(H, params)
    N = params.N;
    M = params.M;
    
    % Step 1: Transform H into delay-Doppler domain
    Fn = dftmtx(N)./sqrt(N);
    Fm = dftmtx(M)./sqrt(M);
    DD_map_complex = Fn' * H * Fm;
    DD_map = abs(DD_map_complex).^2;
    
    % Step 2: Define CFAR parameters
    guard_cells = [2, 2];  % [delay, Doppler]
    training_cells = [8, 8];
    Pfa = 1e-4;
    
    % Step 3: Apply 2D CA-CFAR
    peaks = [];
    for i = 1+training_cells(1)+guard_cells(1) : N-training_cells(1)-guard_cells(1)
        for j = 1+training_cells(2)+guard_cells(2) : M-training_cells(2)-guard_cells(2)
            % Extract training cells (excluding guard and CUT)
            cut = DD_map(i, j);
            training_area = DD_map(i-training_cells(1)-guard_cells(1):i+training_cells(1)+guard_cells(1), ...
                                   j-training_cells(2)-guard_cells(2):j+training_cells(2)+guard_cells(2));
            training_area( training_cells(1)+1:end-training_cells(1), ...
                           training_cells(2)+1:end-training_cells(2)) = 0;  % zero out guard+cut
    
            noise_est = sum(training_area(:)) / ( (2*training_cells(1)+2*guard_cells(1)+1) * ...
                                                  (2*training_cells(2)+2*guard_cells(2)+1) - ...
                                                  (2*guard_cells(1)+1)*(2*guard_cells(2)+1) );
    
            threshold = noise_est * (-log(Pfa));  % for exponential distribution
    
            if cut > threshold
                peaks = [peaks; i, j];
            end
        end
    end
    
    % Step 4: Estimate delay and Doppler
    % delay_bins = (0:N-1) * params.delay_res;
    % doppler_bins = (-M/2:M/2-1) * params.doppler_res;
    % 
    % est_delays = delay_bins(peaks(:,1));
    % est_dopplers = doppler_bins(peaks(:,2) - M/2);

end

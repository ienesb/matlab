function [tau_hats, nu_hats] = detect_targets_omp(H_hat, phis, params)
    v2struct(params);

    I = repmat(eye(M), 1, 1, 9417);
    for idx = 1:5
        temp = I .* pagemtimes(H_hat', phis);
        temp = sum(temp, [1, 2]);
        temp = abs(temp(:));
        [m, max_idx] = max(temp)
        max_idx1 = floor(max_idx/43) + 62;
        max_idx2 = max_idx - (max_idx1-62)*43-1 + 108;
        tau_hat = idx2tau(max_idx1, Tsym, N_fft);
        nu_hat = idx2nu(max_idx2, delta_f, M_fft);
        b = getb(tau_hat, params);
        c = getc(nu_hat, params);
        alpha_hat = abs(b' * H_hat * conj(c)/(N * M * pilot_ratio))
        H_hat = H_hat - alpha_hat * b * c.';
    end



    
    cfar2D = phased.CFARDetector2D('GuardBandSize',guardsize,'TrainingBandSize',trainingsize,...
      'ProbabilityFalseAlarm',pfa, ThresholdOutputPort=true);
    
    [columnInds,rowInds] = meshgrid(108:150, 62:280);
    CUTIdx = [rowInds(:) columnInds(:)]';

    [detections, threshold] = cfar2D(HDD, CUTIdx);
    detection_map = zeros(size(HDD));
    
    r = CUTIdx(1, detections);
    c = CUTIdx(2, detections);
    
    lin = sub2ind(size(HDD), r, c);
    
    detection_map(lin) = HDD(lin);

    epsilon = 5;
    minpts = 10;

    detection_map(lin) = dbscan([r.' c.'], epsilon, minpts);
    detection_map(detection_map < 0) = 0;
    Ntargets = max(detection_map(:));
    
    tau_hats_idx = zeros(Ntargets, 1);
    nu_hats_idx = zeros(Ntargets, 1);
    for k = 1:Ntargets
        temp = zeros(size(HDD));
        temp(detection_map == k) = HDD(detection_map == k);
        [~, target_idx] = max(temp(:));
        [tau_hat_idx, nu_hat_idx] = ind2sub(size(HDD), target_idx);
        tau_hats_idx(k) = tau_hat_idx;
        nu_hats_idx(k) = nu_hat_idx;
    end  
    
    [~, linIdx] = max(detection_map(:));
    [row, col] = ind2sub(size(detection_map), linIdx);
    taus = getDelayArray(Tsym, N_fft);
    nus = getDopplerArray(delta_f, M_fft, is_fftshifted);

    tau_hats = taus(tau_hats_idx);
    nu_hats = nus(nu_hats_idx);
end
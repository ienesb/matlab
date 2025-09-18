function [tau_hats, nu_hats] = detect_targets(HDD, params)
    v2struct(params);
    guardsize = [20, 12];
    trainingsize = [40, 42];
    pfa = 1e-4;
    
    cfar2D = phased.CFARDetector2D('GuardBandSize',guardsize,'TrainingBandSize',trainingsize,...
      'ProbabilityFalseAlarm',pfa, ThresholdOutputPort=true);
    
    [columnInds,rowInds] = meshgrid(108:150, 62:280);
    CUTIdx = [rowInds(:) columnInds(:)]';

    [detections, threshold] = cfar2D(HDD, CUTIdx); %#ok<ASGLU>
    detection_map = zeros(size(HDD));
    
    row = CUTIdx(1, detections);
    col = CUTIdx(2, detections);

    rangesCUT = range_array(row);
    velocitiesCUT = velocity_array(col);
    
    lin = sub2ind(size(HDD), row, col);
    
    detection_map(lin) = HDD(lin);

    epsilon = 10.4897;
    minpts = 1;

    detection_map(lin) = dbscan([rangesCUT.' velocitiesCUT.'], epsilon, minpts);
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
    
    % [~, linIdx] = max(detection_map(:));
    % [row, col] = ind2sub(size(detection_map), linIdx);

    tau_hats = delay_array(tau_hats_idx);
    nu_hats = doppler_array(nu_hats_idx);
end
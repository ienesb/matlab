function [tau, nu] = detect_targets(HDD, params) % assuming there is only 1 target detected. For multiple targets, add clustering !!!!!!!!
    v2struct(params);
    guardsize = [20, 12];
    trainingsize = [40, 42];
    pfa = 1e-5;
    
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
    
    [~, linIdx] = max(detection_map(:));
    [row, col] = ind2sub(size(detection_map), linIdx);
    taus = getDelayArray(Tsym, N_fft);
    nus = getDopplerArray(delta_f, M_fft, is_fftshifted);

    tau = taus(row);
    nu = nus(col);
end
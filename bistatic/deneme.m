main_bistatic_isac;

guardsize = [20, 12];
trainingsize = [40, 42];
paddingsize = [guardsize(1) + trainingsize(1), guardsize(2) + trainingsize(2)];
pfa = 1e-5;

cfar2D = phased.CFARDetector2D('GuardBandSize',guardsize,'TrainingBandSize',trainingsize,...
  'ProbabilityFalseAlarm',pfa, ThresholdOutputPort=true);

% [columnInds,rowInds] = meshgrid(paddingsize(2)+1:M_fft+paddingsize(2), paddingsize(1)+1:N_fft+paddingsize(1));
[columnInds,rowInds] = meshgrid(108:150, 62:280);
CUTIdx = [rowInds(:) columnInds(:)]';

% HDD_padded = padding(HDD, paddingsize(1), paddingsize(2));
HDD_padded = HDD;

[detections, threshold] = cfar2D(HDD_padded,CUTIdx);

detection_map = zeros(size(HDD_padded));

r = CUTIdx(1, detections);
c = CUTIdx(2, detections);

lin = sub2ind(size(HDD_padded), r, c);

detection_map(lin) = HDD(lin);

plotDDMap(pow2db(HDD_padded), params, 1);

plotDDMap(detection_map, params, 1);

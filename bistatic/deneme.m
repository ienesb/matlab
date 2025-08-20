main_bistatic_isac;

N_fft = params.N_fft;
M_fft = params.M_fft;

% DD_map(9,1) = 0;
% DD_map(12,1) = 0;
% DD_map(19,1) = 0;

guardsize = 2;
trainingsize = 10;
paddingsize = guardsize + trainingsize;
pfa = 1e-4;

cfar2D = phased.CFARDetector2D('GuardBandSize',guardsize,'TrainingBandSize',trainingsize,...
  'ProbabilityFalseAlarm',pfa, ThresholdOutputPort=true);

[columnInds,rowInds] = meshgrid(paddingsize+1:M_fft+paddingsize, paddingsize+1:N_fft+paddingsize);
CUTIdx = [rowInds(:) columnInds(:)]';

DD_map_padded = padding(HDD, paddingsize);

[detections, threshold] = cfar2D(DD_map_padded,CUTIdx);

detections = reshape(detections, size(HDD));
threshold = reshape(threshold, size(HDD));

% figure;
% imagesc(mag2db(HDD));

plotDDMap(mag2db(HDD), N_fft, M_fft, params.delta_f);

figure;
imagesc(detections);

figure; hold on;
plot(HDD(:, 1));
plot(threshold(:, 1));
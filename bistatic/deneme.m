main_bistatic_isac;

N_fft = params.N_fft;
M_fft = params.M_fft;

% DD_map(9,1) = 0;
% DD_map(12,1) = 0;
% DD_map(19,1) = 0;

guardsize = [20, 12];
trainingsize = [40, 42];
paddingsize = [guardsize(1) + trainingsize(1), guardsize(2) + trainingsize(2)];
pfa = 1e-4;

cfar2D = phased.CFARDetector2D('GuardBandSize',guardsize,'TrainingBandSize',trainingsize,...
  'ProbabilityFalseAlarm',pfa, ThresholdOutputPort=true);

[columnInds,rowInds] = meshgrid(paddingsize(2)+1:M_fft+paddingsize(2), paddingsize(1)+1:N_fft+paddingsize(1));
CUTIdx = [rowInds(:) columnInds(:)]';

DD_map_padded = padding(HDD, paddingsize(1), paddingsize(2));

[detections, threshold] = cfar2D(DD_map_padded,CUTIdx);

detections = reshape(detections, size(HDD));
threshold = reshape(threshold, size(HDD));

% figure;
% imagesc(mag2db(HDD));

plotDDMap(mag2db(HDD), N_fft, M_fft, params.delta_f);

figure;
imagesc(detections);

delays = linspace(0, params.Tsym*(N_fft-1)/N_fft, N_fft) .* 1e6;
figure; hold on;
plot(delays, HDD(:, 1));
plot(delays, threshold(:, 1));
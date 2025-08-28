clear;
close all;
clc;

load("H.mat")
load("H_hat1.mat")
load("H_hat2.mat")
H_hat = H_hat1;

params = init_simulation_params();
v2struct(params);

%% CFAR Params
guardsize = [20, 12];
trainingsize = [40, 42];
pfa = 1e-4;

cfar2D = phased.CFARDetector2D('GuardBandSize',guardsize,'TrainingBandSize',trainingsize,...
  'ProbabilityFalseAlarm',pfa, ThresholdOutputPort=true);

[columnInds,rowInds] = meshgrid(108:150, 62:280);
CUTIdx = [rowInds(:) columnInds(:)]';

%% Detection

HDD = generate_dd_map(H_hat, params);

[detections, threshold] = cfar2D(HDD, CUTIdx);

detection_map = zeros(size(HDD));

r = CUTIdx(1, detections);
c = CUTIdx(2, detections);

lin = sub2ind(size(HDD), r, c);

% detection_map(lin) = HDD(lin);
% detection_map(lin) = 1;

epsilon = 5;
minpts = 10;
detection_map(lin) = dbscan([r.' c.'], epsilon, minpts);

detection_map_reduced = detection_map(62:280, 108:150);
max(detection_map_reduced(:))

plotDDMap(HDD, params, 1, 1)
plotDDMap(detection_map, params, 1)
hold on;
line([108, 150], [62, 62], 'Color', 'r');
line([108, 150], [280, 280], 'Color', 'r');
line([108, 108], [62, 280], 'Color', 'r');
line([150, 150], [62, 280], 'Color', 'r');

figure;
imagesc(detection_map_reduced);

% close all;
Nc = 100;
Pfa = 1e-4;
release(cfar);
cfar = phased.CFARDetector;
cfar.ThresholdFactor = "Custom";
% cfar.CustomThresholdFactor = Nc*(Pfa^(-1/Nc) - 1);
cfar.CustomThresholdFactor = 1.1;
cfar.ThresholdOutputPort = true;
cfar.NumTrainingCells = Nc;
% cfar.ProbabilityFalseAlarm = Pfa;

[x_detected,th] = cfar(range_profile.',1:length(range_profile));

figure;
stem(x_detected);

figure; hold on;
plot(range_profile);
plot(th);

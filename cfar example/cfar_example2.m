clear;
close all;

npower = db2pow(-10);  % Assume 10dB SNR ratio

cfar = phased.CFARDetector();

cfar.ThresholdOutputPort = true;
cfar.NumTrainingCells = 200;
cfar.ProbabilityFalseAlarm = 1e-3;

cfar

rs = RandStream('mt19937ar','Seed',2010);
Npoints = 1e4;
rsamp = randn(rs,Npoints,1)+1i*randn(rs,Npoints,1);
ramp = linspace(1,10,Npoints)';
xRamp = abs(sqrt(npower*ramp./2).*rsamp).^2;

[x_detected,th] = cfar(xRamp,1:length(xRamp));

figure;
plot(1:length(xRamp), xRamp, 1:length(xRamp), th, find(x_detected), xRamp(x_detected), 'o')
legend('Signal','Threshold','Detections','Location','Northwest')
xlabel('Time Index')
ylabel('Level')

[my_detected, my_threshold] = my_cfar(xRamp);

figure;
plot(1:length(xRamp), xRamp, 1:length(xRamp), my_threshold);
legend('Signal','Threshold','Location','Northwest')
xlabel('Time Index')
ylabel('Level')
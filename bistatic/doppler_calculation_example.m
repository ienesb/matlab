clear;

c = 3e8;
fc = 28e9;
lambda = c / fc;
delta_t = 1/fc; % time difference between two peak points at transmitter 

p1 = [56.9; 10]; % target position and velocity
v = [1.4;-2.2];
% p1 = [79.4; 7]; % target position and velocity
% v = [2.2; -13.7];

tx = [0; 0]; % static transmitter and receiver positions
rx = [50; 0];

d1 = norm(p1 - tx) + norm(p1 - rx);
tau1 = d1/c;

p2 = p1 + delta_t*v;
d2 = norm(p2 - tx) + norm(p2 - rx);
tau2 = d2 / c;

delta_tp = delta_t + tau2 - tau1; % time difference between two peak points at receiver

fcp = 1 / delta_tp;
rho = fcp - fc % doppler shift

nu = velocity2nu(v, p1, tx, rx, lambda)

clear;
close all;
load('Ss_new.mat');

Ndoppler = 30;
Ndelay = 64;
Nangle = 24;

Pfa = 0.15;
Nc = 64;
a = Nc*(Pfa^(-1/Nc)-1);
Ts = ones(Nangle, Ndoppler, Ndelay).*10000;

for idx_angle = 4:Nangle-3
    for idx_doppler = 4:Ndoppler-3 
        for idx_delay = 4:Ndelay-3
            c1 = Ss_new(idx_angle-3:idx_angle+3, idx_doppler-3:idx_doppler+3, idx_delay-3:idx_delay+3);
            c2 = Ss_new(idx_angle-1:idx_angle+1, idx_doppler-1:idx_doppler+1, idx_delay-1:idx_delay+1);
            T = sum(sum(sum(c1))) - sum(sum(sum(c2)));
            T = a*T/Nc;
            Ts(idx_angle,idx_doppler, idx_delay) = T;
            % Ss_new(idx_angle, idx_doppler, idx_delay)
            % T
        end
    end
end

detections = Ss_new > T;

figure;
stem(detections(:));

linearIndices = find(detections == 1); % Get linear indices of ones
[row, col, page] = ind2sub(size(detections), linearIndices); % Convert to 3D indices

disp('Indices of ones:');
disp([row, col, page]); % Display the indices
% Ss_new = zeros(24,30,64);
% idx_delay = 1;
% idx_delayv2 = 1;
% idx_doppler = 1;
% for nu = 1:30
%     for tau = 1:64
%         Ss_new(:,idx_doppler,idx_delayv2) = Ss(:,idx_doppler,idx_delay);
%         idx_delay = idx_delay + 1;
%         idx_delayv2 = idx_delayv2 + 1;
%     end
%     idx_delayv2 = 1;
%     idx_doppler = idx_doppler + 1;
% end

% Ss_new = permute(Ss_new, [2,3,1]);
% for idx = 1:24
%     figure;
%     surf(abs(Ss_new(:,:,idx)));
% end

[maxValue, linearIndex] = max(Ss_new(:)); % Get max value and linear index
[row, col, page] = ind2sub(size(Ss_new), linearIndex); % Convert to 3D indices

disp(['Max value: ', num2str(maxValue)]);
disp(['Index: (', num2str(row), ', ', num2str(col), ', ', num2str(page), ')']);
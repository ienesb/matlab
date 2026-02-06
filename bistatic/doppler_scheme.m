% Coordinates
tx = [0, 0];              % Transmitter
rx = [50, 0];             % Receiver
target = [56.9, 10];      % Target
v_target = [1.4, -2.2];   % Target velocity

% Differences
d_tx = target - tx;  % from TX to target
d_rx = target - rx;  % from RX to target

% Plot
figure; hold on; grid on; axis equal;

% Plot positions
plot(tx(1), tx(2), 'ro', 'MarkerSize', 10, 'DisplayName', 'Transmitter');
plot(rx(1), rx(2), 'bo', 'MarkerSize', 10, 'DisplayName', 'Receiver');
plot(target(1), target(2), 'ks', 'MarkerSize', 10, 'DisplayName', 'Target');

% Draw vectors
quiver(tx(1), tx(2), d_tx(1), d_tx(2), 0, 'r--', 'LineWidth', 1.5, 'DisplayName', 'd_{tx}');
quiver(rx(1), rx(2), d_rx(1), d_rx(2), 0, 'b--', 'LineWidth', 1.5, 'DisplayName', 'd_{rx}');
quiver(target(1), target(2), v_target(1), v_target(2), 0, 'g', 'LineWidth', 1.5, 'DisplayName', 'Velocity');

% Labels
xlabel('x [m]');
ylabel('y [m]');
legend('Location', 'best');
title('Transmitter, Receiver, and Target Scenario');
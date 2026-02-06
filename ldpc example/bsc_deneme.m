clear;
close all;
clc;

numBits = 100000;
p = 0.1;

input_data = randi([0 1], numBits, 1);

[output_data, errors] = bsc(input_data, p);

numErrors = sum(errors); % Total number of errors
errorRate = numErrors / numBits; % Bit error rate

% Display results
fprintf('Number of errors: %d\n', numErrors);
fprintf('Bit error rate: %.2f%%\n', errorRate * 100);



numBits = 100000;
p = 0.1;

data = randi([0 1], numBits, 1);
input_data = repmat(data, 1, 3);
input_data = input_data.';
input_data = input_data(:);



[output_data, errors] = bsc(input_data, p);

output_data = reshape(output_data, 3, numBits);
output_data = output_data.';
output_data = sum(output_data, 2);
output_data(output_data < 2) = 0;
output_data(output_data > 0) = 1;

numErrors = sum(data ~= output_data); % Total number of errors
errorRate = numErrors / numBits; % Bit error rate

% Display results
fprintf('Number of errors: %d\n', numErrors);
fprintf('Bit error rate: %.2f%%\n', errorRate * 100);
% % Define two small matrices
% A = [1 2; 3 4];
% B = [5 6; 7 8];

% Create sample data
N = 50; M = 70;
A = rand(N, M);
B = rand(N, M);

% % Manual Circular Convolution
% result_manual = circular_conv2d(A, B);
% 
% % FFT-based Circular Convolution for verification
% result_fft = real(ifft2(fft2(A) .* fft2(B)));
% 
% disp('Manual Result:');
% disp(result_manual);
% disp('FFT Verification:');
% disp(result_fft);

% FT element-wise multiplication
res1 = sfft(A .* B);

% DD circular convolution
% res2 = circular_conv2d(sfft(A), sfft(B));
tic
res2 = circular_conv2d(sfft(A), sfft(B));
toc
abs(res1 ./ res2);

figure;
mesh(abs(res1));

figure;
mesh(abs(res2));
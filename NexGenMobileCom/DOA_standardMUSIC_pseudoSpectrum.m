% DOA 估計：MUSIC pseudo-spectrum
% 10 根天線 ULA、4 個訊號源、50 個 snapshot、SNR = 10 dB，由雜訊子空間計算 pseudo-spectrum 並找出前 4 個峰值
% 輸出：pseudo-spectrum，並印出估計角度與 MSE
clear;
clc;

M = 10;  
K = 4;  
t = 50;  
theta_true = [-43, -28, 17, 45]; 
angles = -90:90;  
% lambda = 1;    
% d = lambda/2; 
SNR_dB = 10;
SNR = 10^(SNR_dB / 10);

% steer_vec = @(theta) exp(1j*2*pi*d*(0:M-1)'*sin(theta*pi/180)/lambda);
steer_vec = @(theta) exp(1j*pi*(0:M-1)'*sind(theta));

A = zeros(M, K);
for k = 1:K
    A(:,k) = steer_vec(theta_true(k));
end

Rxx = zeros(M, M);
for i = 1:t
    s = (1/sqrt(2)) * (randn(K, 1) + 1j*randn(K, 1)); 
    signal = A * s;
    noise_power = 1 / SNR;
    N = sqrt(noise_power/2)*(randn(M, 1) + 1j*randn(M, 1)); 
    x = signal + N;
    Rxx = Rxx + x * x';
end
Rxx = Rxx / t;

% 特徵分解，取得雜訊子空間 Vn
[V, D] = eig(Rxx);
[eig_vals, idx] = sort(diag(D), 'descend');
V = V(:, idx);
Vn = V(:, end-(M-K)+1:end);  % 雜訊子空間

Pesudo_Spectrum = zeros(size(angles));
for i = 1:length(angles)
    a = steer_vec(angles(i)); % 每一個搜索角度，都會產生一個a向量(M*1)
    Pesudo_Spectrum(i) = 1 / (a' * Vn * Vn' * a);
end

% 尋找前 K 個主峰
[peaks, locs] = findpeaks(abs(Pesudo_Spectrum), 'SortStr', 'descend');
theta_esti = sort(angles(locs(1:K)));
disp(theta_esti);

MSE = sum((theta_true - theta_esti).^2) / K;
disp(MSE);

figure;
plot(angles, abs(Pesudo_Spectrum));
xlabel('DOA (\theta)');
ylabel('Pseudo Spectrum S(\theta)');
title('MUSIC Spectrum');
grid on;

figure;
semilogy(SNR_dB, MSE, '-o');
xlabel('SNR (dB)');
ylabel('MSE');
grid on;

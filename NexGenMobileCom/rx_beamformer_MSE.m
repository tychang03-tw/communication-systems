% 接收端波束成形：MSE 比較
% 與 rx_beamformer_beamPattern.m 相同情境，掃描 SNR = -20～50 dB；高 SNR 時 MVDR 改用訊號子空間處理 ill-conditioned 共變異矩陣
% 輸出：Delay-and-Sum、MMSE、MVDR 的 MSE vs. SNR
clear; clc;

M = 20; 
K = 4; 
Nt = 50000; 
theta = [0, -40, -45, 38];
theta_d = 0;

% 簽名向量
steerVec = @(theta) exp(1j*pi*(0:M-1)'*sind(theta));

% 簽名矩陣 A
A = zeros(M, K);
for k = 1:K
    A(:,k) = steerVec(theta(k));
end
ad = steerVec(theta_d); % 目標方向

% delay-and-sum beamformer（與 SNR 無關）
w_ds = (1/sqrt(M))*ad / norm(ad);

C = zeros(M, K);
for k = 1:K
    C(:,k) = steerVec(theta(k));
end
f = zeros(K, 1); f(1) = 1;

SNR_dB_range = -20:3:50;
MSE_ds = zeros(size(SNR_dB_range));
MSE_mmse = zeros(size(SNR_dB_range));
MSE_mvdr = zeros(size(SNR_dB_range));

for idx = 1:length(SNR_dB_range)
    SNR_dB = SNR_dB_range(idx);
    SNR = 10^(SNR_dB/10);
    noise_power = 1 / SNR;

    Ry = zeros(M, M);
    pd = zeros(M, 1);

    for i = 1:Nt
        s = (1/sqrt(2)) * (randn(K,1) + 1j * randn(K,1));
        td = s(1);  % 目標信號
        signal = A * s;
        N = sqrt(noise_power/2) * (randn(M,1) + 1j * randn(M,1));
        y = signal + N;
        Ry = Ry + y * y';
        pd = pd + y .* conj(td);
    end
    Ry = Ry / Nt;
    pd = pd / Nt;

    % MMSE beamformer
    w_mmse = inv(Ry) * pd;

    % MVDR beamformer 
    if SNR_dB >= 13
        % 高 SNR 時：子空間特徵值近似處理 ill-conditioned
        [U, D] = eig(Ry);
        [ev, idx_ev] = sort(diag(D), 'descend');
        Us = U(:, idx_ev(1:K));
        Ds = diag(ev(1:K));
        Ry_ill = Us * inv(Ds) * Us';
        w_mvdr = Ry_ill * C * inv(C' * Ry_ill * C) * f;
    else
        % 低 SNR 直接反矩陣
        w_mvdr = inv(Ry) * C * inv(C' * inv(Ry) * C) * f;
    end

    % MSE
    mse_ds = 0;
    mse_mmse = 0;
    mse_mvdr = 0;
    for i = 1:Nt
        s = (1/sqrt(2)) * (randn(K,1) + 1j * randn(K,1));
        td = s(1);
        signal = A * s;
        N = sqrt(noise_power/2) * (randn(M,1) + 1j * randn(M,1));
        y = signal + N;

        y_ds = w_ds' * y;
        y_mmse = w_mmse' * y;
        y_mvdr = w_mvdr' * y;

        mse_ds = mse_ds + abs(y_ds - td)^2;
        mse_mmse = mse_mmse + abs(y_mmse - td)^2;
        mse_mvdr = mse_mvdr + abs(y_mvdr - td)^2;
    end

    MSE_ds(idx) = mse_ds / Nt;
    MSE_mmse(idx) = mse_mmse / Nt;
    MSE_mvdr(idx) = mse_mvdr / Nt;
end

figure;
semilogy(SNR_dB_range, MSE_ds, 'k--x',...
         SNR_dB_range, MSE_mmse, 'r-.o',...
         SNR_dB_range, MSE_mvdr, 'b-s');
xlabel('SNR');
ylabel('MSE');
legend('Delay-and-Sum', 'MMSE', 'MVDR');
grid on;


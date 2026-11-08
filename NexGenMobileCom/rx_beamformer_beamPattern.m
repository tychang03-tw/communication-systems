% 接收端波束成形：波束場型
% 20 根天線 ULA，目標在 0 度、干擾在 -40 / -45 / 38 度，SNR = 30 dB，計算 Delay-and-Sum、MMSE、MVDR 的權重
% 輸出：三種波束成形器的波束場型
clear; clc;

M = 20;
K = 4;
Nt = 50000;
theta = [0, -40, -45, 38];
theta_d = 0;
SNR_dB = 30;
SNR = 10^(SNR_dB/10);

% 簽名向量
steerVec = @(theta) exp(1j*pi*(0:M-1)'*sind(theta));

% 簽名矩陣
A = zeros(M, K);
for k = 1:K
    A(:,k) = steerVec(theta(k));
end
ad = steerVec(theta_d);

Ry = zeros(M, M);
pd = zeros(M, 1);
for i = 1:Nt
    s = (1/sqrt(2)) * (randn(K,1) + 1j * randn(K,1));
    td = s(1,:);
    signal = A * s;
    noise_power = var(s) / SNR;
    N = sqrt(noise_power/2) * (randn(M,1) + 1j * randn(M,1));
    y = signal + N;

    Ry = Ry + y * y';
    pd = pd + y .* conj(td);
end
Ry = Ry / Nt;
pd = pd / Nt;

% delay and sum beamformer
w_ds = (1/sqrt(M))*ad / norm(ad);

% MMSE beamformer
w_mmse = inv(Ry)*pd;

% MVDR beamformer
C = zeros(M, K);
for k = 1:K
    C(:,k) = steerVec(theta(k));
end
f = zeros(K, 1); f(1) = 1;
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

% beamPattern
angles = -90:0.01:90;
for i = 1:length(angles)
    a_theta = steerVec(angles(i));
    Bp_ds(i) = abs(w_ds' * a_theta);
    Bp_mmse(i) = abs(w_mmse' * a_theta);
    Bp_mvdr(i) = abs(w_mvdr' * a_theta);
end
figure;
semilogy(angles,Bp_mvdr,'b-',angles,Bp_mmse, 'k--',angles,Bp_ds, 'r--');
hold on;
for i = 2:length(theta)
    xline(theta(i), 'r');
end
xline(theta(1), 'b');
xlabel('(\theta) (degree)');
ylabel('Beampattern');
xlim([-90 90]);
ylim([1e-8 10]);
legend('MVDR','MMSE','Delay and Sum');
grid on;

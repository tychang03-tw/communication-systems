% DOA 估計：ESPRIT 單次估計
% 10 根天線、SNR = 3 dB，利用兩個子陣列訊號子空間之間的 rotational invariance 求出 DOA
% 輸出：印出估計角度與 MSE
clear;
clc;

M = 10;
K = 4;
t = 50;
theta_true = [43, 27, -17,-28] ;

SNR_dB = 3;
SNR = 10^(SNR_dB / 10);

steer_vec = @(theta) exp(1j*pi*(0:M-1)'*sind(theta));

A = zeros(M, K);
for k = 1:K
    A(:,k) = steer_vec(theta_true(k));
end

Rxx = zeros(M, M);
for i = 1:t
    s = (1/sqrt(2)) * (randn(K,1) + 1j * randn(K,1));
    signal = A * s;
    noise_power = 1 / SNR;
    N = sqrt(noise_power/2) * (randn(M,1) + 1j * randn(M,1));
    x = signal + N;
    Rxx = Rxx + x * x';
end
Rxx = Rxx / t;

[V, D] = eig(Rxx);
[eig_vals, idx] = sort(diag(D), 'descend');
V = V(:, idx);
Vs = V(:, 1:K);

J1 = [eye(M-1), zeros(M-1,1)];
J2 = [zeros(M-1,1), eye(M-1)];

V1 = J1 * Vs;
V2 = J2 * Vs;

[T , psi] = eig(pinv(V1) * V2);
psi_d = diag(psi);
theta_esti = asind(angle(psi_d) / pi).';
disp(theta_esti);

MSE = sum(abs(sort(theta_true) - sort(theta_esti)).^2) / K;
disp(MSE);



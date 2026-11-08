% DOA 估計：Root-MUSIC 單次估計
% 10 根天線、SNR = -3 dB，由雜訊投影矩陣建立多項式並求根，取最靠近單位圓且角度不重複的 4 個根換算 DOA
% 輸出：根在複數平面上的分布，並印出估計角度與 MSE
clear;
clc;

M = 10;                    
K = 4;           
t = 50;                   
theta_true = [43, 27, -17,-28] ; 
SNR_dB = -3;
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

% 特徵值分解
[V, D] = eig(Rxx);
[eig_vals, idx] = sort(diag(D), 'descend');
V = V(:, idx);
Vn = V(:, K+1:end);  % 雜訊子空間

% 雜訊投影矩陣
P = Vn * Vn';

% 多項式係數
poly_coeff = zeros(2*M - 1, 1);
for m = 1:M
    for n = 1:M
        poly_coeff(m - n + M) = poly_coeff(m - n + M) + (P(m, n));
    end
end
poly_coeff = flipud(poly_coeff);  % 調整次數順序從高次到低次

% 求根並畫圖
roots_poly = roots(poly_coeff);
figure;
plot(real(roots_poly), imag(roots_poly), 'ko'); % 所有根
hold on;
theta = linspace(0, 2*pi, 1000);
plot(cos(theta), sin(theta), 'r--'); % 單位圓
axis equal;
grid on;
xlabel('Re(z)');
ylabel('Im(z)');
xlim([-1.5 1.5]);
ylim([-1.5 1.5]);

% 根據模長接近1排序
[~, sort_idx] = sort(abs(abs(roots_poly) - 1));

% 只挑選角度唯一的根
selected_roots = [];
selected_angles = [];

for i = 1:length(sort_idx)
    root = roots_poly(sort_idx(i));
    angle_root = angle(root);

    is_duplicate = any(abs(mod(angle_root - selected_angles + pi, 2*pi) - pi) < 1e-2);
    
    if ~is_duplicate
        selected_roots = [selected_roots; root];
        selected_angles = [selected_angles; angle_root];
    end

    if length(selected_roots) >= K
        break;
    end
end

% 畫出選中的根
plot(real(selected_roots), imag(selected_roots), 'bx');

% DOA估計
theta_esti = -asind(angle(selected_roots) / pi).';
disp(theta_esti); 

% MSE計算
[~, idx_true] = sort(theta_true);
[~, idx_est] = sort(theta_esti);
MSE = sum(abs(theta_true(idx_true) - theta_esti(idx_est)).^2) / K;
disp(MSE);



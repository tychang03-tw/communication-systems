% DOA 估計：MUSIC / Root-MUSIC / ESPRIT 比較
% 三種方法在每次 trial 使用同一組接收資料，天線數 M = 10 / 15，每個 SNR 跑 500 次 trial
% 輸出：三種方法在不同天線數下的 MSE vs. SNR
clear;
clc;
rng(0); % 固定亂數種子，結果可重現

K = 4;                              % 訊號源數
t = 50;                             % snapshot 數
theta_true = [43, 27, -17, -28];    % 真實 DOA (degree)
SNR_dB_range = -4:2:14;
num_trials = 500;
M_array = [10, 15];                 % 天線數

angles = -90:0.01:90;               % MUSIC 搜尋網格

methods = {'MUSIC', 'root-MUSIC', 'ESPRIT'};
MSE_all = zeros(length(methods), length(M_array), length(SNR_dB_range));

for m_idx = 1:length(M_array)
    M = M_array(m_idx);
    steer_vec = @(theta) exp(1j*pi*(0:M-1)'*sind(theta)); % d = lambda/2

    A = steer_vec(theta_true);      % M x K 簽名矩陣
    A_grid = steer_vec(angles);     % M x length(angles)，MUSIC 用，先算好加速

    J1 = [eye(M-1), zeros(M-1,1)];  % ESPRIT 子陣列選擇矩陣
    J2 = [zeros(M-1,1), eye(M-1)];

    for snr_idx = 1:length(SNR_dB_range)
        SNR = 10^(SNR_dB_range(snr_idx) / 10);
        noise_power = 1 / SNR;
        se = zeros(length(methods), num_trials); % 每次 trial 的 squared error

        for trial = 1:num_trials
            % 接收資料與樣本共變異矩陣（三種方法共用）----------------
            S = (1/sqrt(2)) * (randn(K,t) + 1j*randn(K,t));
            N = sqrt(noise_power/2) * (randn(M,t) + 1j*randn(M,t));
            X = A*S + N;
            Rxx = (X*X') / t;

            [V, D] = eig(Rxx);
            [~, idx] = sort(real(diag(D)), 'descend');
            V = V(:, idx);
            Vs = V(:, 1:K);         % 訊號子空間
            Vn = V(:, K+1:end);     % 雜訊子空間

            % standard MUSIC -----------------------------------------
            P_music = 1 ./ sum(abs(Vn' * A_grid).^2, 1);
            theta_music = find_K_peaks(P_music, angles, K);

            % root-MUSIC ---------------------------------------------
            P = Vn * Vn';
            poly_coeff = zeros(2*M-1, 1);
            for m = 1:M
                for n = 1:M
                    poly_coeff(m-n+M) = poly_coeff(m-n+M) + P(m,n);
                end
            end
            poly_coeff = flipud(poly_coeff);
            roots_poly = roots(poly_coeff);
            [~, sort_idx] = sort(abs(abs(roots_poly) - 1)); % 越靠近單位圓越前面
            selected_angles = [];
            for i = 1:length(sort_idx)
                ang = angle(roots_poly(sort_idx(i)));
                % 根以 z 與 1/z* 成對出現（角度相同），只取角度不重複的
                if ~any(abs(mod(ang - selected_angles + pi, 2*pi) - pi) < 1e-2)
                    selected_angles = [selected_angles; ang];
                end
                if length(selected_angles) >= K
                    break;
                end
            end
            theta_root = -asind(selected_angles / pi).';

            % ESPRIT -------------------------------------------------
            psi = eig(pinv(J1*Vs) * (J2*Vs));
            theta_esprit = real(asind(angle(psi) / pi)).';

            % MSE ----------------------------------------------------
            ref = sort(theta_true);
            se(1, trial) = mean((ref - sort(theta_music)).^2);
            se(2, trial) = mean((ref - sort(theta_root)).^2);
            se(3, trial) = mean((ref - sort(theta_esprit)).^2);
        end
        MSE_all(:, m_idx, snr_idx) = mean(se, 2);
        fprintf('M = %d, SNR = %3d dB done\n', M, SNR_dB_range(snr_idx));
    end
end

% 作圖 ------------------------------------------------------------------
figure;
% 方法用顏色 + marker 區分；MUSIC 與 root-MUSIC 結果幾乎重疊，
% 所以 MUSIC 用大的空心方形、root-MUSIC 用小的圓，兩者都看得到
colors  = {'k', 'r', 'b'};              % MUSIC / root-MUSIC / ESPRIT
markers = {'s', 'o', '^'};
lines   = {'--', '-'};                  % 天線數：M=10 / M=15
legend_str = {};
for m_idx = 1:length(M_array)
    for k = 1:length(methods)
        h = semilogy(SNR_dB_range, squeeze(MSE_all(k, m_idx, :)), ...
            [colors{k} lines{m_idx} markers{k}]);
        hold on;
        legend_str{end+1} = sprintf('%s, M=%d', methods{k}, M_array(m_idx));
    end
end
hold off;
xlabel('SNR (dB)');
ylabel('MSE');
title(sprintf('DOA Estimation Comparison (K=%d, snapshots=%d)', K, t));
legend(legend_str, 'Location', 'northeast');
grid on;

% -----------------------------------------------------------------------
function theta_est = find_K_peaks(spectrum, angles, K)
% 找出 pseudo-spectrum 中最大的 K 個局部峰值對應的角度
% （不用 findpeaks，避免依賴 Signal Processing Toolbox）
    s = spectrum(:).';
    is_peak = [false, s(2:end-1) > s(1:end-2) & s(2:end-1) >= s(3:end), false];
    peak_idx = find(is_peak);
    [~, order] = sort(s(peak_idx), 'descend');
    peak_idx = peak_idx(order);
    if length(peak_idx) < K
        % 低 SNR 時峰值可能不足 K 個（兩個訊號源被合併），以最強峰值補齊
        if isempty(peak_idx)
            [~, peak_idx] = max(s);
        end
        peak_idx = [peak_idx, repmat(peak_idx(1), 1, K - length(peak_idx))];
    end
    theta_est = angles(peak_idx(1:K));
end

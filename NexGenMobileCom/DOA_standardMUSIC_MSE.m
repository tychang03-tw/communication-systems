% DOA 估計：MUSIC 的 MSE
% 天線數 M = 10 / 15 / 20，每個 SNR 跑 100 次 trial，角度搜尋網格 0.01 度
% 輸出：MSE vs. SNR
clear;

angles = -90:0.01:90;
K = 4;
t = 50;
theta_true = [-43, -28, 17, 45];

SNR_dB_range = -4:2:14;
num_trials = 100;

M_array = [10, 15, 20];
MSE_all = zeros(length(M_array), length(SNR_dB_range));

for m_idx = 1:length(M_array)
    M = M_array(m_idx);
    steer_vec = @(theta) exp(1j*pi*(0:M-1)'*sind(theta));

    for snr_idx = 1:length(SNR_dB_range)
        SNR_dB = SNR_dB_range(snr_idx);
        SNR = 10^(SNR_dB / 10);
        mse_temp = zeros(num_trials, 1);

        for trial = 1:num_trials
            clc;
            disp([M,snr_idx,trial]);
            
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
            Vn = V(:, K+1:end);

            Pesudo_Spectrum = zeros(size(angles));
            for i = 1:length(angles)
                a = steer_vec(angles(i));
                Pesudo_Spectrum(i) = 1 / (a' * Vn * Vn' * a);
            end

            [peaks, locs] = findpeaks(abs(Pesudo_Spectrum), 'SortStr', 'descend');
            theta_esti = angles(locs(1:K));

            [~, idx_true] = sort(theta_true);
            [~, idx_est] = sort(theta_esti);
            mse_temp(trial) = sum(abs(theta_true(idx_true) - theta_esti(idx_est)).^2) / K;
        end
        MSE_all(m_idx, snr_idx) = mean(mse_temp);
    end
end

figure;
markers = {'k--x', 'k-.o','k-s'};
for m_idx = 1:length(M_array)
    semilogy(SNR_dB_range, MSE_all(m_idx,:), markers{m_idx});
    hold on;
end
xlabel('SNR (dB)');
ylabel('MSE');
legend('standard MUSIC M=10','standard MUSIC M=15','standard MUSIC M=20')
grid on;



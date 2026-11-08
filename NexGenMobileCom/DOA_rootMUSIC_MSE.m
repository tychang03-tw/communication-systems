% DOA 估計：Root-MUSIC 的 MSE
% 天線數 M = 10 / 15，每個 SNR 跑 500 次 trial
% 輸出：MSE vs. SNR
clear;


K = 4;
t = 50;
theta_true = [43, 27, -17, -28];

SNR_dB_range = -4:2:14;
num_trials = 500;

M_array = [10, 15];
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

            P = Vn * Vn';

            poly_coeff = zeros(2*M - 1, 1);
            for m = 1:M
                for n = 1:M
                    poly_coeff(m - n + M) = poly_coeff(m - n + M) + (P(m, n));
                end
            end
            poly_coeff = flipud(poly_coeff);

            % 求根
            roots_poly = roots(poly_coeff);

            % 根據絕對值接近1排序
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

            % DOA估計
            theta_esti = -asind(angle(selected_roots) / pi).';
            

            % MSE計算
            [~, idx_true] = sort(theta_true);
            [~, idx_est] = sort(theta_esti);
            mse_temp(trial) = sum(abs(theta_true(idx_true) - theta_esti(idx_est)).^2) / K;
        end
        MSE_all(m_idx, snr_idx) = mean(mse_temp);
    end
end

figure;
markers = {'k--x', 'k-.o'};
for m_idx = 1:length(M_array)
    semilogy(SNR_dB_range, MSE_all(m_idx,:), markers{m_idx});
    hold on;
end
xlabel('SNR (dB)');
ylabel('MSE');
legend('root MUSIC M=10','root MUSIC M=15')
grid on;


% DOA 估計：ESPRIT 的 MSE
% 天線數 M = 10 / 15，每個 SNR 跑 500 次 trial
% 輸出：MSE vs. SNR
clear;

K = 4;
t = 50;
theta_true = [43, -28, -17, 27];
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
            mse_temp(trial) = sum(abs(sort(theta_true) - sort(theta_esti)).^2) / K;
        end
        MSE_all(m_idx, snr_idx) = mean(mse_temp);
    end
end

figure;
markers = {'r--x', 'r-.o', 'r-s'};
for m_idx = 1:length(M_array)
    semilogy(SNR_dB_range, MSE_all(m_idx,:), markers{m_idx});
    hold on;
end
xlabel('SNR (dB)');
ylabel('MSE');
legend('ESPRIT M=10','ESPRIT M=15')
grid on;


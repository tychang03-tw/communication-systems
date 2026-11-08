% MIMO 偵測與空間調變
% 5x8 空間多工比較 ZF、MMSE、V-BLAST（ordered SIC）；4x8 空間調變以 2 bit 選天線 + 3 bit 8PSK，接收端先估天線再做 MRC
% 輸出：BER vs. SNR
clear;
clc;

% BPSK
M = 5;
N = 8;
BPSK_s = 2 * randi([0 1], M, 1) - 1;
times = 100000;

% SM
Nt = 4;
Nr = 8;
M_psk = 8;
k_sm = log2(Nt) + log2(M_psk); % bits=5

sm = randi([0 1], k_sm, 1);
ant_idx = bi2de(sm(1:2)', 'left-msb') + 1;  % 2 bits → 天線
sym_idx = bi2de(sm(3:end)', 'left-msb');  % 3 bits → 0~7
MPSK_s = exp(1j * 2 * pi * sym_idx / M_psk);  % 8PSK

x_sm = zeros(Nt, 1);
x_sm(ant_idx) = MPSK_s;

SNR = -10:2:20;
N0 = 1./((10.^(SNR/10)));
for j = 1:times
    clc
    disp(j);

    H = randn(N, M);
    r = H * (BPSK_s./sqrt(M));

    H_sm = randn(Nr, Nt);
    r_sm = H_sm * x_sm;

    for i = 1:length(SNR)
        noise = sqrt(N0(i)/2) * (randn(size(r)) + 1i * randn(size(r)));
        y = r + noise;

        % Zero-Forcing ------------------------
        s_hat1 = sign(real(pinv(H) * y));
        BER(j, i) = sum(s_hat1 ~= BPSK_s) / M;

        % V-BLAST -----------------------------
        s_hat2 = zeros(M, 1);
        H2 = H;
        y2 = y;
        index = 1:M;

        for k = 1:M
            Q = pinv(H2);
            rowNorms = vecnorm(Q, 2, 2);  % 計算列的範數
            [~, minIdx] = min(rowNorms);  % 找到最小範數的列
            q = Q(minIdx, :);  % 提取對應列向量
            s_hat2(index(minIdx)) = sign(real(conj(q) * y2));

            % 更新 y2，去除已解的部分
            y2 = y2 - (s_hat2(index(minIdx))./sqrt(M)) * H2(:,minIdx);
            H2(:, minIdx) = [];  % 刪除已解的列
            index(minIdx) = [];  % 更新索引
        end
        BER2(j, i) = sum(s_hat2 ~= (BPSK_s)) / M;

        % MMSE -------------------------------
        p(i) =  M./(10.^(SNR(i)/10));
        W_MMSE = inv(H * H' + p(i) * eye(N)) * H;
        s_hat3 = sign(real(W_MMSE' * y));
        BER3(j, i) = sum(s_hat3 ~= BPSK_s) / M;

        % SM ---------------------------------
        noise_sm = sqrt(N0(i)/2) * (randn(Nr,1) + 1j*randn(Nr,1));
        y_sm = r_sm + noise_sm;

        p = zeros(Nt, 1);
        for l = 1:Nt
            h_l = H_sm(:, l);
            p(l) = abs(h_l' * y_sm) / norm(h_l);
        end
        [~, est_ant] = max(p);

        h_hat = H_sm(:, est_ant);
        s_hat = (h_hat' * y_sm) / (norm(h_hat)^2);  % MRC

        constellation = exp(1j * 2 * pi * (0:M_psk-1) / M_psk);

        [~, est_sym] = min(abs(s_hat - constellation));

        rec_bits = [de2bi(est_ant - 1, 2, 'left-msb') de2bi(est_sym - 1, 3, 'left-msb')];
        BER4(j, i) = sum(sm' ~= rec_bits) / k_sm;

    end
end

BER_av = sum(BER) / times;
BER2_av = sum(BER2) / times;
BER3_av = sum(BER3) / times;
BER4_av = sum(BER4) / times;

% 繪製結果
semilogy(SNR , BER_av, 'k--x', SNR , BER2_av, 'b-o', SNR , BER3_av, 'r--s', SNR , BER4_av, 'g--p');
xlabel('SNR');
ylabel('BER');
% xlim([-10 15]);
ylim([1e-6 1]);
legend('Zero Forcing', 'V-BLAST', 'MMSE', 'Spatial Modulation');
grid on;

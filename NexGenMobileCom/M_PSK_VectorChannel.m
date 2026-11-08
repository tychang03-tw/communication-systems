% M-PSK 錯誤率（向量通道）
% BPSK / QPSK / 8PSK / 16PSK 用 pskmod 調變、awgn 加雜訊、pskdemod 解調
% 輸出：錯誤率 vs. SNR
clear; clc;
Mtype = [2 4 8 16];
colors = {'r-o','g-^','b-.','m--'};
figure;
for i = 1:4
    M = Mtype(i);

    % k = log2(M);                       
    % bits = randi([0 1],1,10*k);     
    % s = bi2de(reshape(bits,k,[]).','left-msb');
    % x = pskmod(s,M,0);

    s = randi([0 M-1],1,10000);   
    x = pskmod(s,M,0);

    BER = zeros(1,21);
    for SNR = 0:20
        y = awgn(x,SNR); % awgn(signal,Es/N0)，Es=1(0dB)的時候對應到要加的雜訊
        s_hat = pskdemod(y,M,0);
        [~,ratio] = biterr(s,s_hat);
        BER(SNR+1) = ratio;
    end
    semilogy(0:20,BER,colors{i});
    hold on;
end
ylim([1e-7 1])
xlabel('SNR (dB)');
ylabel('BER');
grid;
legend('BPSK','QPSK','8PSK','16PSK');

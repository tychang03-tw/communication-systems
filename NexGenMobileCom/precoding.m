% ISI 通道的等化與 precoding
% BPSK 通過 5-tap ISI 通道（另有 h2 可切換），比較不處理、DFE、傳送端 precoding（預先抵銷 ISI）
% 輸出：BER vs. Eb/N0
clear;
clc;

h1 = [1 -0.75 0.43 0.21 -0.11];
h2 = [-0.83 0.75 -0.66 0.48 -0.31 0.17];
h = h1;
KK = length(h);
M = 2;
binary_source = 100000;

EbN0 = -10:20;
Eb = 2;
N0 = Eb ./ (10.^(EbN0 / 10));
m = randi([0 M-1], 1, binary_source);
c = real(exp(1i * (2 * pi * m / M)));

% Precoding ---------------------------------------------------
s = zeros(1, binary_source);
for k = 1:binary_source
    if k < KK
        prec = sum(s(k-1:-1:1) .* h(2:k) / h(1));
    else
        prec = sum(s(k-1:-1:k-KK+1) .* h(2:KK) / h(1));
    end
    s(k) = c(k) / h(1) - prec;
end

% Channel (No Precoding Signal) ------------------------------
r = zeros(1, binary_source);
for k = 1:binary_source
    for l = 1:KK
        if k-l+1 > 0
            r(k) = r(k) + h(l) * c(k-l+1);
        end
    end
end

% Channel (Precoded Signal) ----------------------------------
r2 = zeros(1, binary_source);
for k = 1:binary_source
    for l = 1:KK
        if k-l+1 > 0
            r2(k) = r2(k) + h(l) * s(k-l+1);
        end
    end
end

numErrors = zeros(1, length(EbN0));
numErrors2 = zeros(1, length(EbN0));
numErrors3 = zeros(1, length(EbN0));

for i = 1:length(EbN0)
    noise = sqrt(N0(i)/2) * (randn(size(r)) + 1i * randn(size(r)));

    % No Eq -----------------------------------------------------
    r_noisy = r + noise;
    s_hat1 = sign(real(r_noisy));
    numErrors(i) = sum(s_hat1 ~= c);

    % DFE Eq.------------------------------------------------------
    s_hat2 = zeros(1, binary_source);
    for k = 1:binary_source
        ISI = 0;
        for l = 1:KK-1
            if k-l > 0
                ISI = ISI + h(l+1) * s_hat2(k-l);
            end
        end
        s_hat2(k) = sign(real((r_noisy(k) - ISI) / h(1)));
    end
    numErrors2(i) = sum(s_hat2 ~= c);

    % With Precoder ---------------------------------------------
    r2_noisy = r2 + noise;
    s_hat3 = sign(real(r2_noisy));
    numErrors3(i) = sum(s_hat3 ~= c);
    
end

BER = numErrors / binary_source;
BER2 = numErrors2 / binary_source;
BER3 = numErrors3 / binary_source;

semilogy(EbN0, BER, 'r', EbN0, BER2, 'b', EbN0, BER3, 'k--');
xlabel('E_b/N_0 (dB)');
ylabel('BER');
legend('Symbol-by-Symbol', 'DFE', 'Precoder', 'Location', 'southwest');
grid on;

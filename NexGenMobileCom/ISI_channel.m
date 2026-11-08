% ISI 通道：DFE / precoding / OFDM 比較
% BPSK 通過 5-tap ISI 通道，比較逐符號判決、DFE、傳送端 precoding，以及 OFDM（64 個子載波、CP = 16）
% 輸出：BER vs. SNR
clear;

h = [1 -0.52 0.41 -0.28 0.13];

KK = length(h);
M = 2;
N = 64;
CP = 16;
symbol_num = 100000;
binary_source = N * symbol_num;
SNR = -10:2:20;
m = randi([0 M-1], 1, binary_source);
c = real(exp(1i * (2 * pi * m / M)));

% OFDM mod. ---------------------------------------------------
c_parallel = reshape(c, N, []);
c_ofdm = ifft(c_parallel).*sqrt(N);
ofdm_CP = [c_ofdm(end-CP+1:end, :); c_ofdm];
ofdmTx = ofdm_CP(:);
ofdmTx = ofdmTx.';

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

% Channel (Non-Precoded Signal) ------------------------------
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

% Channel (OFDM) ----------------------------------------------
r3 = zeros(1, (N+CP)*symbol_num);
for k = 1:(N+CP)*symbol_num
    for l = 1:KK
        if k-l+1 > 0
            r3(k) = r3(k) + h(l) * ofdmTx(k-l+1);
        end
    end
end
% --------------------------------------------------------------

for i = 1:length(SNR)
    N0 = 1./ (10.^(SNR / 10));

    noise1 = sqrt(N0(i)/2) * (randn(size(r)));
    noise2 = sqrt(N0(i)/2) * (randn(size(r2)));
    noise3 = sqrt(N0(i)/2) * (randn(size(r3)));

    % No Eq -----------------------------------------------------
    r_noisy = r + noise1;
    s_hat1 = sign(real(r_noisy));
    numErrors(i) = sum(s_hat1 ~= c);

    % DFE ------------------------------------------------------
    s_hat2 = zeros(1, binary_source);
    for k = 1:binary_source
        ISI = 0;
        for l = 1:KK-1
            if k-l >= 1
                ISI = ISI + h(l+1) * s_hat2(k-l);
            end
        end
        s_hat2(k) = sign(real((r_noisy(k) - ISI) / h(1))); % Eq.
    end
    numErrors2(i) = sum(s_hat2 ~= c);

    % With Precoder ---------------------------------------------
    r2_noisy = r2 + noise2;
    s_hat3 = sign(real(r2_noisy));
    numErrors3(i) = sum(s_hat3 ~= c);

    % OFDM demod.------------------------------------------------
    r3_noisy = r3 + noise3;
    r4 = reshape(r3_noisy,N+CP,[]);
    rx = r4(CP+1:end,:);
    rxx = fft(rx)./sqrt(N);
    rxx = rxx(:).';
    s_hat4 = sign(real(rxx));
    numErrors4(i) = sum(s_hat4 ~= c);

    % BPSK BER --------------------------------------------------
    z(i) = sqrt(2 * 10.^(SNR(i)./10));
    Qf(i) = (1/sqrt(2*pi)) * integral(@(x) exp(-x.^2 / 2), z(i), inf);

end
BER = numErrors / binary_source;
BER2 = numErrors2 / binary_source;
BER3 = numErrors3 / binary_source;
BER4 = numErrors4 / binary_source;

semilogy(SNR, BER, 'r', SNR, BER2, 'b', SNR, BER3, 'k-', SNR , BER4 ,'k-^');
xlabel('SNR');
ylabel('BER');
ylim([1e-7 1]);
legend('BPSK Symbol-by-Symbol', 'DFE', 'Precoder','OFDM' ,'Location', 'southwest');
grid on;

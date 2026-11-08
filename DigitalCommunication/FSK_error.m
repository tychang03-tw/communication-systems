% FSK 經 matched filter 解調的 BER
% f1 = 1 kHz、f2 = 2 kHz，matched filter 取兩個載波之差，門檻為 0；比較 T_b = 10 ms 與 1 ms 下的模擬 BER 與理論 Q 函數
% 輸出：調變波形、BER vs. Eb/N0、matched filter 輸出波形
clear;
A = 0.8;
f1 = 1000;
f2 = 2000;
fs = 1e6;
T = 10^-2;
T2 = 10^-3;
num = 8;
m = randi([0 1],1,num);

t = 0:1/fs:T-1/fs;
t_total = 0:1/fs:(num*T)-1/fs;
t2 = 0:1/fs:T2-1/fs;
t2_total = 0:1/fs:(num*T2)-1/fs;

data = [];
for i = 1:num
    data = [data, m(i)*ones(1, length(t))];
end
fsk_signal = [];
for i = 1:num
    if m(i) == 1
        fsk_signal = [fsk_signal, A*sin(2*pi*f2*t)];
    else
        fsk_signal = [fsk_signal, A*sin(2*pi*f1*t)];
    end
end

data2 = [];
for i = 1:num
    data2 = [data2, m(i)*ones(1, length(t2))];
end
fsk_signal2 = [];
for i = 1:num
    if m(i) == 1
        fsk_signal2 = [fsk_signal2, A*sin(2*pi*f2*t2)];
    else
        fsk_signal2 = [fsk_signal2, A*sin(2*pi*f1*t2)];
    end
end
figure;
subplot(2,1,1);
plot(t_total, data,t_total, fsk_signal);
title('FSK T=1e-2')
xlim([0 8*T])
ylim([-1.5 1.5])
grid on;
subplot(2,1,2);
plot(t2_total, data2,t2_total, fsk_signal2);
title('FSK T=1e-3')
xlim([0 8*T2])
ylim([-1.5 1.5])
grid on;

EbN0 = 22:-3:-32; EbN0w = 10.^(EbN0./10);
Eb = (A^2)*length(t)/2; N0w = Eb./EbN0w; N0 = 10.*log10(N0w);
Eb2 = (A^2)*length(t2)/2; N0w2 = Eb2./EbN0w; N02 = 10.*log10(N0w2);
for i = 1:length(N02)
    %T1 = 5*10^-3 ---------------------------------------------------------
    z1(i) = sqrt((A.^2).*length(t)./(2.* N0w(i)));
    Q1(i) = (1/sqrt(2*pi)) * integral(@(x) exp(-x.^2 / 2), z1(i), 70);
    
    nois = sqrt((0.5.*10^((-EbN0(i)) ./ 10))*Eb);
    n = nois * randn(size(fsk_signal));
    v = fsk_signal + n;
    matchFilter = A *(sin(2 * pi * f2 * (T-t))-sin(2 * pi * f1 * (T-t)));
    y = conv(v, matchFilter);

    sampledSignal = [];
    t_samples =length(t):length(t):length(y);
    for k = 1:length(t_samples)
        sampledSignal = [sampledSignal, y(t_samples(k))];
    end
    detectedBits = [];
    decisionThreshold = 0; %(sum((A*sin(2*pi*f2*t).^2))-sum(A*sin(2*pi*f1*t).^2))/2
    detectedBits = sampledSignal > decisionThreshold;
    numErrors(i) = sum(m ~= detectedBits);

    %T2 = 5*10^-4 ---------------------------------------------------------
    z2(i) = sqrt((A.^2).*length(t2)./(2.* N0w2(i)));
    Q2(i) = (1/sqrt(2*pi)) * integral(@(x) exp(-x.^2 / 2), z2(i), 70);
    nois2 = sqrt((0.5.*10^((-EbN0(i)) ./ 10))*Eb2);
    n2 = nois2 * randn(size(fsk_signal2));
    v2 = fsk_signal2 + n2;
    matchFilter2 = A *(sin(2 * pi * f2 * (T2-t2))-sin(2 * pi * f1 * (T2-t2)));
    y2 = conv(v2, matchFilter2);

    sampledSignal2 = [];
    t_samples2 =length(t2):length(t2):length(y2);
    for kk = 1:length(t_samples2)
        sampledSignal2 = [sampledSignal2, y2(t_samples2(kk))];
    end
    detectedBits2 = [];
    decisionThreshold2 = 0; %(sum((A*sin(2*pi*f2*t).^2))-sum(A*sin(2*pi*f1*t).^2))/2;
    detectedBits2 = sampledSignal2 > decisionThreshold2;
    numErrors2(i) = sum(m ~= detectedBits2);

end
BER = numErrors / num;
BER2 = numErrors2 / num;

figure;
semilogy(EbN0,BER,'b-o',EbN0,Q1);
xlabel('E_b/N_0 (dB)')
ylabel('BER')
legend('實驗值 P_e ', '理論值 P_e', 'location', 'southwest')
xlim([-30 20])
ylim([1e-4 1])
title('A=0.8,T=1e-2')
grid on;
figure;
semilogy(EbN0,BER2,'m-o',EbN0,Q2);
xlabel('E_b/N_0 (dB)')
ylabel('BER')
legend('實驗值 P_e ', '理論值 P_e', 'location', 'southwest')
title('A=0.8,T=1e-3')
xlim([-30 20])
ylim([1e-4 1])
grid on;

figure;
subplot(411)
sgtitle('FSK Modulation')
plot(t_total, data)
title('原始數位波 T=1e-2')
ylim([-2 2])
grid on
subplot(412)
plot(t_total, y(1:length(t_total)),t_total, data);
title('經過匹配濾波器之傳輸信號 Eb/N0=3db')
grid on;
subplot(413)
plot(t2_total, data2)
title('原始數位波 T=1e-3')
ylim([-2 2])
grid on
subplot(414)
plot(t2_total, y2(1:length(t2_total)),t2_total, data2);
title('經過匹配濾波器之傳輸信號 Eb/N0=3db')
grid on;




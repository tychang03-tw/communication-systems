% ASK 經 matched filter 解調的 BER
% 在每個 bit 結束處取樣，以最佳門檻判決；比較 T_b = 10 ms 與 1 ms 下的模擬 BER 與理論 Q 函數
% 輸出：調變波形、BER vs. Eb/N0、matched filter 輸出波形
clear
clc
A = 0.8;
f = 2000;
fs = 1e5;
num = 8;

T1 = 10^-2;
m = randi([0 1],1,num);
t1 = 0:1/fs:T1-1/fs;
t1_total = 0:1/fs:(num*T1)-1/fs;

T2 = 10^-3;
t2 = 0:1/fs:T2-1/fs;
t2_total = 0:1/fs:(num*T2)-1/fs;

data1 = [];
for i = 1:num
    data1 = [data1, m(i)*ones(1, length(t1))];
end
ask_Tx1 = [];
for i = 1:num
    if m(i) == 1
        ask_Tx1 = [ask_Tx1, A*sin(2*pi*f*t1)];
    else
        ask_Tx1 = [ask_Tx1, zeros(1, length(t1))];
    end
end

data2 = [];
for i = 1:num
    data2 = [data2, m(i)*ones(1, length(t2))];
end
ask_Tx2 = [];
for i = 1:num
    if m(i) == 1
        ask_Tx2 = [ask_Tx2, A*sin(2*pi*f*t2)];
    else
        ask_Tx2 = [ask_Tx2, zeros(1, length(t2))];
    end
end
figure;
subplot(2,1,1);
plot(t1_total, data1,t1_total, ask_Tx1);
title('ASK T=1e-2')
xlim([0 8*T1])
ylim([-2 3])
grid on;
subplot(2,1,2);
plot(t2_total, data2,t2_total, ask_Tx2);
title('ASK T=1e-3')
xlim([0 8*T2])
ylim([-2 3])
grid on;

EbN0db = 16:-3:-35; EbN0w = 10.^(EbN0db./10);

Eb1 = (A^2)*length(t1)/2; N0w1 = Eb1./EbN0w; N01 = 10.*log10(N0w1);
Eb2 = (A^2)*length(t2)/2; N0w2 = Eb2./EbN0w; N02 = 10.*log10(N0w2); 
for i = 1:length(N01)
    %T1
    %理論值
    z1(i) = sqrt((A.^2).*length(t1)./(4.* N0w1(i)));
    Q1(i) = (1/sqrt(2*pi)) * integral(@(x) exp(-x.^2 / 2), z1(i), 70);
    %雜訊
    nois1 = sqrt((0.5.*10^((-EbN0db(i)) ./ 10))*Eb1);
    n1 = nois1 * randn(size(ask_Tx1));
    v1 = ask_Tx1 + n1;
    matchFilter1 = A * sin(2 * pi * f * (T1-t1));

    Rx1 = conv(v1, matchFilter1);
    %位元結束處取樣
    sampledSignal1 = [];
    t_samples1 =length(t1):length(t1):length(Rx1);
    for k = 1:length(t_samples1)
        sampledSignal1 = [sampledSignal1, Rx1(t_samples1(k))];
    end
    %比對
    detectedBits1 = [];
    kopt1 =sum((A*sin(2*pi*f*t1)).^2)/2;
    detectedBits1 = sampledSignal1 > kopt1;
    numErrors1(i) = sum(m ~= detectedBits1);

    %T2
    z2(i) = sqrt((A.^2).*length(t2)./(4.* N0w2(i)));
    Q2(i) = (1/sqrt(2*pi)) * integral(@(x) exp(-x.^2 / 2), z2(i), 70);

    nois2 = sqrt((0.5.*10^((-EbN0db(i)) ./ 10))*Eb2);
    n2 = nois2 * randn(size(ask_Tx2));
    v2 = ask_Tx2 + n2;

    matchFilter2 = A * sin(2 * pi * f * (T2-t2));
    Rx2 = conv(v2, matchFilter2);

    sampledSignal2 = [];
    t_samples2 =length(t2):length(t2):length(Rx2);
    for kk = 1:length(t_samples2)
        sampledSignal2 = [sampledSignal2, Rx2(t_samples2(kk))];
    end
    detectedBits2 = [];
    kopt2 =sum((A*sin(2*pi*f*t2)).^2)/2;
    detectedBits2 = sampledSignal2 > kopt2;
    numErrors2(i) = sum(m ~= detectedBits2);
end
BER1 = numErrors1 / num;
BER2 = numErrors2 / num;

figure;
semilogy(EbN0db,BER1,'b-o',EbN0db,Q1,'m-x');
xlabel('E_b/N_0 (dB)')
ylabel('BER')
legend('實驗值 P_e ', '理論值 P_e', 'location', 'southwest')
title('A=0.8,T=1e-2');
grid on;
figure;
semilogy(EbN0db,BER2,'m-x',EbN0db,Q2);
xlabel('E_b/N_0 (dB)')
ylabel('BER')
legend('實驗值 P_e ', '理論值 P_e', 'location', 'southwest')
title('A=0.8,T=1e-3');
grid on;

figure;
subplot(211)
plot(t1_total, data1)
ylim([-2 3])
grid on
subplot(212)
plot(t1_total, Rx1(1:length(t1_total)),t1_total,data1);
grid on;
hold off;



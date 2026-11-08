% 16-QAM 符號錯誤率
% passband 波形加 AWGN → correlator 取得 I/Q → 最小距離判決，與理論 3Q - (9/4)Q^2 比較
% 輸出：16-QAM 波形、SER vs. Es/N0、星座圖
clear;
clc;
format long
f = 10e6;
T = 1e-6;
fs = 1e8;
num = 1000;
t = 0 : 1/fs : T-1/fs;
t_total = 0 : 1/fs : (num*T)-1/fs;
a = 1e-3;
EsN0 = -6:2:20;
EsN0w= 10.^(EsN0./10);
Es = 10*a^2;
N0w = Es./EsN0w;

QAMw = zeros(1, length(t_total));
m = (sign(randn(1, num)) + 2*sign(randn(1, num))) + 1i*((sign(randn(1, num)) + 2*sign(randn(1, num))));
for i = 1:num
    QAMw(1,(i-1)*length(t)+1 : i*length(t))...
        = real(m(i)).*a.*sqrt(2/T).*cos(2*pi*f*t)+imag(m(i)).*a.*sqrt(2/T).*sin(2*pi*f*t);
end

E = (QAMw*QAMw')./num;
for k=1:length(EsN0)
    noise_std = sqrt((0.5.*10.^(-EsN0(k) ./ 10))*E);
    n = noise_std .* randn(size(QAMw));
    recieve_QAMw = QAMw + n;

    z(k) = sqrt(2*a^2./N0w(k));
    Qf(k) = 3.*((1/sqrt(2*pi)) * integral(@(x) exp(-x.^2 / 2), z(k), 100)) - (9/4).*((1/sqrt(2*pi)) * integral(@(x) exp(-x.^2 / 2), z(k), 100)).^2;

    I = zeros(num, length(t));
    Q = zeros(num, length(t));
    for i = 1:num
        r = recieve_QAMw(1,(i-1)*length(t)+1 : i*length(t));
        I(i,:) = sum(r .* sqrt(2/T) .* cos(2*pi*f*t)) * (1/fs);
        Q(i,:) = sum(r .* sqrt(2/T) .* sin(2*pi*f*t)) * (1/fs);
    end
    y = I(:,1) + 1i*Q(:,1);
    yx = [I(:,1) Q(:,1)];
    for j = 1:num %M APrule => '等事前機率' => MLrule => 'AWGN' => MDrule 歐幾里得距離
        r1(1,j) = norm(yx(j,:)-1e-3.*[-3 3]);
        r2(1,j) = norm(yx(j,:)-1e-3.*[-1 3]);
        r3(1,j) = norm(yx(j,:)-1e-3.*[1 3]);
        r4(1,j) = norm(yx(j,:)-1e-3.*[3 3]);
        r5(1,j) = norm(yx(j,:)-1e-3.*[-3 1]);
        r6(1,j) = norm(yx(j,:)-1e-3.*[-1 1]);
        r7(1,j) = norm(yx(j,:)-1e-3.*[1 1]);
        r8(1,j) = norm(yx(j,:)-1e-3.*[3 1]);
        r9(1,j) = norm(yx(j,:)-1e-3.*[-3 -1]);
        r10(1,j) = norm(yx(j,:)-1e-3.*[-1 -1]);
        r11(1,j) = norm(yx(j,:)-1e-3.*[1 -1]);
        r12(1,j) = norm(yx(j,:)-1e-3.*[3 -1]);
        r13(1,j) = norm(yx(j,:)-1e-3.*[-3 -3]);
        r14(1,j) = norm(yx(j,:)-1e-3.*[-1 -3]);
        r15(1,j) = norm(yx(j,:)-1e-3.*[1 -3]);
        r16(1,j) = norm(yx(j,:)-1e-3.*[3 -3]);
        [~, minIndex] = min([r1(1,j),r2(1,j),r3(1,j),r4(1,j),r5(1,j),r6(1,j),r7(1,j),r8(1,j),r9(1,j),r10(1,j),r11(1,j),r12(1,j),r13(1,j),r14(1,j),r15(1,j),r16(1,j)]);
        %找到最小的返回到第幾個case
        switch minIndex
            case 1
                rx(1,j) = (-3 + 3i);
            case 2
                rx(1,j) = (-1 + 3i);
            case 3
                rx(1,j) = (1 + 3i);
            case 4
                rx(1,j) = (3 + 3i);
            case 5
                rx(1,j) = (-3 + 1i);
            case 6
                rx(1,j) = (-1 + 1i);
            case 7
                rx(1,j) = (1 + 1i);
            case 8
                rx(1,j) = (3 + 1i);
            case 9
                rx(1,j) = (-3 - 1i);
            case 10
                rx(1,j) = (-1 - 1i);
            case 11
                rx(1,j) = (1 - 1i);
            case 12
                rx(1,j) = (3 - 1i);
            case 13
                rx(1,j) = (-3 - 3i);
            case 14
                rx(1,j) = (-1 - 3i);
            case 15
                rx(1,j) = (1 - 3i);
            case 16
                rx(1,j) = (3 - 3i);
        end
    end
numErrors(k) = sum(m ~= rx);
end
SER = numErrors / num;
figure;
plot(t_total, QAMw);
xlim([0 5*T]);
ylim([-7 7]);
grid on;
figure;
semilogy(EsN0,SER,'b-x',EsN0,Qf,'m-o');
xlabel('Es/N0(dB)');
ylabel('SER');
legend('實驗值 P_e ', '理論值 P_e', 'location', 'southwest')
title('16 QAM')
grid on;
figure;
scatter(real(y), imag(y), 'p');
title('Constellation');
grid on;



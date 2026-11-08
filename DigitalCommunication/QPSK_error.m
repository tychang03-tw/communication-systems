% QPSK 符號錯誤率
% passband 波形加 AWGN → correlator 取得 I/Q → 內積判決（等能量時等同最小距離判決），與理論 2Q - Q^2 比較
% 輸出：QPSK 波形、SER vs. Es/N0、星座圖
clear;
clc;
format long
f = 3e6;
T = 1e-6;
fs = 1e8;
num = 800;
t = 0 : 1/fs : T-1/fs;
t_total = 0 : 1/fs : (num*T)-1/fs;
Es = 1e-6;
EsN0 = -4:2:16;
QPSKw = zeros(1, length(t_total));
m = sign(randn(1,num)) + 1i * sign(randn(1,num));
for i = 1:num
    QPSKw(1,(i-1)*length(t)+1 : i*length(t))...
        = real(m(i)).*sqrt(Es/2).*sqrt(2/T).*cos(2*pi*f*t) + imag(m(i)).*sqrt(Es/2).*sqrt(2/T).*sin(2*pi*f*t);
end
E = (QPSKw*QPSKw')./num;
for k=1:length(EsN0)
    noise_std = sqrt((0.5.*10.^(-EsN0(k) ./ 10))*E);
    n = noise_std .* randn(size(QPSKw));
    recieve_QAMw = QPSKw + n;
    z(k) = sqrt(10.^(EsN0(k)./10));
    Qf(k) = 2*(1/sqrt(2*pi)) * integral(@(x) exp(-x.^2 / 2), z(k), 70) - ((1/sqrt(2*pi)) * integral(@(x) exp(-x.^2 / 2), z(k), 70)).^2;
    I = zeros(num, length(t));
    Q = zeros(num, length(t));
    for i = 1:num
        r = recieve_QAMw(1,(i-1)*length(t)+1 : i*length(t));
        I(i,:) = sum(r .* sqrt(2/T) .* cos(2*pi*f*t)) * (1/fs);
        Q(i,:) = sum(r .* sqrt(2/T) .* sin(2*pi*f*t)) * (1/fs);
    end
    y = I(:,1) + 1i*Q(:,1);
    yx = [I(:,1) Q(:,1)];
    for j = 1:num % MDrule-->E同--->內積rule 
        r1(1,j) = dot(yx(j,:), 1e-3.*[1 1] );
        r2(1,j) = dot(yx(j,:), 1e-3.*[-1 1]);
        r3(1,j) = dot(yx(j,:), 1e-3.*[-1 -1]);
        r4(1,j) = dot(yx(j,:), 1e-3.*[1 -1]);

        if r1(1,j)==max([r1(1,j), r2(1,j), r3(1,j), r4(1,j)])
            rx(1,j)= (1+1i) ;
        elseif r2(1,j)==max([r1(1,j), r2(1,j), r3(1,j), r4(1,j)])
            rx(1,j)= (-1+1i);
        elseif r3(1,j)==max([r1(1,j), r2(1,j), r3(1,j), r4(1,j)])
            rx(1,j)= (-1-1i) ;
        elseif r4(1,j)==max([r1(1,j), r2(1,j), r3(1,j), r4(1,j)])
            rx(1,j)= (1-1i);
        end
    end
    numErrors(k) = sum(m ~= rx);
end
SER = numErrors / num;
figure;
plot(t_total, QPSKw);
xlim([0 8*T]);
ylim([-3 3]);
grid on;
figure;
semilogy(EsN0,SER,'b-x',EsN0,Qf,'m-o');
xlabel('Es/N0(dB)');
ylabel('SER');
title('QPSK')
legend('實驗值 P_e ', '理論值 P_e', 'location', 'southwest')
grid on;
figure;
scatter(real(y), imag(y), 'p');
title('Constellation');
grid on;

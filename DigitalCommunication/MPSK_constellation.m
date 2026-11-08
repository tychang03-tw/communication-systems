% M-PSK 星座圖
% 8-PSK passband 波形加 AWGN（Es/N0 = 20 dB），以 correlator 取得 I/Q；改 M 可換成其他階數
% 輸出：接收端星座圖
clear;
clc;
f = 3e6;
T = 1e-6;
fs = 1e8;
num = 1e3;
t = 0 : 1/fs : T-1/fs;
t_total = 0 : 1/fs : (num*T)-1/fs;
Es = 1e-6;
EsN0 = 20;
QAMw = zeros(1, length(t_total));

M = 8; %8個相位
m = randi([0 M-1], 1, num); %產生0~7數字
c = exp(1i * (2 * pi * m / M)); %信號用極座標表示
for i = 1:num
    QAMw(1,(i-1)*length(t)+1 : i*length(t))...
        = real(c(i)).*sqrt(Es).*sqrt(2/T).*cos(2*pi*f*t)+imag(c(i)).*sqrt(Es).*sqrt(2/T).*sin(2*pi*f*t);
end
E = (QAMw*QAMw')./num;
noise = sqrt((0.5.*10^(-EsN0 ./ 10))*E);
n = noise .* randn(size(QAMw));
recieve_QAMw = QAMw + n ;

I = zeros(num,length(t));
Q = zeros(num,length(t));
for i = 1:num
    r = recieve_QAMw(1,(i-1)*length(t)+1 : i*length(t));
    I(i,:) = sum(r .* sqrt(2/T) .* cos(2*pi*f*t)).*(1/fs); 
    Q(i,:) = sum(r .* sqrt(2/T) .* sin(2*pi*f*t)).*(1/fs);  
end
y = I(:,1)*1 + Q(:,1)*1i ;

figure;
scatter(real(y),imag(y),'p')
title('constellation');
grid on;





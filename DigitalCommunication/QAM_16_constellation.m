% 16-QAM 星座圖
% 16-QAM passband 波形加 AWGN（Es/N0 = 20 dB），以 correlator 取得 I/Q
% 輸出：接收端星座圖
clear;
clc;
f = 3e6;
T = 1e-6;
fs = 1e8;
num = 1000;
t = 0 : 1/fs : T-1/fs;
t_total = 0 : 1/fs : (num*T)-1/fs;
Es = 1e-6;
EsN0 = 20;
QAMw = zeros(1, length(t_total));
m = (sign(randn(1, num)) + 2*sign(randn(1, num))) + 1i*(1*(sign(randn(1, num)) + 2*sign(randn(1, num))));

for i = 1:num
    QAMw(1,(i-1)*length(t)+1 : i*length(t))...
        = real(m(i)).*sqrt(Es).*sqrt(2/T).*cos(2*pi*f*t)+imag(m(i)).*sqrt(Es).*sqrt(2/T).*sin(2*pi*f*t);
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





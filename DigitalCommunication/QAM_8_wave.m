% 8-QAM 調變波形
% I 軸 4 個準位（±1、±3）、Q 軸 2 個準位（±1），產生 8 個隨機符號的 passband 波形
% 輸出：8-QAM 調變波形
clear
clc
f = 3e6;
T = 1e-6;
fs = 1e8;
num = 8;
t = 0 : 1/fs : T-1/fs;
t_total = 0 : 1/fs : (num*T)-1/fs;
a = 1e-3;

m = (sign(randn(1,num))+2*sign(randn(1,num))) + 1i * sign(randn(1,num));

QAMw = zeros(1, length(t_total));

for i = 1:num
    QAMw(1,(i-1)*length(t) +1 : i*length(t))...
       = real(m(i)) .* a/2.*sqrt(2/T).* cos(2*pi*f*t) + imag(m(i)) .* a/2.*sqrt(2/T).* sin(2*pi*f*t);
end
figure;
plot(t_total, QAMw);
title('QAM');
grid on;



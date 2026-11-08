% 訊號與系統複習：三角脈衝通過理想低通濾波器
% 三角脈衝（tau = 0.1）通過理想 LPF，以時域摺積求輸出
% 輸出：H(f)、x(t)、X(f)、y(t)
clear
clc
fs = 10;  
W = 6;
tau = 0.1; 
f = -30:1/fs:30;
t = -30:1/fs:30; 
tconv = -60:1/fs:60;
x = tripuls(t, 2*tau); 
xf = tau*sinc(tau*f).^2;

Hf = 3*rectpuls(f, W);
ht = 3*tau *(sinc(tau*t));
y = conv(ht,x)./fs;

figure;
subplot(221)
plot(f,Hf);
xlabel('Frequency (Hz)');
xlim([-10 10]);
ylim([-1 4])
title('H(f)');
grid on

subplot(222)
plot(t,x);
xlabel('Time (s)');
xlim([-1 1]);
ylim([-1 2])
title('x(t)');
grid on

subplot(223)
plot(f,xf);
xlabel('Frequency (Hz)');
xlim([-30 30]);
title('X(f)');
grid on

subplot(224)
plot(tconv,y);
xlabel('Time (s)');
xlim([-30 30]);
title('y(t)');
grid on

sgtitle('W=6,tau=0.1');

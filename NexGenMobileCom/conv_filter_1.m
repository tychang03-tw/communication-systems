% 訊號與系統複習：弦波通過理想低通濾波器
% x(t) = 2cos(2*pi*5t) 通過 H(f) = 3*rect(f/6)，以時域摺積求輸出
% 輸出：H(f)、x(t)、X(f)、y(t)
clear
clc
f0 = 5;
fs = 20;
W = 6;
t = -6:1/fs:6-1/fs;
tconv = -12:1/fs:12-2/fs;

x = 2*cos(2*pi*f0*t);
xf = fftshift(fft(x))./240;
f = linspace(-fs/2, fs/2, length(xf));

f_h = -10:1/fs:10-1/fs;
Hf = 3*rectpuls(f_h, W);
h = 3*W*sinc(W*t);

y = conv(h, x) ./ fs;

figure;
subplot(221)
plot(f_h,Hf);
xlabel('Frequency (Hz)');
ylim([-1 4])
title('H(f)');
grid on
subplot(222)
plot(t,x);
xlabel('Time (s)');
xlim([-3 3])
ylim([-3 3])
title('x(t)');
grid on
subplot(223)
plot(f, abs(xf)); 
xlabel('Frequency (Hz)');
ylim([0 1.2])
title('X(f)');
grid on
subplot(224)
plot(tconv,y);
xlabel('Time (s)');
xlim([-3 3])
ylim([-7 7])
title('y(t)');
grid on






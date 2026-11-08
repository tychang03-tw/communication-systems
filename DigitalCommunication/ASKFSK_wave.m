% ASK / FSK 調變波形
% 產生 10 個隨機 bit，分別做 ASK（on-off）與 FSK（f1 = 200 Hz、f2 = 600 Hz）調變
% 輸出：載波、bit stream、調變後波形
clc;
clear;
A = 1;
f = 500;
T = 5*10^-3;
fs = 1e6;
num = 10;
m = randi([0 1],1,num);
t = 0 : 1/fs : T-1/fs;
t_total = 0 : 1/fs : (num*T)-1/fs;

digital_signal = [];
for i = 1:num
    digital_signal = [digital_signal, m(i)*ones(1, length(t))];
end
ask_signal = [];
for i = 1:num
    if m(i) == 1
        ask_signal = [ask_signal, A*sin(2*pi*f*t)];
    else
        ask_signal = [ask_signal, zeros(1, length(t))];
    end
end
f1 = 200;
f2 = 600;
t = 0:1/fs:T-1/fs;
t_total = 0:1/fs:(num*T)-1/fs;

digital_signal = [];
for i = 1:num
    digital_signal = [digital_signal, m(i)*ones(1, length(t))];
end

fsk_signal = [];
for i = 1:num
    if m(i) == 1
        fsk_signal = [fsk_signal, A*sin(2*pi*f2*t)];
    else
        fsk_signal = [fsk_signal, A*sin(2*pi*f1*t)];
    end
end

figure;
subplot(311);
plot(t_total, A*sin(2*pi*f*t_total));
ylim([-1.5 1.5]);
title('ASK');
grid on;
subplot(312);
plot(t_total, digital_signal);
ylim([-1.2 2.2]);
grid on;
subplot(313);
plot(t_total, ask_signal);
ylim([-1.5 1.5]);
grid on;

figure;
subplot(411);
plot(t_total, A*sin(2*pi*f1*t_total));
ylim([-1.5 1.5]);
grid on;
title('FSK');
subplot(412);
plot(t_total, A*sin(2*pi*f2*t_total));
ylim([-1.5 1.5]);
grid on;
subplot(413);
plot(t_total, digital_signal);
ylim([-1.2 2.2]);
grid on;
subplot(414);
plot(t_total, fsk_signal);
ylim([-1.5 1.5]);
grid on;





% 16-QAM 功率頻譜密度（periodogram 法）
% 每組波形直接做 |FFT|^2，200 次平均；比較 T = 1 us 與 4 us 和理論 PSD
% 輸出：平均 PSD、所有次數重疊的 PSD、平均 PSD 與理論曲線比較
clear;
clc;
f = 3e6;       
T = 1e-6;    
T2 = 4e-6; 
fs = 1e8;   
num = 200;    
t = 0 : 1/fs : T-1/fs;     
t2 = 0 : 1/fs : T2-1/fs;  
t_total = 0 : 1/fs : (num*T)-1/fs;  
t_total2 = 0 : 1/fs : (num*T2)-1/fs;  
a = 1e-3;      
times = 200; 

%case1
psd_all = zeros(times, length(t_total)); 
QAMw = zeros(1, length(t_total)); 
for j = 1:times
    m = (sign(randn(1, num)) + 2*sign(randn(1, num))) + 1i*((sign(randn(1, num)) + 2*sign(randn(1, num))));
    for i = 1:num
        QAMw(1,(i-1)*length(t)+1 : i*length(t))...
        = real(m(i)).*a.*sqrt(2/T).*cos(2*pi*f*t)+imag(m(i)).*a.*sqrt(2/T).*sin(2*pi*f*t);
    end
    y = fft(QAMw);
    % f_fft = (-length(QAMw)/2:length(QAMw)/2-1) * (fs/length(QAMw));
    psd = abs(fftshift(y)).^2 / (fs * length(QAMw));
    psd_all(j,:) = psd;
end
f_fft = linspace(-fs/2, fs/2, length(psd));
psd_av = sum(psd_all) / times;


%case2
psd_all2 = zeros(times, length(t_total2)); 
QAMw2 = zeros(1, length(t_total2)); 
for j = 1:times
    m = (sign(randn(1, num)) + 2*sign(randn(1, num))) + 1i*(1*(sign(randn(1, num)) + 2*sign(randn(1, num))));
    for i = 1:num
        QAMw2(1,(i-1)*length(t2)+1 : i*length(t2))...
        = real(m(i)).*a.*sqrt(2/T2).*cos(2*pi*f*t2)+imag(m(i)).*a.*sqrt(2/T2).*sin(2*pi*f*t2);
    end
    y2 = fft(QAMw2);
    % f_fft2 = (-length(QAMw2)/2:length(QAMw2)/2-1) * (fs/length(QAMw2));
    psd2 = abs(fftshift(y2)).^2 / (fs * length(QAMw2));
    psd_all2(j,:) = psd2;
end
f_fft2 = linspace(-fs/2, fs/2, length(psd2));
psd_av2 = sum(psd_all2) / times;

figure;
subplot(211);
plot(f_fft, psd_av);
title('T=1e-6 psd');
xlim([-2*1e7 2*1e7]);
grid on;
subplot(212);
plot(f_fft, psd_all);
title('200 times overlapped psd');
xlim([-2*1e7 2*1e7]);
grid on;

figure;
subplot(211);
plot(f_fft2, psd_av2);
title('T=4e-6 psd');
xlim([-2*1e7 2*1e7]);
grid on;
subplot(212);
plot(f_fft2, psd_all2);
title('200 times overlapped psd');
xlim([-2*1e7 2*1e7]);
grid on;

%ideal--------------------------------------------- 
fc = 3e6;
a = 1e-3;

P11 = sqrt(T/2).*( sinc(T.*(f_fft-fc)) + sinc(T.*(f_fft+fc)) );
P12 = 1i.*sqrt(T/2).*( sinc(T.*(f_fft-fc)) - sinc(T.*(f_fft+fc)) );
P21 = sqrt(T2/2).*( sinc(T2.*(f_fft2-fc)) + sinc(T2.*(f_fft2+fc)) );
P22 = 1i.*sqrt(T2/2).*( sinc(T2.*(f_fft2-fc)) - sinc(T2.*(f_fft2+fc)) );

Esi = 10*a.^2;
S1 = (Esi/(2*T)).*(abs(P11).^2 + abs(P12).^2);
S2 = (Esi/(2*T2)).*(abs(P21).^2 + abs(P22).^2);
%--------------------------------------------------
figure;
subplot(211);
plot(f_fft, psd_av,f_fft, S1,'k');
xlim([-2*1e7 2*1e7]);
title('T=1e-6 psd');
legend('psd_a_v','ideal')
grid on;
subplot(212);
plot(f_fft2, psd_av2,f_fft2, S2,'k');
xlim([-2*1e7 2*1e7]);
title('T=4e-6 psd');
legend('psd_a_v','ideal')
grid on;





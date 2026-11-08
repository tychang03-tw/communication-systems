% 16-QAM 功率頻譜密度（自相關法）
% 每組波形先算時間自相關 R(tau) 再做 FFT，100 次平均；比較 T = 1 us 與 4 us 和理論 PSD
% 輸出：平均 PSD 與理論曲線、所有次數重疊的 PSD
clear;
clc;

fc = 3e6;
T = 1e-6;
T2 = 4e-6;
fs = 1e8;
num = 200;
t = 0 : 1/fs : T-1/fs;
t_total = 0 : 1/fs : (num*T)-1/fs;
t2 = 0 : 1/fs : T2-1/fs;
t_total2 = 0 : 1/fs : (num*T2)-1/fs;
a = 1e-3;
times = 100;
QAMw = zeros(1, length(t_total));
QAMw2 = zeros(1, length(t_total2));
%case 1
for j = 1:times
    m = (sign(randn(1, num)) + 2*sign(randn(1, num))) + 1i*((sign(randn(1, num)) + 2*sign(randn(1, num))));
    for i = 1:num
        QAMw(1,(i-1)*length(t)+1 : i*length(t))...
            = real(m(i)).*a.*sqrt(2/T).*cos(2*pi*fc*t)+imag(m(i)).*a.*sqrt(2/T).*sin(2*pi*fc*t);
    end
    R_tau(j,:) = conv(QAMw,fliplr(QAMw))./ length(QAMw);
    S_f(j,:) = abs(fftshift(fft(R_tau(j,:))))./fs;
    % R_tau(j,:) = xcorr(QAMw, 'biased');
    % S_f(j,:) = abs(fftshift(fft(R_tau(j,:)))) ./ fs; 
end
psd_av = sum(S_f) ./ times;
f_fft = linspace(-fs/2, fs/2, length(S_f));
%case 2
for j = 1:times
    m = (sign(randn(1, num)) + 2*sign(randn(1, num))) + 1i*((sign(randn(1, num)) + 2*sign(randn(1, num))));
    for i = 1:num
        QAMw2(1,(i-1)*length(t2)+1 : i*length(t2))...
            = real(m(i)).*a.*sqrt(2/T2).*cos(2*pi*fc*t2)+imag(m(i)).*a.*sqrt(2/T2).*sin(2*pi*fc*t2);
    end
    R_tau2(j,:) = conv(QAMw2,fliplr(QAMw2))./ length(QAMw2);
    S_f2(j,:) = abs(fftshift(fft(R_tau2(j,:))))./fs;
end
psd_av2 = sum(S_f2) ./ times;
f_fft2 = linspace(-fs/2, fs/2, length(S_f2));
%ideal--------------------------------------------- 
P11 = sqrt(T/2).*( sinc(T.*(f_fft-fc)) + sinc(T.*(f_fft+fc)) );
P12 = 1i.*sqrt(T/2).*( sinc(T.*(f_fft-fc)) - sinc(T.*(f_fft+fc)) );
P21 = sqrt(T2/2).*( sinc(T2.*(f_fft2-fc)) + sinc(T2.*(f_fft2+fc)) );
P22 = 1i.*sqrt(T2/2).*( sinc(T2.*(f_fft2-fc)) - sinc(T2.*(f_fft2+fc)) );
Esi = 10*a.^2;
S1 = (Esi/(2*T)).*(abs(P11).^2 + abs(P12).^2);
S2 = (Esi/(2*T2)).*(abs(P21).^2 + abs(P22).^2);
%--------------------------------------------------
figure;
subplot(411)
plot(f_fft,psd_av,f_fft,S1,'k')
title('T=1e-6 psd');
legend('psd_a_v','ideal')
xlim([-1.5*1e7 1.5*1e7]);
grid on
subplot(412)
plot(f_fft,S_f)
title('200 times overlapped psd');
xlim([-1.5*1e7 1.5*1e7]);
grid on
subplot(413)
plot(f_fft2,psd_av2,f_fft2,S2,'k')
title('T=4e-6 psd');
legend('psd_a_v','ideal')
xlim([-1.5*1e7 1.5*1e7]);
grid on
subplot(414)
plot(f_fft2,S_f2)
title('200 times overlapped psd');
xlim([-1.5*1e7 1.5*1e7]);
grid on
sgtitle('16-QAM psd')





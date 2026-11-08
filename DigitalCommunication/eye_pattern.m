% 眼圖
% NRZ 訊號（T = 1 ms）通過頻寬 W 的理想低通通道 h(t) = W*sinc(W*t)，W = 700 / 900 / 1100 / 150k Hz，各重疊 100 組隨機序列
% 輸出：四種頻寬下的眼圖
clear;
clc;

num = 10;
T = 1e-3;
fs = 1e6;
t = 0:1/fs:T-1/fs;
t_total = -5*T:1/fs:5*T-1/fs;
t_conv = -10*T:1/fs:10*T -2/fs;
W_values = [700, 900, 1100, 150e3];
times = 100;

figure;
for k = 1:length(W_values)
    W = W_values(k);
    h_c = W * sinc(W * t_total);
    subplot(2,2,k);
    hold on;
    for j = 1:times
        m = sign(randn(1, num));
        Tx = zeros(1, length(t_total));
        for i = 1:num
            Tx((i-1)*length(t) + 1 : i*length(t)) = m(i);
        end
        Rx(j,:) = conv(Tx,h_c)./fs;
        plot(t_conv, Rx(j,:));
    end
    title(['W = ', num2str(W), ' Hz']);
    xlim([-2*T  2*T])
    grid on;
    hold off;
end

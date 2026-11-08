% QPSK 的 MAP 決策區域
% 依 D_i(r) = -||r - s_i||^2 + N0*ln P(m_i) 計算平面上每一點的判決；程式含等機率（delta）與不等機率（delta2）兩組，目前畫不等機率
% 輸出：決策區域圖
clear;
clc;
f = 3e6;
T = 1e-6;
fs = 1e8;
num = 1e3;
t = 0 : 1/fs : T-1/fs;
t_total = 0 : 1/fs : (num*T)-1/fs;
Es = 1;
EsN0 = 3;
EsN0w = 10.^(EsN0./10);
N0w = Es./EsN0w;
QAMw = zeros(1, length(t_total));

M = 4;
m = randi([0 M-1], 1, num); % 產生0~M-1數字
c = exp(1i * (2 * pi * m / M)); % 信號用極座標表示
for i = 1:num
    QAMw(1,(i-1)*length(t)+1 : i*length(t)) ...
        = real(c(i)).*sqrt(Es).*sqrt(2/T).*cos(2*pi*f*t) ...
        + imag(c(i)).*sqrt(Es).*sqrt(2/T).*sin(2*pi*f*t);
end
E = (QAMw*QAMw')/num;
noise = sqrt((0.5.*10^(-EsN0 ./ 10))*E);
n = noise .* randn(size(QAMw));
recieve_QAMw = QAMw + n;

I = zeros(num, length(t));
Q = zeros(num, length(t));
for i = 1:num
    r = recieve_QAMw(1, (i-1)*length(t)+1 : i*length(t));
    I(i, :) = sum(r .* sqrt(2/T) .* cos(2*pi*f*t)) .* (1/fs);
    Q(i, :) = sum(r .* sqrt(2/T) .* sin(2*pi*f*t)) .* (1/fs);
end
y = I(:, 1) + 1i*Q(:, 1);
yx = [I(:, 1) Q(:, 1)];

%case1
p0 = 1/4;
p1 = 1/4;
p2 = 1/4;
p3 = 1/4;

p11 = [p0, p1, p2, p3];
delta = N0w .* log(p11);

%case2
p0_1 = 4/24;
p1_1 = 13/24;
p2_1 = 1/24;
p3_1 = 6/24;

p22 = [p0_1, p1_1, p2_1, p3_1];
delta2 = N0w .* log(p22);

s = [1 ,0; 0 ,1; -1 ,0; 0 ,-1];% 參考點
x_grid = linspace(-10, 10, 500); % 從 '哪裡' 到 '哪裡' 生成 '多少網格點'
y_grid = linspace(-10, 10, 500);
[X,Y] = meshgrid(x_grid, y_grid); % 將x_grid跟y_grid兩個一維向量變成二維
I = zeros(size(X));

% 算all網格點決策值並找到最近的參考點
for i = 1:length(x_grid)
    for j = 1:length(y_grid)
        point = [X(i,j), Y(i,j)];
        D = zeros(1, 4);
        for k = 1:4
            D(1,k) = -norm(point - s(k,:)).^2 + delta2(k); % Di(r)=-||r-si||^2+N0*ln(P(mi)
        end
        [~, I(i,j)] = max(D);% Ik={r|Dk(r)>=Di(r),i~=k}，k=0~M-1
    end
end


figure;
imagesc(x_grid, y_grid, I);
% colormap(jet(4)); % 使用 4 種顏色來表示不同的決策區域
set(gca,'YDir','normal'); % imagesc語法，X軸為左到右增加，Y軸為下到上增加(reverse)，所以最後結果將Y軸轉回來(normal)

hold on;
%scatter(real(y), imag(y), 'p');
scatter(real(c), imag(c), 'p');
xlim([-3 3]);
ylim([-3 3]);
title('Decision region');
grid on;
hold off;

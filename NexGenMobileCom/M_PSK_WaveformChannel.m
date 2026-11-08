% M-PSK 錯誤率（波形通道）
% BPSK / QPSK / 8PSK / 16PSK 的 passband 波形加 AWGN，手刻 correlator 取得 I/Q 後做最小距離判決（10^5 個符號，執行時間較長）
% 輸出：BER（以 SER / log2(M) 近似）vs. Es/N0
clear;
f = 3e6;
T = 1e-6;
fs = 1e8;
num = 100000;
t = 0 : 1/fs : T-1/fs;
t_total = 0 : 1/fs : (num*T)-1/fs;
Es = 1e-3;
EsN0 = 0:2:15;
PSKw = zeros(1, length(t_total));
M = [2 4 8 16];

for s = 1:4
    MM = M(s);
    m = randi([0 MM-1], 1, num);
    c = exp(1i * (2 * pi * m / MM));
    cc = exp(1i * (2 * pi * (0:MM-1) / MM));
    cx = zeros(MM, 2);
    b = log2(MM);

    for j = 1:MM
        cx(j, :) = [real(cc(j)) imag(cc(j))];
    end

    for i = 1:num
        PSKw(1,(i-1)*length(t)+1 : i*length(t))...
            = real(c(i)).*sqrt(Es).*sqrt(2/T).*cos(2*pi*f*t)+imag(c(i)).*sqrt(Es).*sqrt(2/T).*sin(2*pi*f*t);
    end

    E = (PSKw*PSKw')./num;

    for k=1:length(EsN0)

        noise = sqrt((0.5.*10^(-EsN0(k) ./ 10))*E);
        n = noise .* randn(size(PSKw));
        recieve_PSKw = PSKw + n ;

        I = zeros(num,length(t));
        Q = zeros(num,length(t));
        
        for i = 1:num
            r = recieve_PSKw(1,(i-1)*length(t)+1 : i*length(t));
            I(i,:) = sum(r .* sqrt(2/T) .* cos(2*pi*f*t)).*(1/fs);
            Q(i,:) = sum(r .* sqrt(2/T) .* sin(2*pi*f*t)).*(1/fs);
        end
        y = I(:,1) + 1i*Q(:,1);
        yx = [I(:,1) Q(:,1)];

        d = zeros(num, 1);
        for i = 1:num
            distances = sqrt(sum((yx(i,:) - cx).^2,2));% 語法中 yx(i,:)-cx 雖然是 1行 減 16 行
            [~, d(i)] = min(distances);                     % 但會產生 1行 減掉 16行的每一行 的結果
        end                                                 % (yx(i,:) - cx).^2 每個元素平方
        detected = (d-1)';                                  % sum((yx(i,:) - cx).^2,2) 兩個維度的值相加 
        numErrors(k) = sum(detected ~= m);                  % 等效於平方和
    end
    SER = numErrors / num;
    BER(s,:) = SER./ b;
end

figure;
semilogy(EsN0,BER(1,:),'k-o',EsN0,BER(2,:),'g-s',EsN0,BER(3,:),'b-^',EsN0,BER(4,:),'m-d');
xlabel('Es/N0(dB)');
ylabel('BER');
title('M-PSK compare')
legend('B-PSK','Q-PSK','8-PSK','16-PSK')
grid;


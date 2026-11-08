% 室內射線追蹤（ray tracing）
% 50 m x 50 m 房間，Tx 與 Rx 距離 2～20 m，用 image method 計算 LOS 加上 1～3 次牆面反射（反射係數 -0.7）的接收功率
% 輸出：接收功率 vs. 距離，並以線性擬合估計路徑損耗斜率
clc;
clear;
c = 3e8;
f = 1e9;
lambda = c/f;
Pt = 1;
Gt = 1.6;
Gr = 1.6;
a = -0.7;
P0 = Pt*Gt*Gr*(lambda/(4*pi))^2;
room = [0,50,0,50]; % [x_min x_max y_min y_max]
rx = [25, 25];
d = 2:0.1:20;
Pr_all = zeros(length(d), 4);
for i = 1:length(d)
    tx = [25 + d(i), 25];
    % LOS--------------------------------------------------
    d_los = norm(tx - rx);
    Pr_los = 1/d_los * exp(-1i*2*pi*d_los/lambda);
    Pr_all(i,1) = P0*abs(Pr_los)^2;
    % 反射一次--------------------------------------------------
    Pr_ref1 = 0;
    walls = 1:4; % 左右上下牆
    for w1 = walls
        switch w1
            case 1, tx1 = [-tx(1), tx(2)];
            case 2, tx1 = [2*room(2)-tx(1), tx(2)];
            case 3, tx1 = [tx(1), -tx(2)];
            case 4, tx1 = [tx(1), 2*room(4)-tx(2)];
        end
        d_ref1 = norm(tx1 - rx);
        Pr_ref1 = Pr_ref1 + a/d_ref1 * exp(-1i*2*pi*d_ref1/lambda);
    end
    Pr_all(i,2) = P0*abs(Pr_los + Pr_ref1)^2;
    % 反射兩次--------------------------------------------------
    Pr_ref2 = 0;
    walls = 1:4;
    for w1 = walls
        for w2 = walls
            if w2 == w1
                continue; % 跳過同一面牆兩次反射
            end
            switch w1
                case 1, tx1 = [-tx(1), tx(2)];
                case 2, tx1 = [2*room(2)-tx(1), tx(2)];
                case 3, tx1 = [tx(1), -tx(2)];
                case 4, tx1 = [tx(1), 2*room(4)-tx(2)];
            end
            switch w2
                case 1, tx2 = [-tx1(1), tx1(2)];
                case 2, tx2 = [2*room(2)-tx1(1), tx1(2)];
                case 3, tx2 = [tx1(1), -tx1(2)];
                case 4, tx2 = [tx1(1), 2*room(4)-tx1(2)];
            end
            d_ref2 = norm(tx2 - rx);
            Pr_ref2 = Pr_ref2 + a^2/d_ref2 * exp(-1i*2*pi*d_ref2/lambda);
        end
    end
    Pr_all(i,3) = P0*abs(Pr_los + Pr_ref1 + Pr_ref2)^2;
    % 反射三次--------------------------------------------------
    Pr_ref3 = 0;
    walls = 1:4;
    for w1 = walls
        for w2 = walls
            for w3 = walls
                % 避免同一牆連續三次
                if w2 == w1 || w3 == w2
                    continue;
                end
                switch w1
                    case 1, tx1 = [-tx(1), tx(2)];
                    case 2, tx1 = [2*room(2)-tx(1), tx(2)];
                    case 3, tx1 = [tx(1), -tx(2)];
                    case 4, tx1 = [tx(1), 2*room(4)-tx(2)];
                end
                switch w2
                    case 1, tx2 = [-tx1(1), tx1(2)];
                    case 2, tx2 = [2*room(2)-tx1(1), tx1(2)];
                    case 3, tx2 = [tx1(1), -tx1(2)];
                    case 4, tx2 = [tx1(1), 2*room(4)-tx1(2)];
                end
                switch w3
                    case 1, tx3 = [-tx2(1), tx2(2)];
                    case 2, tx3 = [2*room(2)-tx2(1), tx2(2)];
                    case 3, tx3 = [tx2(1), -tx2(2)];
                    case 4, tx3 = [tx2(1), 2*room(4)-tx2(2)];
                end
                d_ref3 = norm(tx3 - rx);
                Pr_ref3 = Pr_ref3 + a^3/d_ref3 * exp(-1i*2*pi*d_ref3/lambda);
            end
        end
    end
    Pr_all(i,4) = P0*abs(Pr_los + Pr_ref1 + Pr_ref2 + Pr_ref3)^2;
end

figure;
for k = 1:4
    Pr_dB = 10*log10(Pr_all(:,k));
    % 產生線性擬合斜線
    p = polyfit(10*log10(d), Pr_dB, 1);    
    P_fit = polyval(p, 10*log10(d));
    slope = p(1); % p = [斜率, 截距]
    
    subplot(2,2,k);
    h1 = semilogx(d, Pr_dB); hold on;
    h2 = semilogx(d, P_fit);
    legend(h2,sprintf('Slope=%.2f',slope),'Location','northwest');
    
    grid on;
    xlabel('Distance (m)');
    ylabel('Relative Power (dB)');
    set(gca,'XTick',[2 5 10 20]);
    ylim([-70 -25]);
   
    switch k
        case 1, title_txt = 'LOS';
        case 2, title_txt = 'LOS + FirstOrderReflection';
        case 3, title_txt = 'LOS + FirstOrderReflection + SecondOrderReflection';
        case 4, title_txt = 'LOS + FirstOrderReflection + SecondOrderReflection + ThirdOrderReflection';
    end
    title(title_txt);
end

% IIR 濾波器分析
% 以 filterDesigner 設計出的 4 階 IIR 係數，計算頻率響應與 group delay，並印出轉移函數 H(z)
% 輸出：振幅響應、group delay、H(z)
format long;
% 分母
Den =[1.000000000000000  -0.649546903894243   0.859609692401431  -0.195456561685657   0.071282237383512];
% 分子
Num =[0.225681883788483  -0.685207604739853   0.950644377680610  -0.685207604739853   0.225681883788483];
[H,W] = freqz(Num,Den,1000) ;
[g,w] = grpdelay(Num, Den, 1000);
figure;
subplot(211)
plot(W/pi,abs(H));
ylim([0 1.1])
xlabel('W/pi');
ylabel('Magnitude response');
grid on
subplot(212)
plot(W/pi,abs(g));
xlabel('W/pi');
ylabel('group delay');
grid on
sgtitle('IIR filter')

N_num = length(Num); 
N_den = length(Den); 
Num_str = "";
Den_str = ""; 

for i = 1:N_num
    coeff = Num(i);
    power = i - 1; % z^{-power}
    if coeff > 0 && i > 1
        Num_str = strcat(Num_str, " + "); 
    elseif coeff < 0
        Num_str = strcat(Num_str, " - "); 
    end
    Num_str = strcat(Num_str, sprintf("%.6f", abs(coeff))); % 添加係數
    if power > 0
        Num_str = strcat(Num_str, sprintf(" z^-%d", power)); % 添加z次數
    end
end
for i = 1:N_den
    coeff = Den(i);
    power = i - 1; % z^{-power}
    if coeff > 0 && i > 1
        Den_str = strcat(Den_str, " + ");
    elseif coeff < 0
        Den_str = strcat(Den_str, " - ");
    end
    Den_str = strcat(Den_str, sprintf("%.6f", abs(coeff)));
    if power > 0
        Den_str = strcat(Den_str, sprintf(" z^-%d", power));
    end
end
disp("H(z) = (" + Num_str + ") / (" + Den_str + ")");

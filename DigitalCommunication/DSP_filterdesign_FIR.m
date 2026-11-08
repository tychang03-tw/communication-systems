% FIR 濾波器分析
% 以 filterDesigner 設計出的 13-tap FIR 係數，計算頻率響應與 group delay，並印出轉移函數 H(z)
% 輸出：振幅響應、group delay、H(z)
format long;
Num = ...
   [-0.034469549768619  -0.014546035001740   0.059595885386943   0.059605406772247  -0.072051470491630...
   -0.306542463302543   0.589264662951166  -0.306542463302543  -0.072051470491630   0.059605406772247...
    0.059595885386943  -0.014546035001740  -0.034469549768619];
[H,W] = freqz(Num,1000) ;
[g,w] = grpdelay(Num,1000);
figure;
subplot(211)
plot(W/pi,abs(H));
xlabel('W/pi');
ylabel('Magnitude response');
grid on
subplot(212)
plot(W/pi,abs(g));
xlabel('W/pi');
ylabel('group delay');
grid on
sgtitle('FIR filter')

N = length(Num); 
Hz_str = "";
for i = 1:N
    coeff = Num(i);
    power = i - 1; % z^{-power}
    if coeff > 0 && i > 1
        Hz_str = strcat(Hz_str, " + "); % 正號
    elseif coeff < 0
        Hz_str = strcat(Hz_str, " - "); % 負號
    end
    Hz_str = strcat(Hz_str, sprintf("%.6f", abs(coeff))); % 添加係數
    if power > 0
        Hz_str = strcat(Hz_str, sprintf(" z^-%d", power)); % 添加次數
    end
end
disp("H(z) = " + Hz_str);

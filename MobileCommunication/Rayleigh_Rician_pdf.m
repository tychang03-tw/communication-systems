% Rayleigh / Rician 分布
% 由二路高斯隨機變數產生 envelope，比較模擬 histogram 與理論 PDF
% 輸出：Rayleigh 理論 / 模擬、Rician 模擬，以及 Rician 理論與模擬比較
clear
clc
% Rayleigh theory --------------------------------
z = 0:0.01:10-0.01; % A + sigma*randn(1,num)
f1 = raylpdf(z,1);
f2 = raylpdf(z,2); 
figure;
subplot(221);
plot(z,f1,z,f2);
axis([0,10,0,0.8]);
title('Rayleigh theory');
legend('mean=0,Var=1','mean=0,Var=4');
xlabel('z');
ylabel('f(z)');
grid;
% Rayleigh simulation ---------------------------------
num = 1e6;
x1 = randn(1,num);        
y1 = randn(1,num);
z1 = sqrt(x1.^2 + y1.^2);
[Y1,~] = hist(z1,z);
x2 = 2*randn(1,num);   
y2 = 2*randn(1,num);
z2 = sqrt(x2.^2 + y2.^2);
[Y2,~] = hist(z2,z);
subplot(222);
bar(z,Y1/(num*0.01));
hold on;
bar(z,Y2/(num*0.01));
axis([0,10,0,0.8]);
title('Rayleigh simulation');
legend('mean=0,Var=1','mean=0,Var=4');
xlabel('z');
ylabel('f(z)');
grid;
% Rician theory ---------------------------------
sigma = 1;
A = 1; 
f_Rician = (z/sigma^2).*exp(-(z.^2+A^2)/(2*sigma^2)).*besseli(0,A*z./(sigma^2)); % order 0
% Rician simulation -----------------------------------
a = 0.6; % LOS權重
x3 = randn(1,num); % NLOS 
y3 = randn(1,num);
z3 = sqrt(x3.^2 + y3.^2);
[Y3,~] = hist(z3,z);
x4 = 1 + randn(1,num); % LOS
% y4 = randn(1,num);
% z4 = sqrt(x4.^2 + y4.^2);
[Y4,~] = hist(x4,z);
r1 = (1-a)*Y3 + a*Y4;
subplot(223);
bar(z,r1/(num*0.01),'b');
axis([0,10,0,0.8]);
title('Rician simulation');
legend('alpha = 0.6')
xlabel('z');
ylabel('f(z)');
grid;
% compare------------------------------------
subplot(224);
plot(z,f_Rician,'m');
hold on;
bar(z,r1/(num*0.01),'b');
axis([0,10,0,0.8]);
title('Compare');
legend('Rician theory','Rician simulation (a=0.6)');
xlabel('z');
ylabel('f(z)');
grid;



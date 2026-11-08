% 積分傾卸接收機（integrate-and-dump receiver）
% 用 10^5 個 bit 模擬 antipodal 訊號，比較振幅 20 / 30 mV 與 T_b = 1 / 0.1 ms 下的 BER 與理論值
% 輸出：BER vs. 雜訊功率頻譜密度 N0
clear
clc
num=10.^5;
A1 = 0.02;
A2 = 0.03;
T1 = 1/1000;
T2 = 0.1/1000;
N0 = -20:-2:-86;
N0w=10.^(N0./10);

for i = 1:length(N0)
    x1(i) = sqrt(2.*(A1.^2).*T1./N0w(i));
    x2(i) = sqrt(2.*(A2.^2).*T1./N0w(i));
    x3(i) = sqrt(2.*(A1.^2).*T2./N0w(i));
    Q1(i) = (1/sqrt(2*pi)) * integral(@(t) exp(-t.^2 / 2), x1(i), 60);% A1 T1
    Q2(i) = (1/sqrt(2*pi)) * integral(@(t) exp(-t.^2 / 2), x2(i), 60);% A2 T1
    Q3(i) = (1/sqrt(2*pi)) * integral(@(t) exp(-t.^2 / 2), x3(i), 60);% A1 T2
end
semilogy(N0, Q1,'r-o');
xlim([-85 -45])
ylim([10.^-7 1])
xlabel('Noise Power Spectral(dB)')
ylabel('BER')
title('ideal')
grid on
set(gca,'XDir','reverse')

m=round(rand(1, num)); 
s1=(2*m-ones(size(m)))*A1; 
s2=(2*m-ones(size(m)))*A2; 

for ii = 1:length(N0)
    nois1=10^(N0(ii)/20)*sqrt(T1/2);
    nois2=10^(N0(ii)/20)*sqrt(T2/2);
    n1=nois1*randn(size(s1)); 
    n2=n1./nois1.*nois2; 

    v1 = s1.*T1 + n1;
    v2 = s2.*T1 + n1;
    v3 = s1.*T2 + n2;

    err1(ii,:) = sign(s1)~=sign(v1); % A1 T1
    err2(ii,:) = sign(s2)~=sign(v2); % A2 T1
    err3(ii,:) = sign(s1)~=sign(v3); % A1 T2

end
error1 = sum(err1.')/num;
error2 = sum(err2.')/num;
error3 = sum(err3.')/num;

%第一題 T1 A1 A2
figure;
title('第一題');
xlabel('Noise Power Spectral(dB)');
ylabel('BER');
semilogy(N0,error1,'r-o',N0,error2,'g-o',N0,Q1,'b-o',N0,Q2,'k-o')
xlim([-85 -45])
ylim([10.^-7 1])
grid on
legend('T1 = 10^-^3 sec, A1 = 20mV', 'T1 = 10^-^3 sec, A2 = 30mV','true P_e_1','true P_e_2','location','southwest');
set(gca,'XDir','reverse')

%第二題 A1 T1 T2
figure
title(['第二題']);
xlabel('Noise Power Spectral(dB)');
ylabel('BER');
semilogy(N0,error1,'r-o',N0,error3,'g-o',N0,Q1,'b-o',N0,Q3,'k-o')
xlim([-85 -45])
ylim([10.^-7 1])
grid on
legend('T1 = 10^-^3 sec, A1 = 20mV', 'T1 = 10^-^4 sec, A2 = 20mV','true P_e_1','true P_e_3','location','southwest');
set(gca,'XDir','reverse')



% 衰落通道下的 BPSK 錯誤率
% Rayleigh fading 用三種方式產生（raylrnd、兩路高斯、circular symmetric complex Gaussian），另以 LOS 權重 a 模擬 Rician 通道
% 輸出：AWGN 理論值與各通道的 BER vs. Eb/N0
clear
clc
bitnum = 1e5;
s = randi([0 1],1,bitnum);
x = pskmod(s,2,pi/2);

% Rayleigh 語法
h = raylrnd(1,1,bitnum); % 振幅為Rayleigh
r = h.*x;

% 二路高斯
x1 = pskmod(s,2,0);
x2 = pskmod(s,2,pi/2);
nor1 = abs(normrnd(0,1,1,bitnum)).*x1;
nor2 = abs(normrnd(0,1,1,bitnum)).*x2;

% Circular Symmetric Guassian 寫法
h2 = abs(normrnd(0,1,1,bitnum)+1i*normrnd(0,1,1,bitnum)); % abs : 振幅為Rayleigh
r2 = h2.*x;

for Eb_N0 = 0:20
    BER_BpskTheory(1,Eb_N0+1) = qfunc(sqrt(2*10.^(Eb_N0/10)));

    % Rayleigh 語法 --------------------------------
    y = awgn(r,Eb_N0,'measured'); % 預設為 snr = Es/N0
    s_hat = pskdemod(y,2,pi/2);
    [num,ratio] = biterr(s,s_hat);
    BER_BpskRayleigh(1,Eb_N0+1) = ratio;

    % 二路高斯 ---------------------------------
    y2 = awgn(nor1+nor2,Eb_N0,'measured');
    s_hat2 = pskdemod(y2,2,(0 + pi/2)/2);
    [num2,ratio2] = biterr(s,s_hat2);
    BER_BpskRayleigh2(1,Eb_N0+1) = ratio2;

    % Circular Symmetric Guassian 寫法 --------------------------------
    y3 = awgn(r2,Eb_N0,'measured');
    s_hat3 = pskdemod(y3,2,pi/2);
    [num3,ratio3] = biterr(s,s_hat3);
    BER_BpskRayleigh3(1,Eb_N0+1) = ratio3;

    % Rician LOS權重為a ---------------------------------
    a = 0.5;
    y4 = awgn(a*r+(1-a)*x,Eb_N0,'measured');
    s_hat4 = pskdemod(y4,2,pi/2);
    [num4,ratio4] = biterr(s,s_hat4);
    BER_BpskRayleigh4(1,Eb_N0+1) = ratio4;

end
semilogy(0:20,BER_BpskTheory,'m-.', ...
    0:20,BER_BpskRayleigh,'k-p', ...
    0:20,BER_BpskRayleigh2,'b-^', ...
    0:20,BER_BpskRayleigh3,'r-o', ...
    0:20,BER_BpskRayleigh4,'g--')
xlabel('E_b/N_0 (dB)'); ylabel('BER');
legend('AWGN','Rayleigh','two path Gaussian','Circular Symmetric Gaussian','Rician');
ylim([1e-7 1])
grid


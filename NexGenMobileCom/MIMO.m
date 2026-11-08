% MIMO 分集
% BPSK 比較 SISO、MISO（4 根 Tx，Tx beamforming）、SIMO（4 根 Rx，MRC）與 2x2 Alamouti code，對 10^5 次通道實現取平均
% 輸出：BER vs. SNR
clear;
clc;
data_num = 10;
binary_source = 2*randi([0 1],1,data_num) - 1;
times = 100000;
M = 4;

for j = 1 : times
    clc;
    disp(j);
    % SISO tx --------------------------------
    h = randn(1,1);
    r_siso = h.*binary_source;

    % MISO tx --------------------------------
    h2 = randn(M,1);
    W_Tx = (conj(h2) ./ norm(h2)); % tx beamformer
    r_miso = binary_source.*W_Tx.*h2;

    % SIMO tx ---------------------------------
    h3 = randn(M,1);
    r_simo = binary_source.*h3;
    W_Rx = h3; % rx beamformer (MRC)

    % Alamouti code tx -------------------------
    h4 = randn(2,2);
    s = reshape(sqrt(2)^-1.*binary_source,2,[]);
    s_Alamouti = [];
    for i = 1:(data_num/2)
        ss = s(:,i);
        sss = [ss(1,1) -ss(2,1)';
            ss(2,1) ss(1,1)'];
        s_Alamouti = [s_Alamouti sss];
    end
    r_Alamouti(1,:) = h4(1,1).*s_Alamouti(1,:) + h4(1,2).*s_Alamouti(2,:);
    r_Alamouti(2,:) = h4(2,1).*s_Alamouti(1,:) + h4(2,2).*s_Alamouti(2,:);
    % ------------------------------------------
    SNR = 0:3:30;
    for i = 1:length(SNR)
        N0 = M ./ (10.^(SNR / 10));
        N02 = 4 ./ (10.^(SNR / 10));
        % SISO rx ------------------------------
        noise1 = sqrt(N0(i)/2) * randn(size(r_siso));
        y_siso = r_siso + noise1;
        s_hat1 = sign(conj(h).*y_siso);
        numErrors = sum(s_hat1 ~= binary_source);
        BER(j,i) = numErrors / data_num;
        % MISO rx ------------------------------
        r_miso = sum(r_miso,1);
        noise2 = sqrt(N0(i)/2) * (randn(size(r_miso)));
        y_miso = r_miso + noise2;
        s_hat2 = sign(y_miso);
        numErrors2 = sum(s_hat2 ~= binary_source);
        BER2(j,i) = numErrors2 / data_num;
        % SIMO rx -------------------------------
        for k = 1:M
            noise3(k,:) = sqrt(N0(i)/2) * randn(size(r_simo(k,:)));
            y_simo(k,:) = r_simo(k,:) + noise3(k,:);
        end
        s_hat3 = sign(W_Rx'*y_simo);
        numErrors3 = sum(s_hat3 ~= binary_source);
        BER3(j,i) = numErrors3 / data_num;
        % Alamouti code rx ----------------------
        for k = 1:2
            noise4(k,:) = sqrt(N02(i)/2) * randn(size(r_Alamouti(k,:)));
            y_Alamouti(k,:) = r_Alamouti(k,:) + noise4(k,:);
        end
        for k = 1:(data_num/2)
            s_Rx(1,2*k-1) = h4(:,1)'*y_Alamouti(:,2*k-1) + y_Alamouti(:,2*k)'*h4(:,2);
            s_Rx(1,2*k) = h4(:,2)'*y_Alamouti(:,2*k-1) - y_Alamouti(:,2*k)'*h4(:,1);
        end
        s_hat4 = sign(s_Rx);
        numErrors4 = sum(s_hat4 ~= binary_source);
        BER4(j,i) = numErrors4 / data_num;
    end
end
BER_av = sum(BER)./times;
BER_av2 = sum(BER2)./times;
BER_av3 = sum(BER3)./times;
BER_av4 = sum(BER4)./times;

semilogy(SNR,BER_av,'k--x',SNR,BER_av2,'r:s',SNR,BER_av3,'k-.^',SNR,BER_av4,'b-o');
xlabel('SNR');
ylabel('BER');
xlim([0 25])
ylim([1e-6 1]);
legend('SISO', 'MISO', 'SIMO','Alamouti code','Location', 'southwest');
grid on;

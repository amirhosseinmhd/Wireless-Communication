
clc;clear;close all;

%% Theoritical 
SNR_db=-10:1:10;
SNR_mag = 10.^(SNR_db./10);
figure
for L = 1 : 5
    Pe = zeros(1, length(SNR_mag));
    for i = 1 : length(SNR_mag)
        mu = sqrt(SNR_mag(i)/(1+SNR_mag(i)));
        temp = 0;
        for k = 0 : L - 1
            temp = temp + nchoosek(L -1 +k , k) * ((1+mu)/2)^k;
        end
        Pe(i) = ((1-mu)/2).^L .*temp;
    end
    semilogy(SNR_db ,Pe,'LineWidth',1.5);
    hold on;
end
xlabel("SNR");
ylabel("BER");
legend('L = 1','L = 2','L = 3','L = 4','L = 5')
title("BPSK With Time Diversitiy - Theoritical")
grid on
hold off

%% Simulation

M = 2;
N= 10^5; % For smoother Figure Please Change Here
data_ =randi([0,1],1,N);
real_data = real(pskmod(data_,M));
SNR_dB = -10:0.2:10; % And Also Step here
SNR_mag = 10.^(SNR_dB/10);

len_snr = length(SNR_mag);
noise = zeros(len_snr, length(real_data));
rxSig = zeros(5, len_snr, length(real_data));
figure;
for L = 1 : 5
    h = zeros(L, N);
    for l=1:1:L
        h(l,:)= sqrt(1/2)* (wgn(1, N, 0)+ 1i *wgn(1, N, 0));
        hstar_abs_h = conj(h)./(abs(h));
        for i=1:1:len_snr
            noise(i,:)=((SNR_mag(i)./sqrt(2)).^(0.5))*(randn(1,length(real_data))...
                +1i*randn(1,length(real_data)));
            rxSig(l,i,:)=(noise(i,:)+real_data.*h(l,:)).*hstar_abs_h(l,:);
       end
    end
    rxSig_acc = zeros(1,len_snr,N);
    for i=1:1:L
        rxSig_acc = rxSig_acc+rxSig(i,:,:);
    end
    BER = zeros(len_snr, 1);
    for i=1:1:len_snr
        demodulatedSig(i,:) = pskdemod(rxSig_acc(1,i,:) , 2);
        BER(i) = 1 - sum(demodulatedSig(i,:)==data_)/N;
    end
    semilogy(flip(SNR_dB) , BER, 'LineWidth',1.5);
    hold on;
   
end
xlabel("SNR");
ylabel("BER");
legend('L = 1','L = 2','L = 3','L = 4','L = 5')
title("BPSK With Time Diversitiy")
grid on
hold off;
%% 

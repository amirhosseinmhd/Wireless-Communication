%% Wireless Communication Part 1 Question 6

clc;clear;close all;

M = 2;
N = 10^6;
data_ = randi([0,1],1,N);
real_data = real(pskmod(data_,M));
h = 1/sqrt(2) * ( (wgn(1,N,0)+1i*wgn(1,N,0)));
noise_ = sqrt(1/2) * ( randn(1,N) + 1i*randn(1,N) );
SNR_dB = -10:1:10;
SNR_mag = 10.^(SNR_dB./10);
len_snr = length(SNR_mag);

BER_Theory = (1+2*(1-(1/2 - 1/2*(1+2./SNR_mag)... 
    .^(-1/2)))).*(1/2 - 1/2*(1+2./SNR_mag).^(-1/2)).^2;

s = zeros(2, N);
s(:, 1:2:end) = (1/sqrt(2))*reshape(real_data,2,N/2);
s(:, 2:2:end) = (1/sqrt(2))*(kron(ones(1, N/2),[-1;1]).*flipud(reshape(conj(real_data), 2, N/2)));
channel = zeros(2,N);
channel(:,1:2:end) = reshape(h,2,N/2);
channel(:,2:2:end) = kron(ones(1,N/2),[1;-1]).*flipud(reshape(h,2,N/2));
channel(1,:) = conj(channel(1,:));
channel_P = sum(channel.*conj(channel),1);
H = kron(reshape(h,2,N/2), ones(1,2));
BER = zeros(1, len_snr);
for i = 1:length(SNR_dB)
    ym = sum(H.*s,1) + 10^(-SNR_dB(i)/20)*noise_;
    temp = kron(reshape(ym,2,N/2),ones(1,2));
    temp(2,:) = conj(temp(2,:));
    rxSig = sum(channel.*temp,1)./channel_P;
    rxSig(2:2:end) = conj(rxSig(2:2:end));
    final(i,:) = real(rxSig)<0;
    BER(i) = 1 - sum(final(i,:)==data_)/N;
end


semilogy(SNR_dB,BER,SNR_dB,BER_Theory);
legend('Simulation ','Theoretical');
xlabel("SNR");ylabel("BER")
title('Alamouti Scheme');
grid on



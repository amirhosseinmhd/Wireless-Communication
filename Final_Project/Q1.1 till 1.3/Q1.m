%% Wireless Commuincation CA3;

clc;clear all; close all;

        %% 1. Part A
% M = 2;
SNR_dB = -20:1:20;
SNR_mag = 10.^(SNR_dB/10);
len_snr = length(SNR_dB);

N = 1e6;
M = 2;
data_ = randi([0 M-1],N, 1);
txSig = pskmod(data_, M);
noise = zeros(N, len_snr);
h = (normrnd(0,1,[N,1]) + 1i*normrnd(0,1,[N,1])) /sqrt(2);
BER = zeros(len_snr, 1);
rxSig = zeros(size(txSig));
for ii = 1:len_snr
     noise(:,ii) = (((1/SNR_mag(ii))./2).^(0.5))*...
                (randn(N, 1)+1i*randn(N, 1));
     rxSig(:,ii) = txSig .* h + noise(:,ii);
     data = pskdemod(rxSig(:,ii), M);
     BER(ii) = 1 - sum(data == data_)/N;
    
end

my_semilogy(SNR_dB, BER, "Theory BER for Rayleigh Fading Channel - Non Coherent Detection", "SNR dB",...
            "BER", "BER")
ylim([0.45, 0.55])
        %% 1. Part B

BER_BPSK = qfunc(sqrt(2*SNR_mag));
my_semilogy(SNR_dB, BER_BPSK, "Theory BER - without Fading", "SNR dB",...
            "BER", "BER")


        %% 1. Part C
clc
P_e = 10^(-6);
SNR_mag = (qfuncinv(P_e)^2)/2;
SNR_dB = 10*log10(SNR_mag);
fprintf("SNR Requeired for reaching P_e = 1e-6 is %.3f", SNR_dB);

%% 2
clc;clear ; close all;
        %% 2. Part A
SNR_dB = (-20:1:20);
SNR_mag = 10.^(SNR_dB/10);
BER_theory_coherent = 1./(2.*SNR_mag + 2);
my_semilogy(SNR_dB, BER_theory_coherent,"Theory BER Non-Coherent Detection", "SNR dB", "BER");

        %% 2. Part B
clear; close all;
% Creating Data
SNR_dB = (-20:1:20);
SNR_mag = 10.^(SNR_dB/10);
N = 1e6;
M = 2;
data_ = randi([0 M-1],N, 1);
txSig = zeros(2*N,1);

for i =1 : N
    if data_(i)==1
       txSig(2*i-1:2*i) = [1 ;0];
    else
       txSig(2*i-1:2*i) = [0 ;1];
    end
end

% Transmitting Data
h = (normrnd(0,1,[2*N,1]) + 1i*normrnd(0,1,[2*N,1])) /sqrt(2);
rxSig_ = txSig .* h*(sqrt(2));
len_snr = length(SNR_dB);
rxSig = zeros(2*N, len_snr);
noise_ = zeros(2*N, len_snr);
for i = 1:len_snr
    noise_(:,i)=(((1/SNR_mag(i))./2).^(0.5))*(randn(1,2*N)+1i*randn(1,2*N));
    rxSig(:,i) =  rxSig_ + noise_(:,i) ;
end
% Reciving Signal
demodulated_rx = zeros(N, len_snr);

for j = 1 : len_snr
    for i = 1 : N
        if abs(rxSig(2*i-1 ,j))>abs(rxSig(2*i, j))
            demodulated_rx(i, j) = 1;
        else
            demodulated_rx(i, j) = 0;
        end
    end
end
BER_sim = zeros(1, len_snr);
for i =1 :len_snr
    BER_sim(i) = sum(demodulated_rx(:,i)~=data_)/(N);
end
BER_theory_coherent = 1./(2.*SNR_mag + 2);
figure;
semilogy(SNR_dB, BER_theory_coherent,'mx-')
title("BER Non-Coherent Detection")
grid on;
hold on ;
semilogy(SNR_dB ,BER_sim, 'bp-')
legend('Theory BER','Simulation BER')
xlabel("SNR dB");
ylabel("BER")

        %% 2. Part C

P_e = 10^-6;
SNR_mag = (1/(2*P_e) - 1);
SNR_dB = 10*log10(SNR_mag);
text = sprintf("SNR Requeired for reaching P_e = 1e-6 is %.3f", SNR_dB);
disp(text);

%% 3
clc; clear all; close all;


        %% 3. Part A
SNR_dB = (-10:1:10)';
SNR_mag = 10.^(SNR_dB/10);
BER_Rayleigh_FlatFading = 0.5.*(1-sqrt(SNR_mag./(SNR_mag+1)));

my_semilogy(SNR_dB, BER_Rayleigh_FlatFading, "Theory Coherent BER for Rayleigh Fading Channel", "SNR dB",...
            "BER", "BER")

        %% 3. Part B

% Creating Data
SNR_dB = (-20:1:25);
SNR_mag = 10.^(SNR_dB/10);
N = 1e6;
M = 2;
data_ = randi([0 M-1], N, 1);
txSig = pskmod(data_, M);

% Transmitting Data
h = sqrt(1/2)*(normrnd(0,1,[N,1]) + 1i*normrnd(0,1,[N,1])) ;
rxSig_ = txSig .* h;
len_snr = length(SNR_dB);
rxSig = zeros(N, len_snr);
noise_ = zeros(N, len_snr);
for i = 1:len_snr
    noise_(:,i)=(((1/SNR_mag(i))./2).^(0.5))*(randn(1,N)+1i*randn(1,N));
    rxSig(:,i) =  rxSig_ + noise_(:,i) ;
end
% Reciving Signal

demodulated_rx = zeros(N, len_snr);
for j = 1:len_snr
    for i = 1:N
        if real(rxSig(i,j)/(h(i))) < 0
            demodulated_rx(i,j) = 1;
        else
            demodulated_rx(i,j) = 0;
        end
    end
end

BER_sim = zeros(1, len_snr);
for i =1 :len_snr
    BER_sim(i) = sum(demodulated_rx(:,i)~=data_)/(N);
end

BER_theory_coherent = 0.5.*(1-sqrt(SNR_mag./(SNR_mag+1)));
BER_theory_nonCoherent = 1./(2.*SNR_mag + 2);

figure;

semilogy(SNR_dB, BER_theory_coherent,'mx-')
title("BER Coherent Detection")
grid on;
hold on ;
semilogy(SNR_dB ,BER_sim, 'bp-')
legend('Theory BER Coherent','Simulatin Coherent BER')
xlabel("SNR dB");
ylabel("BER")

figure;
semilogy(SNR_dB, BER_theory_nonCoherent,'cd-')
title("BER Coherent and NonCoherent Detection Comparison")
grid on;
hold on ;
semilogy(SNR_dB ,BER_sim, 'bp-')
legend('BER Non Coherent','BER Coherent')
xlabel("SNR dB");
ylabel("BER")


%%
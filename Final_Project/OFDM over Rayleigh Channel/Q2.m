%% Wirless Communication Second Part


clear all; clc; close all;
%% 5

T_c = 5 * 10^-3;
T_d = 10 * 10^-6;
N = 10^7;
W = 20*10^6;
L = T_d*W;

N_cp = L;

%% 2 

upper_bound = T_c * W;
lower_bound = T_d * W;


N_c = 8000;
%%


nFFT        = N_c; % fft size
nDSC        = N_c +  L; % number of data subcarriers
nBitPerSym  = N_c + L; % number of bits per OFDM symbol (same as the number of subcarriers for BPSK)
SNR_dB = (-10:1:10);
SNR_mag = 10.^(SNR_dB/10);
M = 2;




%% Creating BPSK symbols
data_ = randi([0 M-1], N, 1);
txBits = pskmod(data_, M);




%% Serial to Parallel
n_blocks = int32(N / N_c);
ifft_in = reshape(txBits, [N_c, n_blocks]);
cp_in = ifft(ifft_in, N_c, 1) ;


%% adding CP
cp_out = [cp_in(N_c-L+1:N_c, :);cp_in ];

                %% Channel:
                
rng(1);
h = sqrt(1/2)*(normrnd(0,1,[L,n_blocks]) + 1i*normrnd(0,1,[L,n_blocks])) ;
% h = ones(int32(L),int32(n_blocks));
len_snr = length(SNR_mag);
% SNR_mag = SNR_mag(10); % To be changed for performing diffrent SNRs
snr_idx = 20;
rxSig = zeros(N_c + 2*L - 1, n_blocks);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CHANNNGNNGGEEE
channel_out = zeros(N_c + 2*L - 1, n_blocks);

% for snr_idx = 1:len_snr
    for block_idx = 1:n_blocks
        channel_out(:, block_idx) = conv(cp_out(:, block_idx), h(:, block_idx));
        [noise_size, ~] = size(channel_out);
        noise_ = (((1/15)./2).^(0.5))*(randn(noise_size, 1)+ ...
        1i*randn(noise_size, 1)); % to be changed for performing diffrent SNRs
%         noise_ = zeros(noise_size, 1);
        rxSig(:, block_idx) = channel_out(:, block_idx) + noise_;
    end
    clear channel_out
    
    
    % end


    %% Removing CP:  

    rm_cp = rxSig(int32(L+1):int32(noise_size - L) + 1, :);
    % Checking Size:
    size(rm_cp)
    %% Performing FFT
    fft_out = fft(rm_cp, N_c, 1);

    %% 
    H = fft(h, N_c, 1);
    fft_out = fft_out./H;
    xx = reshape(fft_out, [N, 1]);
    final = pskdemod(xx,M);
    sum(final == data_)/N







%% Preparation part
clear all;clc;
% sys para
gi=1/4;
fftlen = 64;
gilen = gi*fftlen;
fs = 20e6;

% training seq
ShortTrain = sqrt(13/6)*[0,0,1+1i,0,0,0,-1-1i,0,0,0,1+1i,0,0,0,-1-1i,0,0,0,-1-1i,0,0,0,1+1i ...
            0,0,0,0,0,0 ,-1-1i,0,0,0,-1-1i,0,0,0,1+1i,0,0,0,1+1i,0,0,0,1+1i,0,0,0,1+1i,0,0].';
NumberSTBlks = 10;

short_train = tx_freqd_to_timed(ShortTrain);
%plot(abs(short_train));
short_train_blk = short_train(1:16);
short_train_blks = repmat(short_train_blk,NumberSTBlks,1);

LongTrain = [1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,1,1,-1,-1,1,1 ,-1,1,-1,1,1,1,1,1,-1,-1,1,1,-1,1 ...
                  -1,1,-1,-1,-1,-1,-1,1,1,-1,-1,1,-1,1,-1,1,1,1,1].';
NumberLTBlks = 2;

long_train = tx_freqd_to_timed(LongTrain);
long_train_syms = [long_train(fftlen-2*gilen+1:fftlen,:);long_train;long_train];

preamble = [short_train_blks;long_train_syms];

%% Channel
h = zeros(gilen,1);
h(1) = 1;
h(5) = 0.5;
h(10) = 0.3;
cfo = 0.2*fs/fftlen

%% Loop start
snr = 5:5:30;
mse = zeros(1,length(snr));
pkt_num = 1000;
for snr_idx = 1:length(snr)
    snr_idx
    est_err = zeros(1,pkt_num);
    for pkt_idx = 1:pkt_num
        %% Transmitter
        tx = preamble;
        %% Channel
        rx_signal = filter(h,1,tx);
        noise_var = 1/(10^(snr(snr_idx)/10))/2;
        len = length(rx_signal);
        noise = sqrt(noise_var) * (randn(len,1) + 1i*randn(len,1));
         % add noise
        rx_signal = rx_signal + noise;
         % add CFO
         total_length = length(rx_signal);
         t = [0:total_length-1]/fs;
         phase_shift = exp(1i*2*pi*cfo*t).';
         rx_signal = rx_signal.*phase_shift;
         
        %% Receiver
        pkt_det_offset = 10;
        rlen = length(short_train_blks)-pkt_det_offset;
        D = 16;
        auto_corr = rx_signal(pkt_det_offset:pkt_det_offset+rlen-D).*conj(rx_signal(pkt_det_offset+D:pkt_det_offset+rlen));
        
        % add all esti
        auto_corr_mean = sum(auto_corr);
        freq_est = -angle(auto_corr_mean)/(2*D*pi/fs);
        est_err(pkt_idx) = (freq_est-cfo)/cfo;
        %radians_per_sample = 2*pi*freq_est/fs;
        %time_base = 0:length(rx_signal)-1;
        %correction = exp(-1i*(radians_per_sample)*time_base);
        %out_signal = rx_signal.*correction.';
    end
    mse(snr_idx) = mean(abs(est_err.^2));
end
semilogy(snr,mse,'-bo');
xlabel('SNR/dB');
ylabel('MSE');

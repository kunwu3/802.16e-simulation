%% Preparation part
clear all;clc;
% sys para
gi=1/4;
fftlen = 64;
gilen = gi*fftlen;

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
%% Loop start
snr = -10:5:30;
mse = zeros(1,length(snr));
pkt_num = 1000;
for snr_idx = 1:length(snr)
    snr_idx
    est_err = zeros(1,pkt_num);
    for pkt_idx = 1:pkt_num
        %% Transmitter
        tx = preamble;
       %% Channel
        h = zeros(gilen,1);
        h(1) = 1;
        h(5) = 0.5;
        h(10) = 0.3;
        % index define
        UsedSubcIdx = [7:32 34:59];
        reorder = [33:64 1:32];

        channel = fft(h,64);
        channel(reorder) = channel;
        channel = channel(UsedSubcIdx);

        %% Transmitter
        tx = preamble;

        %% Pass channel
        rx_signal = filter(h,1,tx);
        noise_var = 1/(10^(snr(snr_idx)/10))/2;
        len = length(rx_signal);
        noise = sqrt(noise_var) * (randn(len,1) + 1i*randn(len,1));
        % add noise
        rx_signal = rx_signal + noise;
        
        %% Channel estimation
        rx_long = rx_signal(161:end);
        long_tr_syms = rx_long(33:end);
        long_tr_syms = reshape(long_tr_syms,64,2);

        % to frequency domain
        freq_long_tr = fft(long_tr_syms)/(64/sqrt(52));
        freq_long_tr(reorder,:) = freq_long_tr;
        freq_tr_syms = freq_long_tr(UsedSubcIdx,:);

        channel_est = mean(freq_tr_syms,2).*conj(LongTrain);
        err_est(pkt_idx) = mean(abs(channel_est-channel).^2)/mean(abs(channel).^2);
    end
    mse(snr_idx) = mean(abs(err_est).^2);
end
semilogy(snr,mse);
xlabel('SNR/dB');
ylabel('MSE');
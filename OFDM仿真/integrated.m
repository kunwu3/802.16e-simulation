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

%% Loop start
snr = 8:2:30;
pkt_num = 1000;

freq_mse = zeros(1,length(snr));
time_mse = zeros(1,length(snr));
channel_mse = zeros(1,length(snr));

for snr_idx = 1:length(snr)
    snr_idx
    freq_est_err = zeros(1,pkt_num);
    time_est_err = zeros(1,pkt_num);
    channel_est_err = zeros(1,pkt_num);
    for pkt_idx = 1:pkt_num
        %% Transmitter
        tx = preamble;
        %% Channel
        h = zeros(gilen,1);
        h(1) = 1;
        h(5) = 0.5;
        h(10) = 0.3;
        start_noise_len = 500;
        cfo = 0.2*fs/fftlen;
        ideal_start = 193;

        % index define
        UsedSubcIdx = [7:32 34:59];
        reorder = [33:64 1:32];
        channel = fft(h,64);
        channel(reorder) = channel;
        channel = channel(UsedSubcIdx);
        
        rx_signal = filter(h,1,tx);
        noise_var = 1/(10^(snr(snr_idx)/10))/2;
        len = length(rx_signal);
        noise = sqrt(noise_var) * (randn(len,1) + 1i*randn(len,1));
        start_noise = sqrt(noise_var) * (randn(start_noise_len,1) + 1i*randn(start_noise_len,1));
         % add noise
        rx_signal = rx_signal + noise;
        rx_signal = [start_noise;rx_signal];
         % add CFO
         total_length = length(rx_signal);
         t = [0:total_length-1]/fs;
         phase_shift = exp(1i*2*pi*cfo*t).';
         rx_signal = rx_signal.*phase_shift;
         
        %% Receiver
        % Packet detect
        search_win = 700;
        D = 16;
        
        delay_xcorr = rx_signal(1:search_win).*conj(rx_signal(D+1:search_win+D));
        ma_delay_xcorr = abs(filter(ones(1,D),1,delay_xcorr));
        ma_rx_pwr = filter(ones(1,D),1,abs(rx_signal(D+1:search_win+D)).^2);
        delay_len = length(ma_delay_xcorr);
        ma_M = ma_delay_xcorr(1:delay_len)./ma_rx_pwr(1:delay_len);
        ma_M(1:16) = [];

        threshold = 0.75;
        pre_thres_idx = find(ma_M > threshold);
        
        detect_win = 12;
        detect_len = length(pre_thres_idx)-detect_win;
        for detect_idx = 1:detect_len
            detect_thres(detect_idx) = mean(ma_M(pre_thres_idx(detect_idx):pre_thres_idx(detect_idx+detect_win)));
            if detect_thres(detect_idx) > threshold
                break
            end
        end
        thres_idx = pre_thres_idx(detect_idx);

        detected_packet = rx_signal(thres_idx:length(rx_signal));
        % Frequency sync
        pkt_det_offset = 10;
        rlen = length(short_train_blks)-pkt_det_offset;
        
        auto_corr = detected_packet(pkt_det_offset:pkt_det_offset+rlen-D).*conj(detected_packet(pkt_det_offset+D:pkt_det_offset+rlen));
        auto_corr_mean = sum(auto_corr);
        freq_est = -angle(auto_corr_mean)/(2*D*pi/fs);
        freq_est_err(pkt_idx) = (freq_est-cfo)/cfo;
        radians_per_sample = 2*pi*freq_est/fs;
        time_base = 0:length(detected_packet)-1;
        correction = exp(-1i*(radians_per_sample)*time_base);
        out_signal = detected_packet.*correction.';
        % Fine Time sync
        start_search = 150;
        end_search = 200;
        time_corr_long = zeros(1,end_search-start_search+1);
        
        for idx = start_search:end_search
            time_corr_long(idx-start_search+1) = sum((out_signal(idx:idx+63).*conj(long_train)));
        end
        
        [max_corr,long_search_idx] = max(abs(time_corr_long));
        
        fine_time_est = start_search-1 + long_search_idx;
        time_err_est(pkt_idx) = fine_time_est - ideal_start;
        % Channel Estimation
        long_tr_syms = out_signal(fine_time_est:end);
        long_tr_syms = reshape(long_tr_syms,64,2);

        % to frequency domain
        freq_long_tr = fft(long_tr_syms)/(64/sqrt(52));
        freq_long_tr(reorder,:) = freq_long_tr;
        freq_tr_syms = freq_long_tr(UsedSubcIdx,:);

        channel_est = mean(freq_tr_syms,2).*conj(LongTrain);
        channel_err_est(pkt_idx) = mean(abs(channel_est-channel).^2)/mean(abs(channel).^2);
    end
    freq_mse(snr_idx) = mean(abs(freq_est_err.^2));
    time_mse(snr_idx) = mean(abs(time_est_err.^2));
    channel_mse(snr_idx) = mean(abs(channel_est_err.^2));
end
subplot(1,3,1);
xlabel('SNR/dB');
ylabel('MSE'); 
plot(snr,freq_mse)
title('freq mse')
subplot(1,3,2);
xlabel('SNR/dB');
ylabel('MSE');
plot(snr,time_mse)
title('time mse')
subplot(1,3,3);
xlabel('SNR/dB');
ylabel('MSE');
plot(snr,channel_mse)
title('channel mse')
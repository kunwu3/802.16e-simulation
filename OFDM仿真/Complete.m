clear all; clc;

%% ************** System parameters ********************
ml = 2;                      % Modulation level: 2--4QAM; 4--16QAM; 6--64QAM
NormFactor = sqrt(2/3*(ml.^2-1));
gi = 1/4;                   % Guard interval
fftlen = 64;
gilen = gi*fftlen;           % Length of guard interval (samples)
blocklen = fftlen + gilen;   % total length of the block with CP
non_fftlen = 52;
% conv
trellis = poly2trellis(7,[133 171]);
code_rate = 1/2;
tb = 7*5;

bits_per_sym = code_rate*non_fftlen*ml;
NumSyms = 50;
TotalNumBits = bits_per_sym*NumSyms;
fs = 20e6;
% training seq
ShortTrain = sqrt(13/6)*[0,0,1+1i,0,0,0,-1-1i,0,0,0,1+1i,0,0,0,-1-1i,0,0,0,-1-1i,0,0,0,1+1i ...
            0,0,0,0,0,0 ,-1-1i,0,0,0,-1-1i,0,0,0,1+1i,0,0,0,1+1i,0,0,0,1+1i,0,0,0,1+1i,0,0].';
NumberSTBlks = 10;

short_train = tx_freqd_to_timed(ShortTrain);
short_train_blk = short_train(1:16);
short_train_blks = repmat(short_train_blk,NumberSTBlks,1);

LongTrain = [1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,1,1,-1,-1,1,1 ,-1,1,-1,1,1,1,1,1,-1,-1,1,1,-1,1 ...
                  -1,1,-1,-1,-1,-1,-1,1,1,-1,-1,1,-1,1,-1,1,1,1,1].';
NumberLTBlks = 2;

long_train = tx_freqd_to_timed(LongTrain);
long_train_syms = [long_train(fftlen-2*gilen+1:fftlen,:);long_train;long_train];

preamble = [short_train_blks;long_train_syms];
% index define
%DataSubcPatt = [1:5 7:19 21:26 27:32 34:46 48:52]';
%DataSubcIdx = [7:11 13:25 27:32 34:39 41:53 55:59];
UsedSubcIdx = [7:32 34:59];
reorder = [33:64 1:32];

%% *******Loop start***********
snr = 10:5:30;
ber = zeros(1,length(snr));
pkt_num = 1000;
for snr_idx = 1:length(snr)
    err = 0;
    snr_idx
    for pkt_idx = 1:pkt_num
        %% *********Transmitter*********
        % Generate the information bits
        inf_bits = randn(1,TotalNumBits)>0;
        
        % Encoding
        coded_bits = convenc(inf_bits,trellis);
        
        % Interleave
        interleaved_bits = tx_interleaver(coded_bits,non_fftlen, ml);
        
        % Modulate
        paradata = reshape(interleaved_bits,length(interleaved_bits)/ml,ml);
        mod_ofdm_syms = qammod(bi2de(paradata),2^ml)./NormFactor;

        % Map&IFFT
        mod_ofdm_syms = reshape(mod_ofdm_syms,non_fftlen,NumSyms);
        mod_ofdm_syms = [zeros(6,NumSyms);mod_ofdm_syms(1:26,1:NumSyms);zeros(1,NumSyms);mod_ofdm_syms(27:52,1:NumSyms);zeros(5,NumSyms)];
        mod_ofdm_syms = mod_ofdm_syms(reorder,:);
        tx_blks = sqrt(fftlen)*ifft(mod_ofdm_syms);

        % Add CP
        tx_frames = [tx_blks(fftlen-gilen+1:fftlen,:); tx_blks];

        % P/S
        trans_bits = reshape(tx_frames,NumSyms*blocklen,1);
        
        % Transmit
        tx = [preamble;trans_bits];

        %% ************Channel*************
        h = zeros(gilen,1);
        h(1) = 1;
        h(5) = 0.5;
        h(10) = 0.3;
        start_noise_len = 500;
        cfo = 0.2*fs/fftlen;
        ideal_start = 193;
                
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
         
        %% ************Receiver************
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
        % Channel Estimation
        long_tr_syms = out_signal(fine_time_est:fine_time_est+127);
        long_tr_syms = reshape(long_tr_syms,64,2);

        freq_long_tr = fft(long_tr_syms)/(64/sqrt(52)); % to frequency domain
        freq_long_tr(reorder,:) = freq_long_tr;
        freq_tr_syms = freq_long_tr(UsedSubcIdx,:);

        channel_est = mean(freq_tr_syms,2).*conj(LongTrain);
        
        % remove preamble
        data_syms = out_signal(fine_time_est+128:end);
        data_syms = reshape(data_syms,80,NumSyms);
        
        % remove gi
        data_syms(1:16,:) = [];
        freq_data = fft(data_syms)/(64/sqrt(52));
        freq_data(reorder,:) = freq_data;
        
        % select data carriers
        freq_data_syms = freq_data(UsedSubcIdx,:);
        
        % channel equalization
        chan_corr_mat = repmat(channel_est,1,size(freq_data_syms,2));
        freq_data_syms = freq_data_syms.*conj(chan_corr_mat);
        chan_sq_amp = abs(channel_est).^2;
        chan_sq_amp_mtx = repmat(chan_sq_amp,1,size(freq_data_syms,2));
        data_syms_out = freq_data_syms./chan_sq_amp_mtx;
        
        Data_seq = reshape(data_syms_out,52*NumSyms,1).*NormFactor;
        
        % demodulate
        DemodSeq = de2bi(qamdemod(Data_seq,2^ml),ml);
        SerialBits = reshape(DemodSeq,size(DemodSeq,1)*ml,1).';
        % De-interleave
        deint_bits = rx_deinterleave(SerialBits, non_fftlen, ml);
        % decode
        hard_decision = deint_bits > 0;
        DecodedBits = vitdec(hard_decision,trellis,tb,'trunc','hard');
        err = err+sum(abs(DecodedBits-inf_bits));
    end
        ber(snr_idx) = err/(TotalNumBits*pkt_num);
end
semilogy(snr,ber);
xlabel('SNR/dB');
ylabel('BER');
grid;hold on;
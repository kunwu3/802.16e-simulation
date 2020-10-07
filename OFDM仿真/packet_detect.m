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

%% Channel
h = zeros(gilen,1);
h(1) = 1;
h(5) = 0.5;
h(10) = 0.3;
start_noise_len = 500;
snr = 50;

%% Transmitter
tx = [short_train_blks;long_train_syms];

%% Pass channel
rx_signal = filter(h,1,tx);
noise_var = 1/(10^(snr/10))/2;
len = length(rx_signal);
noise = sqrt(noise_var) * (randn(len,1) + 1i*randn(len,1));
% add noise
rx_signal = rx_signal + noise;
start_noise = sqrt(noise_var) * (randn(start_noise_len,1) + 1i*randn(start_noise_len,1));
rx_signal = [start_noise;rx_signal];


%% Receiver
search_win = 700;
D = 16;

% Calculate delay correlation
delay_xcorr = rx_signal(1:search_win).*conj(rx_signal(D+1:search_win+D));
% Moving average of the delayed correlation
ma_delay_xcorr = abs(filter(ones(1,D),1,delay_xcorr));
% Moving average of received power
ma_rx_pwr = filter(ones(1,D),1,abs(rx_signal(D+1:search_win+D)).^2);
% The decision variable
delay_len = length(ma_delay_xcorr);
ma_M = ma_delay_xcorr(1:delay_len)./ma_rx_pwr(1:delay_len);
% Remove delay samples
ma_M(1:16) = [];

plot(ma_M);grid;

threshold = 0.75;
%thres_idx = find(ma_M > threshold);

%if isempty(thres_idx)
    %thres_idx = 1;
%else
    %thres_idx = thres_idx(1)
%end
pre_thres_idx = find(ma_M > threshold);        
detect_win = 12;
detect_len = length(pre_thres_idx)-detect_win;
for detect_idx = 1:detect_len
detect_thres(detect_idx) = mean(ma_M(pre_thres_idx(detect_idx):pre_thres_idx(detect_idx+detect_win)));
    if detect_thres(detect_idx) > threshold
        break
    end
end
thres_idx = pre_thres_idx(detect_idx)

detected_packet = rx_signal(thres_idx:length(rx_signal));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                %
% Uplink and Downlink Synchronization of 802.16                  %
%                                                                %
% by Qing Wang                                                   %
%                                                                %
% San Diego State University                                     %
%                                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% Downlink model %%%%%%%%%%%%%%%

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi_spin = 3/360;   % phase offset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BW = 2.4;           % Bandwidth of channel(GHz)
G = 1/8;            % ratio of CP
N_symbol = 100;     % number of OFDM symbols
dly = zeros(1,80);  % random delay
Nfft = 256;         % Size of FFT/IFFT
chan_flg = 1;       % if flg = 1, channel is activated (=0 AWGN channel)
chan_type = 1;      % if chan_type = 1, multipath channel, =2, SUI channel
N_SUI = 3;          % six SUI channel model
cyc_pfx = 1;        % if flg = 1, cyclic prefix is activated
noise_flg = 1;      % if flg = 1, AWGN is added
spin_flg = 1;       % if flg = 1, phase offset is activated
equ_flg = 1;        % if flg = 1, equalizer is activated
N = 256;            % Number of subcarriers
Ng = N*G;           % Length of cyclic prefix
Ns = Ng + N;        % Symbol length (CP + useful part)
num_bins = 200;     % Number of data subcarriers in OFDM symbol
DC = 0;             % DC bin
c = 1/sqrt(10);     % multiplying with the constellation point
                    % to achieve equal average power (16-QAM)
guard_interval_l = zeros(1, 0.5*(N - num_bins));   % left guard band
guard_interval_r = zeros(1, 0.5*(N - num_bins)-1); % right guard band

%generate pilots that to be inserted;
%here we use DL, if UL, seed = [1,0,1,0,1,0,1,0,1,0,1];
seed = [1,1,1,1,1,1,1,1,1,1,1];
%seed = [1,0,1,0,1,0,1,0,1,0,1];

for i = 1:N_symbol+2
    wk = seed (11);
    next = xor(seed(9),seed(11));
    seed = [next seed(1,1:10)];
end

A = 1 - 2*wk; % Values defined in the standard.
B = 1 - 2*(~wk);
value_carrier = [A,B,A,B,B,B,A,A];
%value_carrier = [A,B,A,B,A,A,A,A];
% if UL, the values should be [A,B,A,B,A,A,A,A]

pilots = value_carrier.*2;
%The factor of "2" is due to the fact that the pilots are
%transmitted to a double power of the information bits.

% in order to distract data subcarriers at receiver, need to know the
% pilot and data subcarriers positions;
v_pilots = [41,66,91,116,142,167,192,217];
v_data = setxor(1:N,v_pilots);

% generating 50 OFDM symbols
xx = []; % xx is the OFDM time series output [50*(32+256)]
tx_data = [];
for nn1 = 1:N_symbol
    % generate binary data and 16-QAM modulation
    rng('shuffle')
    tx = randi([0,1],1,4*(num_bins-length(pilots)));
    qamdata = bi2de(reshape(tx,4,num_bins-length(pilots)).','left-msb').';
    maping = bin2gray(qamdata,'qam',16);
    data = c*qammod(maping,16);
    % combine guard carriers, data carriers and pilot carriers;
    fx = [guard_interval_l,data(1:12),pilots(1),data(13:36),pilots(2), ...
        data(37:60),pilots(3),data(61:84),pilots(4),data(85:96),DC, ...
        data(97:108),pilots(5),data(109:132),pilots(6),data(133:156), ...
        pilots(7),data(157:180),pilots(8),data(181:192),guard_interval_r];
    x = ifft(fx,Nfft);
    if (cyc_pfx == 0)
        x_cyc = [zeros(1,Ng),x];
    else
        x_cyc = [x(N-Ng+1:N),x];
    end
    % original serial binary data to be used for BER calculation
    tx_data = [tx_data,tx];
    xx = [xx x_cyc];
end
xx_cyc_ser = xx;

% add long preamble
P_ALL = [1-1i,1-1i,-1-1i,1+1i,1-1i,1-1i,-1+1i,1-1i,1-1i, ...
    1-1i,1+1i,-1-1i,1+1i,1+1i,-1-1i,1+1i,-1-1i,-1-1i, ...
    1-1i,-1+1i,1-1i,1-1i,-1-1i,1+1i,1-1i,1-1i,-1+1i, ...
    1-1i,1-1i,1-1i,1+1i,-1-1i,1+1i,1+1i,-1-1i,1+1i, ...
    -1-1i,-1-1i,1-1i,-1+1i,1-1i,1-1i,-1-1i,1+1i,1-1i, ...
    1-1i,-1+1i,1-1i,1-1i,1-1i,1+1i,-1-1i,1+1i,1+1i, ...
    -1-1i,1+1i,-1-1i,-1-1i,1-1i,-1+1i,1+1i,1+1i,1-1i, ...
    -1+1i,1+1i,1+1i,-1-1i,1+1i,1+1i,1+1i,-1+1i,1-1i, ...
    -1+1i,-1+1i,1-1i,-1+1i,1-1i,1-1i,1+1i,-1-1i,-1-1i, ...
    -1-1i,-1+1i,1-1i,-1-1i,-1-1i,1+1i,-1-1i,-1-1i, ...
    -1-1i,1-1i,-1+1i,1-1i,1-1i,-1+1i,1-1i,-1+1i,-1+1i, ...
    -1-1i,1+1i,0,-1-1i,1+1i,-1+1i,-1+1i,-1-1i,1+1i, ...
    1+1i,1+1i,-1-1i,1+1i,1-1i,1-1i,1-1i,-1+1i,-1+1i, ...
    -1+1i,-1+1i,1-1i,-1-1i,-1-1i,-1+1i,1-1i,1+1i, ...
    1+1i,-1+1i,1-1i,1-1i,1-1i,-1+1i,1-1i,-1-1i,-1-1i, ...
    -1-1i,1+1i,1+1i,1+1i,1+1i,-1-1i,-1+1i,-1+1i,1+1i, ...
    -1-1i,1-1i,1-1i,1+1i,-1-1i,-1-1i,-1-1i,1+1i,-1-1i, ...
    -1+1i,-1+1i,-1+1i,1-1i,1-1i,1-1i,1-1i,-1+1i,1+1i, ...
    1+1i,-1-1i,1+1i,-1+1i,-1+1i,-1-1i,1+1i,1+1i,1+1i, ...
    -1-1i,1+1i,1-1i,1-1i,1-1i,-1+1i,-1+1i,-1+1i,-1+1i, ...
    1-1i,-1-1i,-1-1i,1-1i,-1+1i,-1-1i,-1-1i,1-1i,-1+1i, ...
    -1+1i,-1+1i,1-1i,-1+1i,1+1i,1+1i,1+1i,-1-1i,-1-1i, ...
    -1-1i,-1-1i,1+1i,1-1i,1-1i];

%%%%%%%% insert 2nd part of long preamble %%%%%%%
for k = 1:num_bins+1 % -100 ~ +100
    if mod(k,2) == 1 % even bins
        dat_l(k) = sqrt(2)*P_ALL(k);
    else
        dat_l(k) = 0;
    end
end
fx_l = [guard_interval_l,dat_l,guard_interval_r];
xll = ifft(fx_l,Nfft); % [128,128] pattern
xl = [xll(N-Ng+1:N),xll]; % add 32 bit cycle prefix
ofdm_tx = [xl,xx_cyc_ser]; % [(32+256)_l,50*(32+256)_data]

%%%%%%%%% insert 1st part of long preamble %%%%%%
for k = 1:num_bins+1
    if mod(k,4) == 1 % every 4th bin
        dat_s(k) = sqrt(2)*sqrt(2)*conj(P_ALL(k));
    else
        dat_s(k) = 0;
    end
end
fx_s = [guard_interval_l,dat_s,guard_interval_r];
xss = ifft(fx_s,Nfft); %[64,64,64,64] pattern
xs = [xss(N-Ng+1:N),xss]; % add 32 bit cycle prefix
ofdm_tx = [dly,xs,ofdm_tx]; % [dly,(32+256)_s,(32+256)_l,50*(32+256)_data]

figure
plot(-127:-64,real(xss(1:64)))
hold on
plot(-64*ones(1,length(-2:2)),(-2:2),'m','linewidth',2)
plot(-63:0,real(xss(65:128)),'r')
plot(0*ones(1,length(-2:2)),(-2:2),'m','linewidth',2)
plot(1:64,real(xss(129:192)))
plot(64*ones(1,length(-2:2)),(-2:2),'m','linewidth',2)
plot(65:128,real(xss(193:256)),'r')
hold off
axis([-127,128,-0.15,0.15]); grid on
title('The 1st symbol of long preamble, 4-replicate')

figure
plot(-127:0,real(xll(1:128)))
hold on
plot(1:128,real(xll(129:256)),'r')
plot(0*ones(1,length(-2:2)),(-2:2),'m','linewidth',2)
hold off
axis([-127,128,-0.15,0.15]); grid on
title('The 2nd symbol of long preamble, 2-replicate')

figure
subplot(211)
plot(real(ofdm_tx(1:1000)))
axis([0,1000,-0.15,0.15]);grid on
title('Real Part of The OFDM Data')
xlabel('Samples / n')
ylabel('Amplitude')
subplot(212)
plot(imag(ofdm_tx(1:1000)))
axis([0,1000,-0.15,0.15]);grid on
title('Imag Part of The OFDM Data')
xlabel('Samples / n')
ylabel('Amplitude')

xx_tr = [];
for tt = 1:Ns:Ns*N_symbol
    pp = xx_cyc_ser(tt:tt+Ns-1);
    pp2 = pp(Ng+1:Ns); % discard cyclic prefix
    fpp = fft(pp2,Nfft);
    xx_tr = [xx_tr,fpp];
end

xx_par = reshape(xx_tr,N,N_symbol).';
xx_par_data_with_guard = xx_par(:,v_data);
xx_par_data = xx_par_data_with_guard(:,29:221);

figure
subplot(211)
plot(0,0)
hold on
for n = N_symbol-9:N_symbol
plot(xx_par_data(n,:),'ro')
end
hold off
axis([-1.5,1.5,-1.5,1.5])
axis('square')
grid on
title('Signal Constellation -- Before pass the channel')

subplot(212)
plot(0,0)
hold on
for i = N_symbol-9:N_symbol
plot([-100+4:100-4],real(xx_par_data(i,:)),'ro')
end
hold off
axis([-100,100,-2,2]);
grid on
xlabel('Frequency Index'); ylabel('Amplitude')
title('Spectrum Lines -- In Phase -- Before pass the channel')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pass channel
% input xx_channel_in
xx_channel_in = ofdm_tx;
if (chan_flg == 1)
    if (chan_type == 1)
        hc = [1 0 0 0 0 0.5 0 0.2]; % multipath channel
    elseif (chan_type == 2)
        hc = channelSUI(G,BW,N_SUI); % SUI channel
    end
else
    hc = [1 0];
end
channel_out = filter(hc,1,xx_channel_in);

% add phase offset
if(spin_flg == 1)
    ofdm_spin = channel_out.*exp(1i*2*pi*phi_spin*(1:length(channel_out)));
else
    ofdm_spin = channel_out;
end

if chan_flg == 1
    v_EbNo_dB = 4:2:14;
else
    v_EbNo_dB = 2:2:12;
end
snr = v_EbNo_dB + 10*log10(num_bins/N) + 10*log10(4);

no_error = [];
v_ber = [];

for ii = 1:length(snr)
    % noise
    if (noise_flg == 1)
        ofdm_rx = awgn(ofdm_spin, snr(ii), 'measured');
    else
        ofdm_rx = ofdm_spin;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% at receiver
ds_input = ofdm_rx;
delta = .01;
reg = zeros(1,65); % delay register
reg_cross = zeros(1,64); % cross
reg_auto = zeros(1,64); % auto
cic_reg = ones(1,64);
for nn2 = 1:2*N
    reg = [ds_input(nn2),reg(1:64)];
    
    r64 = reg(1)*conj(reg(65)); % cross correlation
    r64_sv(nn2) = r64;
    reg_cross = [r64 reg_cross(1:63)];
    r_64_hat = cic_reg*reg_cross'; % averaged cross correlation
    r_64_hat_sv(nn2) = r_64_hat;
    r0 = reg(65)*conj(reg(65)); % auto correlation
    r0_sv(nn2) = r0;
    reg_auto = [r0 reg_auto(1:63)];
    r_0_hat = cic_reg*reg_auto'; % averaged auto correlation
    r_0_hat_sv(nn2) = r_0_hat;
    c64 = r_64_hat/(r_0_hat+delta); % cross over auto
    c64_sv(nn2) = c64;
end

% de-spinner

phi1 = phase(sum(r_64_hat_sv(N:N+64))/length(r_64_hat_sv(N:N+64)));
phi2 = -(2*pi-phi1)/64;

if(spin_flg == 1)
    ds_output = exp(1i*phi2*(1:length(ds_input))).*ds_input;
else
    ds_output = ds_input;
end

% Frame detection
% correlation between received preamble and local copy of preamble
sync_seq = ds_output;
P = zeros(1,2*Ns);
R = zeros(1,2*Ns);
M = zeros(1,2*Ns);

% d - time index corresponding to the 1st sample in a window of N/4+Ng samples
for d = 1:2*Ns
    for m = 0:1:N/4+Ng-1
        P(d) = P(d) + sync_seq(d+m)*conj(sync_seq(d+m+N/4)) ...
            + sync_seq(d+m)*conj(sync_seq(d+m+2*N/4)) ...
            + sync_seq(d+m)*conj(sync_seq(d+m+3*N/4));
        
        R(d) = R(d) + power(abs(sync_seq(d+m+N/4)),2) ...
            + power(abs(sync_seq(d+m+2*N/4)),2)...
            + power(abs(sync_seq(d+m+3*N/4)),2);
    end
end
Metric = power(abs(P),2)./power(R,2);
agc = max(Metric); % Output normalizing factor for automatic gain control

% AGC (auto gain control)
ds_output = ds_output./agc;

% t - time index corresponding to the start of 2nd symbol of long preamble
local_long_preamble = xll(1:N/2); % local copy of the long preamble at Rx
for t = 1:2*Ns
    for m = 0:1:N/2-1
        M(t) = M(t) + sync_seq(t+m)*conj(local_long_preamble(m+1));
    end
end
ind_ll = find(abs(M) >= 0.95*max(abs(M)));
% excluding cyclic prefix, 1st and 129th in 2nd symbol of long preamble
% are both the peak value, thus, we choose the 1st one as the symbol start
if length(ind_ll) > 1
    ind_l = ind_ll(length(ind_ll)-1);
else
    ind_l = ind_ll;
end
ind_data = ind_l + N;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rx_eq = ds_output(ind_data:end);
zz = []; %% zz is the output of fft
for nn = 1:Ns:Ns*N_symbol
    yy = rx_eq(nn:nn+Ns-1);
    yy2 = yy(Ng+1:Ns); % disgard cyclic prefix
    fyy = fft(yy2,Nfft);
    zz = [zz,fyy];
end

% recovery data -- 50 symbols
rx_par = reshape(zz,N,N_symbol).'; % bin -128 - bin +128;

% FEQ (frequency domain equalizer) (use 2nd symbol of the long preamble)
rx_ll = ds_output(ind_l:ind_l+N-1);
rx_l = fft(rx_ll,Nfft);
% interpolation to get the channel estimation
if (equ_flg == 1)
    chan_res = rx_l;
    % even bins
    for k = 29:2:127 %(-100,-98,...-2)
        chan_res(k) = rx_l(k)./fx_l(k);
    end
    for k = 131:2:229 % (+2,...+98,+100)
        chan_res(k) = rx_l(k)./fx_l(k);
    end
    % note that DC bin is at 129;
    chan_res(129) = 0.5*(chan_res(127)+chan_res(131));
    % odd bins (-99,-97,...-1,+1,...+97,+99)
    for k = 30:2:228
        chan_res(k) = 0.5 * ( chan_res(k-1) + chan_res(k+1) );
    end
else
    chan_res = 1;
end

for n = 1:N_symbol
    ofdm_equ_out(n,:) = rx_par(n,:)./chan_res;
end

% extract the data carriers by removing the pilot carriers;
temp = ofdm_equ_out(:,v_data);
pilot_temp = ofdm_equ_out(:,v_pilots);
rx_par2 = rx_par(:,v_data);
rx_par3 = rx_par2(:,29:221);

% calculating phase factor by using pilot carriers;
for n = 1:N_symbol
    pilot_est(n,:) = pilot_temp(n,:)./pilots;
    phi(n) = phase(sum(pilot_est(n,:))/length(pilot_est(n,:)));
end
% correcting the phase factor caused by the frequency offset residual
if spin_flg == 1
    for n = 1:N_symbol
        ofdm_equ_data(n,:) = temp(n,:)*exp(-1i*phi(n));
    end
else
ofdm_equ_data = temp;
end
% remove the guard carriers
% ofdm_equ_data is the demodulated, de-packed and equalized signal
ofdm_equ_data2 = ofdm_equ_data(:,29:221);

% remove the DC bin
ofdm_equ_data3 = [ofdm_equ_data2(:,1:96),ofdm_equ_data2(:,98:193)];
% demodulate the received data to binary data
rx = reshape(ofdm_equ_data3.',1,numel(ofdm_equ_data3));
rx_demod = qamdemod((1/c)*rx,16);
demaping = gray2bin(rx_demod,'qam',16);
rx_bi = de2bi(demaping,'left-msb');
rx_data = reshape(rx_bi.',1,numel(rx_bi));
% calculate the bit error rate
[no_of_error(ii), v_ber(ii)] = biterr(tx_data, rx_data);

end

%%%%% figures %%%%%
link_flg = 0; % downlink
%v_ber_perfect = perfect_sync(link_flg,dly,fx_l,N_symbol,tx_data,ofdm_tx,chan_flg,hc,snr);
%v_ber_no_sync = no_sync(link_flg,dly,fx_l,phi_spin,N_symbol,tx_data,ofdm_tx,chan_flg,hc,snr);

% power spectrum of the OFDM
win = kaiser(length(ofdm_rx),10);
X1 = fft(ofdm_rx.*win',2048);
miny = min(20*log10(abs(X1)));
maxy = max(20*log10(abs(X1)));
figure
plot([-0.5:1/2048:0.5-1/2048],20*log10(abs(X1)))
grid on
xlabel('Normalized Frequency')
ylabel('Magnitude / dB')
title('OFDM Power Spectrum--Before De-rotation')
axis([-0.5,0.5,miny+10,maxy+5])
hold on
plot(phi_spin*ones(1,length(miny+10:maxy+5)),(miny+10:maxy+5),'r--','linewidth',2)
hold off

figure
ANG = angle(conv(r64_sv,ones(1,64)));
plot(ANG)
miny = min(ANG);
maxy = max(ANG);
title('Estimated Rotation Angle -- phase of the cross-correlation')
grid on
xlabel('Samples / n')
ylabel('Estimated Angle / Radian')
axis([0,Ns,miny-0.5,maxy+0.5])
if spin_flg == 1
    hold on
    plot((0:length(r64_sv)),-phi_spin*360*ones(1,length(0:length(r64_sv))),'r--')
    hold off
    text(180,-phi_spin*360-0.35,'CFO = 3/360','VerticalAlignment','bottom','fontsize',10,'Interpreter','latex');
end

figure
subplot(211)
plot(abs(r_64_hat_sv))
hold on
plot(r_0_hat_sv,'r--')
hold off
grid on
title('Auto and Cross Correlation');
legend('Cross Correlation','Auto Correlation');
set(legend,'Location','NorthWest');
xlabel('Samples / n'); ylabel('Amplitude')
axis([0,1.5*Ns,-0.05,max(r_0_hat_sv)+0.05])

subplot(212)
plot(abs(r_64_hat_sv))
hold on
plot(r_0_hat_sv,'r--')
hold off
grid on
title('Zoom In');
legend('Cross Correlation','Auto Correlation');
set(legend,'Location','NorthWest');
xlabel('Samples / n'); ylabel('Amplitude')
axis([150,220,-0.05,max(r_0_hat_sv)+0.05])

figure
subplot(211)
plot(abs(c64_sv))
x_axis = find(abs(c64_sv)>=0.5);
grid on
title('Ratio of Cross over Auto Correlation')
axis([0,1.5*Ns,-0.1,max(abs(c64_sv))+0.1])
xlabel('Samples / n'); ylabel('Amplitude')
hold on
plot((x_axis(1)-20:x_axis(1)+20),0.5*ones(1,41),'r--','linewidth',2)
hold off
text(165,0.55,'threshold = 0.5','HorizontalAlignment','left','fontsize',10,'Interpreter','latex')

subplot(212)
plot(abs(c64_sv))
grid on
title('Zoom In')
axis([120,170,-0.05,max(abs(c64_sv))+0.05])
xlabel('Samples / n'); ylabel('Amplitude')
hold on
plot((x_axis(1)-3:x_axis(1)+3),0.5*ones(1,7),'r--','linewidth',2)
plot(x_axis(1)*ones(1,length(-0.05:max(abs(c64_sv))+0.05)),(-0.05:max(abs(c64_sv))+0.05),'r--')
hold off
text(x_axis(1)+4,0.55,'threshold = 0.5','HorizontalAlignment','left','fontsize',10,'Interpreter','latex')

% power spectrum of the OFDM
win = kaiser(length(rx_eq),10);
X2 = fft(rx_eq.*win',2048);
miny2 = min(20*log10(abs(X2)));
maxy2 = max(20*log10(abs(X2)));
figure
plot([-0.5:1/2048:0.5-1/2048],20*log10(abs(X2)))
grid on
xlabel('Normalized Frequency')
ylabel('Magnitude / dB')
title('OFDM Power Spectrum--After De-rotation')
axis([-0.5,0.5,miny2+10,maxy2+5])
hold on
plot(DC*ones(1,length(miny2+10:maxy2+5)),(miny2+10:maxy2+5),'r--','linewidth',2)
hold off

figure
d = 1:2*Ns;
plot(d,M);
grid on;
title('Cross Correlation with long preamble');
xlabel('Time (sample)'); ylabel('Timing Metric');

figure
subplot(211)
plot(0,0)
hold on
for i = N_symbol-9:N_symbol
    plot(rx_par3(i,:),'ro')
end
plot(data,'bx','linewidth',2);
hold off
title('Overlaid Demodulated Signal Constellation -- Before Equalizer ')
grid on
axis([-1.5,1.5,-1.5,1.5])
axis('square')

subplot(212)
plot(0,0)
hold on
for i = N_symbol-9:N_symbol
    plot([-100+4:100-4],real(rx_par3(i,:)),'ro')
end
hold off
title('Spectrum Lines -- In Phase -- Before Equalizer')
axis([-100,100,-2,2]);
grid on
xlabel('Frequency Index'); ylabel('Amplitude')

figure % channel response and estimator
plot(-0.5:1/256:0.5-1/256, 20*log10(abs(fft(hc,256))),'b','linewidth',2);
hold on
plot(-0.5:1/256:0.5-1/256, 20*log10(abs(chan_res)),'rx-');
hold off
xlabel('Normalized Frequency'); ylabel('Amplitude /dB')
legend('actual channel','channel estimator')
set(legend,'Interpreter','latex','Location','South'); grid on

figure
subplot(211)
plot(0,0)
hold on
for i = N_symbol-9:N_symbol
    plot(temp(i,29:221),'ro')
end
plot(data,'bx','linewidth',2);
hold off
title('Overlaid Demodulated Signal Constellation -- After Equalizer')
grid on
axis([-1.5,1.5,-1.5,1.5])
axis('square')

subplot(212)
plot(0,0)
hold on
for i = N_symbol-9:N_symbol
    plot([-100+4:100-4],real(temp(i,29:221)),'ro')
end
hold off
title('Spectrum Lines -- In Phase -- After Equalizer')
axis([-100 100 -2 2]);
grid on
xlabel('Frequency Index'); ylabel('Amplitude')

figure
subplot(211)
plot(0,0)
hold on
for n = N_symbol-9:N_symbol
    plot(ofdm_equ_data2(n,:),'ro')
end
plot(data,'bx','linewidth',2);
hold off
axis([-1.5,1.5,-1.5,1.5])
axis('square')
grid on
if chan_flg == 0
title('Overlaid Demodulated Signal Constellation -- E_b/N_o = 12dB')
else
title('Overlaid Demodulated Signal Constellation -- E_b/N_o = 14dB')
end
subplot(212)
plot(0,0)
hold on
for i = N_symbol-9:N_symbol
plot([-100+4:100-4],real(ofdm_equ_data2(i,:)),'ro')
end
hold off
axis([-100,100,-2,2]);
grid on
xlabel('Frequency Index'); ylabel('Amplitude')
title('Spectrum Lines -- In Phase')

figure('Color','w');
semilogy(v_EbNo_dB,v_ber,'x-','linewidth',2);
if chan_flg == 0
    title('16-QAM OFDM in AWGN channel')
else
    if chan_type == 1
        title('16-QAM OFDM in multipath fading channel')
    elseif chan_type == 2
        title('16-QAM OFDM in SUI-3 Channel')
    end
end
legend('Proposed Method')
set(legend,'Interpreter','latex','Location','SouthWest');
axis([min(v_EbNo_dB),max(v_EbNo_dB),min(v_ber) 1]); grid on;
xlabel('Eb/No (dB)','Interpreter','latex'); ylabel('Bit Error Rate','Interpreter','latex');
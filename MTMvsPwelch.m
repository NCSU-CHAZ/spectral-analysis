% figure comparing pwelch and multitaper techniques for thesis
% K. Anarde, 08/11/2018

% pick a time series segment (68 min)
% cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/WAVE_STATS/ADV/Bispectra/ADVsp-68min/
% load('IN-mtmBisp_ADVsp_68min.mat')
% nFs = 16;
% n = 8;

% pick a time series segment (17 min)
cd ~/Research/BARis/COASTRR/HARVEY/DATA/POSTSTORM/WAVE_STATS/ADV/Bispectra/ADVsp-17min/
load('IN-mtmBspd_ADVsp_17minZP.mat')
nFs = 16;
n   = 40; % 35

% choose a ~68 minute window
rData = detrend(ts(:,n));
nLengthT = length(rData);

% segment into 4 blocks
nSeg = 4;
nLengthS = nLengthT/nSeg;

% create Hamming window
rHam = hamming(nLengthT/nSeg);

% multiply each block by the taper
window1 = rHam .* rData(1:nLengthS); 
window2 = rHam .* rData(nLengthS+1:nLengthS*2); 
window3 = rHam .* rData((nLengthS*2)+1:nLengthS*3); 
window4 = rHam .* rData((nLengthS*3)+1:nLengthS*4); 

% calculate the one-sdied psd
psd1 = abs(fft(window1)/nLengthS); psd1 = psd1(1:nLengthS/2+1); psd1(2:end-1) = 2*psd1(2:end-1);
psd2 = abs(fft(window2)/nLengthS); psd2 = psd2(1:nLengthS/2+1); psd2(2:end-1) = 2*psd2(2:end-1);
psd3 = abs(fft(window3)/nLengthS); psd3 = psd3(1:nLengthS/2+1); psd3(2:end-1) = 2*psd3(2:end-1);
psd4 = abs(fft(window4)/nLengthS); psd4 = psd4(1:nLengthS/2+1); psd4(2:end-1) = 2*psd4(2:end-1);
f1 = nFs * (0:(nLengthS/2))/nLengthS;

% plot the original time series
figure; 
subplot(4,4,[1 4])
plot(rData); plotFancyAxis
xlim([0 16384])
% plot the hamming windows
subplot(4,4,5)
plot(rHam); plotFancyAxis
xlim([0 4096])
subplot(4,4,6)
plot(rHam); plotFancyAxis
xlim([0 4096])
subplot(4,4,7)
plot(rHam); plotFancyAxis
xlim([0 4096])
subplot(4,4,8)
plot(rHam); plotFancyAxis
xlim([0 4096])
% plot the tapered time series
subplot(4,4,9)
plot(window1); plotFancyAxis
xlim([0 4096])
subplot(4,4,10)
plot(window2); plotFancyAxis
xlim([0 4096])
subplot(4,4,11)
plot(window3); plotFancyAxis
xlim([0 4096])
subplot(4,4,12)
plot(window4); plotFancyAxis
xlim([0 4096])
% plot the one-sided psd
subplot(4,4,13)
semilogy(f1,psd1); plotFancyAxis
xlim([0 0.04])
subplot(4,4,14)
semilogy(f1,psd2); plotFancyAxis
xlim([0 0.04])
subplot(4,4,15)
semilogy(f1,psd3); plotFancyAxis
xlim([0 0.04])
subplot(4,4,16)
semilogy(f1,psd4); plotFancyAxis
xlim([0 0.04])

% plot the mean psd (i.e. a crude pwelch)
figure
% semilogy(f1,mean([psd1 psd2 psd3 psd4],2))
% hold on
% plot pwelch with ~8 overlapped segments
[psdW, fW] = pwelch(rData, [], [], [], nFs); 
semilogy(fW,psdW)
hold on
% plot snn from mtm in R
semilogy(f(:,n),snn(:,n))
xlim([0 0.04])
plotFancyAxis

% slepian tapers
NW = 3;
K  = 4;
dps = dpss(nLengthT,NW,K);

% taper the time series with the slepian tapers
window1 = dps(:,1) .* rData;
window2 = dps(:,2) .* rData;
window3 = dps(:,3) .* rData;
window4 = dps(:,4) .* rData;

% plot the original time series
figure; 
subplot(4,4,[1 4])
plot(rData); plotFancyAxis
xlim([0 16384])
% plot the slepian tapers
subplot(4,4,[5 8])
plot(dps); plotFancyAxis
xlim([0 16384])
% plot the tapered time series
subplot(4,4,[9 12])
plot(window1); plotFancyAxis
hold on
plot(window2); plotFancyAxis
hold on
plot(window3); plotFancyAxis
hold on
plot(window4); plotFancyAxis
xlim([0 16384])

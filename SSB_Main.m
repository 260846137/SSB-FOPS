clear all;
close all;
clc;

MaxCur = 1.5; %m/s;%%%%Use definded parameter--the maximum current velocity
load('TestData.mat');
figure;
imagesc(SpecHead.doppler_freq,SpecHead.rang,20*log10(RangeDopplerSpectrum));axis xy;hold on
xlabel('Doppler frequency (Hz)');ylabel('Range (km)')
colorbar;colormap(jet);
set(gcf,'color','w');

[F1] = FirstOrderPeaksDefinition_BasedOnSecondPeaks(RangeDopplerSpectrum,SpecHead,MaxCur);

plot(SpecHead.doppler_freq(F1.PosiL(3:end)),SpecHead.rang(3:40),'k','linewidth',1)
plot(SpecHead.doppler_freq(F1.PosiR(3:end)),SpecHead.rang(3:40),'k','linewidth',1)
plot(SpecHead.doppler_freq(F1.NegaL(3:end)),SpecHead.rang(3:40),'k','linewidth',1)
plot(SpecHead.doppler_freq(F1.NegaR(3:end)),SpecHead.rang(3:40),'k','linewidth',1)



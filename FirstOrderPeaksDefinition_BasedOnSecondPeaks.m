function [FirstOrderPeaks] = FirstOrderPeaksDefinition_BasedOnSecondPeaks(RangeDopplerSpectrum,SpecHead,MaxCur)
NoiseF = [2.7,3.2];
noisefact = 6.3;%%SNR threshold
NoiseF = NoiseF*SpecHead.FBragg;
temp = NoiseF(1) - SpecHead.doppler_freq;
[~,I] = min(abs(temp));
NoiseBin(1) = I;
temp = NoiseF(2) - SpecHead.doppler_freq;
[~,I] = min(abs(temp));
NoiseBin(2) = I;
temp = -NoiseF(2) - SpecHead.doppler_freq;
[~,I] = min(abs(temp));
NoiseBin(3) = I;
temp = -NoiseF(1) - SpecHead.doppler_freq;
[~,I] = min(abs(temp));
NoiseBin(4) = I;

for irange = 3:40
   
    Spectrum = RangeDopplerSpectrum(irange,:);
    Spectrum = Spectrum(:).^2;
    Spectrum = smooth(Spectrum,3);
    NoiseLevel(irange) = sum([Spectrum(NoiseBin(1):NoiseBin(2));Spectrum(NoiseBin(3):NoiseBin(4))])/...
        length([[Spectrum(NoiseBin(1):NoiseBin(2));Spectrum(NoiseBin(3):NoiseBin(4))]]);
    
    % step 1: extracting Doppler bins with frequency shift no more than
    % maximum  current velocity
    currmaxBin = round(MaxCur/median(diff(SpecHead.doppler_vel)));
    temp = Spectrum((SpecHead.iFBragg(2)-currmaxBin):(SpecHead.iFBragg(2)+currmaxBin));
    % step 2:找最强峰
    [BraggMaxPower,BraggMaxPowerBin] = max(abs(temp));
    BraggMaxPowerBin = SpecHead.iFBragg(2)-currmaxBin + BraggMaxPowerBin -1;
    % step 3: 计算阈值--正一阶峰
    BraggMaxPowerBinToSecondOrderPeaks = round(SpecHead.FBragg*sqrt(2)/ median(diff(SpecHead.doppler_freq))) ...
        + SpecHead.spectra_length/2;
    BraggMaxPowerBinToSecondOrderPeaks = (BraggMaxPowerBin - SpecHead.iFBragg(2)) +...
        BraggMaxPowerBinToSecondOrderPeaks;
    T = sum(Spectrum(BraggMaxPowerBinToSecondOrderPeaks-3:...
        BraggMaxPowerBinToSecondOrderPeaks+3))/(7);
    T = max([T,NoiseLevel(irange)*noisefact]);
    % 正一阶峰右边界
    temp = Spectrum(BraggMaxPowerBin:SpecHead.iFBragg(2)+currmaxBin);
    temp = temp - T;
    FdownPowerBinR = find(temp<0,1,'first');
     if isempty(FdownPowerBinR)
        PositiveRight(irange) = currmaxBin + SpecHead.iFBragg(2);
     else
        PositiveRight(irange) = BraggMaxPowerBin + FdownPowerBinR -2;
     end
     % 正一阶峰左边界
    temp = Spectrum(SpecHead.iFBragg(2)-currmaxBin:BraggMaxPowerBin);
    temp = flip(temp,1);
    temp = temp - T;
    FdownPowerBinR = find(temp<0,1,'first');
     if isempty(FdownPowerBinR)
        PositiveLeft(irange) = SpecHead.iFBragg(2) - currmaxBin;
     else
        PositiveLeft(irange) = BraggMaxPowerBin - FdownPowerBinR +2;
     end
     % 负一阶峰
     temp = Spectrum((SpecHead.iFBragg(1)-currmaxBin):(SpecHead.iFBragg(1)+currmaxBin));
     %找最强峰
    [BraggMaxPower,BraggMaxPowerBin] = max(temp);
    BraggMaxPowerBin = SpecHead.iFBragg(1)-currmaxBin + BraggMaxPowerBin -1;
     BraggMaxPowerBinToSecondOrderPeaks = SpecHead.spectra_length/2 - ...
         round(SpecHead.FBragg*sqrt(2)/ median(diff(SpecHead.doppler_freq)));
    BraggMaxPowerBinToSecondOrderPeaks = (BraggMaxPowerBin - SpecHead.iFBragg(1)) +...
        BraggMaxPowerBinToSecondOrderPeaks;
    T = sum(Spectrum(BraggMaxPowerBinToSecondOrderPeaks-3:...
        BraggMaxPowerBinToSecondOrderPeaks+3))/(7);
%     T = T*noisefact;
    T = max([T,NoiseLevel(irange)*noisefact]);
    % 负一阶峰右边界
    temp = Spectrum(BraggMaxPowerBin:SpecHead.iFBragg(1)+currmaxBin);
    temp = temp - T;
    FdownPowerBinR = find(temp<0,1,'first');
     if isempty(FdownPowerBinR)
        NegativeRight(irange) = currmaxBin + SpecHead.iFBragg(1);
     else
        NegativeRight(irange) = BraggMaxPowerBin + FdownPowerBinR -2;
     end
     % 负一阶峰左边界
    temp = Spectrum(SpecHead.iFBragg(1)-currmaxBin:BraggMaxPowerBin);
    temp = flip(temp,1);
    temp = temp - T;
    FdownPowerBinR = find(temp<0,1,'first');
     if isempty(FdownPowerBinR)
        NegativeLeft(irange) = SpecHead.iFBragg(1) - currmaxBin;
     else
        NegativeLeft(irange) = BraggMaxPowerBin - FdownPowerBinR +2;
     end
     
     %频点
     FirstOrderBins{irange} = [PositiveLeft(irange):PositiveRight(irange),...
         NegativeLeft(irange):NegativeRight(irange)];
     NumFirstOrderBin(irange) = length(FirstOrderBins{irange});
end
FirstOrderPeaks.PosiL = PositiveLeft;
FirstOrderPeaks.PosiR = PositiveRight;
FirstOrderPeaks.NegaL = NegativeLeft;
FirstOrderPeaks.NegaR = NegativeRight;
FirstOrderPeaks.index = FirstOrderBins;
FirstOrderPeaks.Num = NumFirstOrderBin;


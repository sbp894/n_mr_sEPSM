function [NormPowerMod] = MODENERGY(ModFreq, Freq, Amp)
% Created By: Varsha Hariram Rallapalli
% Date: 07/01/2015
% This function sums energy within a band with a specified modulation
% center frequency
% ModFreq = modulation center frequency
% Freq = Frequencies of the signal
% Amp = Energy at each frequency of the signal
switch(ModFreq)
    case 1
        LB_ModFreq = 0;
    otherwise
    LB_ModFreq = ModFreq./sqrt(2);
end

UB_ModFreq = ModFreq.*sqrt(2);

Lower = find(Freq>LB_ModFreq,1,'first');
Upper = find(Freq<UB_ModFreq,1,'last');

PowerMod = sum(Amp(Lower:Upper));

% NormPowerMod = PowerMod/mean(Amp(1:end)); % Normalize using the mean signal amp

%% MH: This value (Amp(0) seems to be near the noise floor, near 0, since we subtract one from the SUMCOR before taking FFT.
% So, I think we have already normalized in computing the SUMCOR - to make
% 1 the no-correlation value, so let's try it withouth this normalization.
% The PSDs seems to make sense, so then diving by a number that is near 0
% and in the noise floor seems like a bad idea.  we can look into this in
% more detail later on, but right now it is obscuring our results i think.

% NormPowerMod = PowerMod/Amp(1); % Normalize using DC power
NormPowerMod = PowerMod;

end


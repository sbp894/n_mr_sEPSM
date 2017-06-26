function stim = mnoise(cf,dur,ramplength,silencelength,fs)
% function y=gnoise(lf,hf,duration,dBNo,fs) makes a white noise.
%
% lf            Low frequency
% hf            High frequency
% dur           stimulus duration in seconds
% dBNo          Noise spectrum level (not dB SPL)
% ramplength    Onset and offset ramp duration in seconds
% silencelength Duration of silence in seconds (pre and post signal)
% fs            Sampling frequency
%

%% Make signal shape
tdres = 1/fs;
rampts = ramplength * fs;
steadypts = ceil(dur * fs - 2*rampts);
totalpts = steadypts + (rampts*2);
step = pi/(rampts-1);
x=[0:step:pi];
offramp = (1+cos(x))./2;
onramp = (1+cos(fliplr(x)))./2;
o=ones(1,(steadypts));
wholeramp = [onramp o offramp]; % Envelope for stimulus (i.e. on/off ramps)
zerolength = ceil(silencelength * fs);
%% Calculate scale
q = 2.871; % for 1 octave bandwidth 
bw =  cf/q;
lf = cf * (sqrt(1+(1/(4*q^2))) - 1/(2*q));
hf = cf * (sqrt(1+(1/(4*q^2))) + 1/(2*q));
%% Make noise
% Make 10 noise instances to calculate mean RMS, we will use the last noise
% instance
fres = 1./(totalpts * tdres);
randn('state',sum(100*clock)); % seed the random # generator with the clock
temp_rms = [];
for i = 1:10
    stim = randn(1,totalpts);
    stimspectrum = fft(stim);
    nptslofreq = floor(lf/fres);
    nptshifreq = ceil(hf/fres);
    stimspectrum(1:nptslofreq) = 0.;  % zero out all bands that apply (Mickey's bug-fix in here)
    stimspectrum(nptshifreq:(totalpts - nptshifreq+2)) = 0.;
    stimspectrum((totalpts - nptslofreq+2):totalpts) = 0.;
    % So, values in the fft of stim1 (both real & imag parts) have been zeroed
    stim = real(ifft(stimspectrum)); % inverse FFT to create noise waveform
    temp_rms = [temp_rms, rms(stim)];
end
tmp_rms = mean(temp_rms);
stim = stim/tmp_rms .* wholeramp; %normalize by rms and scale up to desired rms
stim = [zeros(zerolength,1); stim'; zeros(zerolength,1)];
return
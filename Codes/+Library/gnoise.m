function stim = gnoise(lf,hf,dur,dBNo,ramplength,silencelength,fs)
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
bw = round(hf-lf);
if bw == 0
    bw = 1;
end
rmsnoise = dBNo + (10 .* log10(bw));
scale_noise = 20e-6*10^((rmsnoise)/20);  % into Pascals

%% Make noise
% Make 10 noise instances to calculate mean RMS, we will use the last noise
% instance
fres = 1./(totalpts * tdres);
randn('state',sum(100*clock)); % seed the random # generator with the clock
temp_rms = [];
for ii = 1:10
    stim = randn(1,totalpts);
    stimspectrum = fft(stim);
    nptslofreq = lf/fres; %
    nptshifreq = hf/fres;  %
    idx_hf = nptshifreq+1;
    idx_lf = nptslofreq+1;
    spec_pos = stimspectrum(1:((fs/2)/fres)+1);  % positive freqs inclusing 0 and fs/2
    spec_pos(1:idx_lf-1) = 0;
    spec_pos(idx_hf+1:end) = 0;
    % Define phases
    % spec_pos(idx_lf:idx_hf) = spec_pos(idx_lf:idx_hf) * exp(1i + 0);
    spec2ifft = [spec_pos conj(fliplr(spec_pos(2:end-1)))];
    stim = ifft(spec2ifft);
    temp_rms = [temp_rms, rms(stim)];
end
tmp_rms = mean(temp_rms);
stim = stim/tmp_rms * scale_noise .* wholeramp; %normalize by rms and scale up to desired rms
stim = [zeros(zerolength,1); stim'; zeros(zerolength,1)];
return
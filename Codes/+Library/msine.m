function stim = msine(mf,dur,phase,ramplength,silencelength,fs)
% function stim = msine(mf,dur,ramplength,silencelength,fs) makes a sine
% wave modulator.
%
% mf            Modulation frequency
% dur           Stimulus duration in seconds
% phase         Starting phase
% ramplength    Onset and offset ramp duration in seconds
% silencelength Duration of silence in seconds (pre and post signal)
% fs            Sampling frequency

%% Make ramps and silence
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

%% Make modulator
nsamples=round(dur*fs);
stim = sin(phase+2*pi*mf.*(1:nsamples)'./fs);

stim = stim.* wholeramp(1:length(stim))'; % Apply ramps
stim = [zeros(zerolength,1); stim; zeros(zerolength,1)]; % Apply silence
return
function [ y ] = getSSN( speech, snr, fs )
%GETSSN load noise file and adjust level to fit signal to noise ratio
    [x,Fs] = audioread(['+Stim/ssn.wav']);
    y = resample(x,fs,Fs);
    nSpeech = length(speech);
    nSeg = floor(length(y)/nSpeech);
    % pick a random segment from the noise file
    startIdx = randi(nSeg-2 ,1)*nSpeech;
    y = y(startIdx:startIdx+nSpeech -1)';
    srms=Library.rms(speech);
    nrms=Library.rms(y);
    y=y.*(srms/(10^(snr/10)))/nrms;
    % check snr
    % snr_ck = 10*log10(Library.rms(speech)/Library.rms(y));
    y=y';
end


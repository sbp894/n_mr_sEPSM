function [ y, fs ] = getWav( fileName, fs )
%GETWAV loads a wav file and resamples it, if necessary

    [x,Fs] = audioread(['+Stim/' fileName]);
    y = resample(x,fs,Fs);

end


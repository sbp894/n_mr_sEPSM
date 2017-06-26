function [ y, fs ] = getWavFromPath( filePath, fileName, fs )
%GETWAV loads a wav file and resamples it, if necessary

    [x,Fs] = audioread([filePath fileName]);
    y = resample(x,fs,Fs);

end


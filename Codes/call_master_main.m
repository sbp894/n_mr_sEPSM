%%
clear;
close all;
clc;

global paramsIfCallingFromOutside

rng('default');
fiberTypes=1:3;
noiseTypes={{'SSN'},{'SAM'}};

%%
for fibVar=1:length(fiberTypes)
    paramsIfCallingFromOutside.fiberType=fiberTypes(fibVar); % L/M/H <--> 1/2/3
    for noiseVar=2%1:length(noiseTypes)
        paramsIfCallingFromOutside.noiseTypes=noiseTypes{noiseVar};
        paramsIfCallingFromOutside.SNR=-21:3:12;
        paramsIfCallingFromOutside.level=80;
        paramsIfCallingFromOutside.nRep=25;
        %     paramsIfCallingFromOutside.noiseTypes={'SSN'};
        paramsIfCallingFromOutside.CF=logspace(log10(125), log10(8e3), 21);
        paramsIfCallingFromOutside.species=2;    % 1 for cat (2 for human with Shera et al. tuning; 3 for human with Glasberg & Moore tuning)
        
        master_main;
    end
end

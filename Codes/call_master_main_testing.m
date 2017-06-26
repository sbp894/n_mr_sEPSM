%%
clear;
close all;
clc;

global paramsIfCallingFromOutside

rng('default');
fiberTypes=2;
noiseTypes={{'SAM'}};

%%
for fibVar=1:length(fiberTypes)
    paramsIfCallingFromOutside.fiberType=fiberTypes(fibVar); % L/M/H <--> 1/2/3
    for noiseVar=1:length(noiseTypes)
        paramsIfCallingFromOutside.noiseTypes=noiseTypes{noiseVar};
        paramsIfCallingFromOutside.SNR=-3:6:3;
        paramsIfCallingFromOutside.level=65;
        paramsIfCallingFromOutside.nRep=25;
        %     paramsIfCallingFromOutside.noiseTypes={'SSN'};
        paramsIfCallingFromOutside.CF=1e3; %logspace(log10(125), log10(8e3), 21);
        paramsIfCallingFromOutside.species=2;    % 1 for cat (2 for human with Shera et al. tuning; 3 for human with Glasberg & Moore tuning)
        
        master_main;
    end
end

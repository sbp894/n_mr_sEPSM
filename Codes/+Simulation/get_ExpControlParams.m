function ExpControlParams=get_ExpControlParams()

ExpControlParams.SNR=-21:12:12;
ExpControlParams.level=65;


ExpControlParams.noiseTypes={'SAM', 'SSN'};

ExpControlParams.noisePrefix = cell(1, length(ExpControlParams.noiseTypes));
for i=1:length(ExpControlParams.noiseTypes)
    ExpControlParams.noisePrefix{i} = ['noise' filesep lower(ExpControlParams.noiseTypes{i}) '_simulation_dtu'];
end

ExpControlParams.fiberType=2:3; %1:3; % L/M/H <--> 1/2/3
ExpControlParams.CF=[.5 2]*1e3; %logspace(log10(125), log10(8e3), 21);
ExpControlParams.species=2;    % 1 for cat (2 for human with Shera et al. tuning; 3 for human with Glasberg & Moore tuning)
ExpControlParams.sentences=3:4; %1:10;
ExpControlParams.ohcLoss_dB=0;
ExpControlParams.ihcLoss_dB=0;

ExpControlParams.nRep=25;
ExpControlParams.BootstrapLoopMax=24;
ExpControlParams.BootstrapLoopReport=60;
ExpControlParams.nPSDs2Avg=12;
ExpControlParams.fixSPL=0;
ExpControlParams.winCorr0Add1Mul=1;
ExpControlParams.fs=100e3; % Sampling Frequency
ExpControlParams.modFreqs=[1 2 4 8 16 32 64];
ExpControlParams.mrWindows=2.^10./ExpControlParams.modFreqs*1e-3;
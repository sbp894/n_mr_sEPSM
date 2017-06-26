function [PSDenv_STRUCT,PSDtfs_STRUCT,PowerMod_STRUCT,PowerTfs_STRUCT] = sumcors_sp(SpikeTrains,paramsIN, resultsDir,resultPostfix)

global figHandles

% This function returns the FFT of the SUMCORs (envelope) for each stimulus derived from the
% spike output of one AN fiber. It also returns the envelope power within a
% certain number of modulation frequency bands.

%%% Input parmeters
% stim_SN: Sentence plus noise
% stim_N: Noise
% CF_kHz: Characteristic frequency of the nerve fiber
% fiberType: 1-LSR, 2-MSR, 3-HSR (Spontaneous rate)
% Cohc: Outer hair cell parameter
% Cihc: Inner hair cell parameter
% OALevel_dBSPL: Stimulus intensity (for both conditions)


%% *MH: HOW TO OPTIMIZE??
BootstrapLoopMax=12;
[SACSCCfunctions,SACSCCmetrics,paramsOUT] = SACSCCanal_SNRenv_parallel(SpikeTrains,paramsIN,1,resultsDir,resultPostfix,BootstrapLoopMax);

Library.parsave([resultsDir 'sacMET' filesep 'sacMET' resultPostfix '.mat'],  SACSCCmetrics );

%%% PSDenv computed in SACSCCanal_SNRenv and passed within SACSCCfunctions
% SACSCCfunctions{end} gives the all the average values
PSDenv_S=SACSCCfunctions{end}.PSDenv_A;
PSDenv_N=SACSCCfunctions{end}.PSDenv_B;
PSDenv_SN=SACSCCfunctions{end}.PSDenv_C;
PSDenv_S_noisefloor=SACSCCfunctions{end}.rand.PSDenv_A;
PSDenv_N_noisefloor=SACSCCfunctions{end}.rand.PSDenv_B;
PSDenv_SN_noisefloor=SACSCCfunctions{end}.rand.PSDenv_C;
PSDfreqVEC_Hz=SACSCCfunctions{end}.PSD_freqVEC;

% tfs
tfsCutOff = 1.5;
CF_A_kHz=paramsIN.CF_A_Hz/1e3;

CF_B_kHz=CF_A_kHz;
CF_C_kHz=CF_A_kHz;

if CF_A_kHz < tfsCutOff && CF_B_kHz < tfsCutOff && CF_C_kHz < tfsCutOff
    PSDtfs_S=SACSCCfunctions{end}.PSDtfs_A;
    PSDtfs_N=SACSCCfunctions{end}.PSDtfs_B;
    PSDtfs_SN=SACSCCfunctions{end}.PSDtfs_C;
    PSDtfs_S_noisefloor=SACSCCfunctions{end}.rand.PSDtfs_A;
    PSDtfs_N_noisefloor=SACSCCfunctions{end}.rand.PSDtfs_B;
    PSDtfs_SN_noisefloor=SACSCCfunctions{end}.rand.PSDtfs_C;
else
    PSDtfs_S=0;
    PSDtfs_N=0;
    PSDtfs_SN=0;
    PSDtfs_S_noisefloor=0;
    PSDtfs_N_noisefloor=0;
    PSDtfs_SN_noisefloor=0;
end

%% *MH - look at spectra
if paramsIN.plt
    figure(figHandles.PSDplot);
    clf;
    set (gcf, 'Units', 'normalized', 'Position', [.1 .1 .8 .8]);
    plot(PSDfreqVEC_Hz,PSDenv_S,'b-');
    hold on; plot(PSDfreqVEC_Hz,PSDenv_N,'k-');
    plot(PSDfreqVEC_Hz,PSDenv_SN,'r-');
    xlim([0 64])
    plot(PSDfreqVEC_Hz,PSDenv_S_noisefloor,'b:');
    plot(PSDfreqVEC_Hz,PSDenv_N_noisefloor,'k:');
    plot(PSDfreqVEC_Hz,PSDenv_SN_noisefloor,'r:');
    legend('S','N','SN','NF[S]','NF[N]','NF[SN]')
    xlabel(sprintf('MODULATION FREQUENCY (Hz)'),'FontSize',12,'HorizontalAlignment','center')
    ylabel(sprintf('POWER SPECTRAL DENSITY'),'FontSize',12)
    
    title(sprintf('CF=%.3f kHz; SRtype=%d (1: LSR; 2: MSR; 3: HSR);\nCohc=%.2f; Cihc=%.2f; Overall (Signal) Level: %.1f dB SPL; SNR = %.1f dB', ...
        CF_C_kHz,paramsIN.FiberType,99,99,paramsIN.level,paramsIN.SNR2use_dB),'FontSize',12);
    
    Library.saveFigureAs([resultsDir 'eps' filesep 'eps' resultPostfix '.eps'])
    Library.saveFigureAs([resultsDir 'png' filesep 'psd' resultPostfix '.tiff'])
end
% Divide the FFT output into modulation filters
ModFreq = [1,2,4,8,16,32,64];
PowerModS = zeros(1,length(ModFreq));
PowerModSN = zeros(1,length(ModFreq));
PowerModN = zeros(1,length(ModFreq));
PowerModS_noisefloor = zeros(1,length(ModFreq));
PowerModSN_noisefloor = zeros(1,length(ModFreq));
PowerModN_noisefloor = zeros(1,length(ModFreq));

%% Look at RAND PSDs as well to get noise floor.
for i = 1:length(ModFreq)
    PowerModS(i) = Library.MODENERGY(ModFreq(i), PSDfreqVEC_Hz, PSDenv_S);
    PowerModN(i) = Library.MODENERGY(ModFreq(i), PSDfreqVEC_Hz, PSDenv_N);
    PowerModSN(i) = Library.MODENERGY(ModFreq(i), PSDfreqVEC_Hz, PSDenv_SN);
    PowerModS_noisefloor(i) = Library.MODENERGY(ModFreq(i), PSDfreqVEC_Hz, PSDenv_S_noisefloor);
    PowerModN_noisefloor(i) = Library.MODENERGY(ModFreq(i), PSDfreqVEC_Hz, PSDenv_N_noisefloor);
    PowerModSN_noisefloor(i) = Library.MODENERGY(ModFreq(i), PSDfreqVEC_Hz, PSDenv_SN_noisefloor);
end

%Pass back structures to make cleaner (to fit in noise floor comps)
PSDenv_STRUCT.PSDenv_S=PSDenv_S;
PSDenv_STRUCT.PSDenv_N=PSDenv_N;
PSDenv_STRUCT.PSDenv_SN=PSDenv_SN;
PSDenv_STRUCT.PSDfreqVEC_Hz=PSDfreqVEC_Hz;
PSDenv_STRUCT.PSDenv_S_noisefloor=PSDenv_S_noisefloor;
PSDenv_STRUCT.PSDenv_N_noisefloor=PSDenv_N_noisefloor;
PSDenv_STRUCT.PSDenv_SN_noisefloor=PSDenv_SN_noisefloor;

PSDtfs_STRUCT.PSDtfs_S=PSDtfs_S;
PSDtfs_STRUCT.PSDtfs_N=PSDtfs_N;
PSDtfs_STRUCT.PSDtfs_SN=PSDtfs_SN;
PSDtfs_STRUCT.PSDfreqVEC_Hz=PSDfreqVEC_Hz;
PSDtfs_STRUCT.PSDtfs_S_noisefloor=PSDtfs_S_noisefloor;
PSDtfs_STRUCT.PSDtfs_N_noisefloor=PSDtfs_N_noisefloor;
PSDtfs_STRUCT.PSDtfs_SN_noisefloor=PSDtfs_SN_noisefloor;

if CF_A_kHz < tfsCutOff && CF_B_kHz < tfsCutOff && CF_C_kHz < tfsCutOff
    NYQ=0.5*(1/paramsOUT.DELAYbinwidth_sec);
    PowerTfs_STRUCT.PowerTfsS = Library.MODENERGY(NYQ, PSDfreqVEC_Hz, PSDtfs_S);
    PowerTfs_STRUCT.PowerTfsN = Library.MODENERGY(NYQ, PSDfreqVEC_Hz, PSDtfs_N);
    PowerTfs_STRUCT.PowerTfsSN = Library.MODENERGY(NYQ, PSDfreqVEC_Hz, PSDtfs_SN);
    PowerTfs_STRUCT.PowerTfsS_noisefloor = Library.MODENERGY(NYQ, PSDfreqVEC_Hz, PSDtfs_S_noisefloor);
    PowerTfs_STRUCT.PowerTfsN_noisefloor = Library.MODENERGY(NYQ, PSDfreqVEC_Hz, PSDtfs_N_noisefloor);
    PowerTfs_STRUCT.PowerTfsSN_noisefloor = Library.MODENERGY(NYQ, PSDfreqVEC_Hz, PSDtfs_SN_noisefloor);
else
    PowerTfs_STRUCT.PowerTfsS = 0;
    PowerTfs_STRUCT.PowerTfsN = 0;
    PowerTfs_STRUCT.PowerTfsSN = 0;
    PowerTfs_STRUCT.PowerTfsS_noisefloor = 0;
    PowerTfs_STRUCT.PowerTfsN_noisefloor = 0;
    PowerTfs_STRUCT.PowerTfsSN_noisefloor = 0;
end

PowerMod_STRUCT.PowerModS=PowerModS;
PowerMod_STRUCT.PowerModN=PowerModN;
PowerMod_STRUCT.PowerModSN=PowerModSN;
PowerMod_STRUCT.PowerModS_noisefloor=PowerModS_noisefloor;
PowerMod_STRUCT.PowerModN_noisefloor=PowerModN_noisefloor;
PowerMod_STRUCT.PowerModSN_noisefloor=PowerModSN_noisefloor;

end
function [PSDenv_STRUCT,PSDtfs_STRUCT,PowerMod_STRUCT,PowerTfs_STRUCT] = sumcors_mod(stim_S, stim_N, stim_SN, A, B, C, AN, anal)

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SET PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CF_kHz = AN.CF(C.cF_i)/1e3;
fiberType = AN.fiberType.sr;
Cohc = AN.Cohc(C.cF_i);
Cihc = AN.Cihc(C.cF_i);
OALevel_dBSPL = A.level(C.level_i);
SNR2use_dB = B.SNR(C.snr_i);
resultDir = anal.resultsDir;
resultPostfix1 = sprintf(anal.resultTxt,AN.CF(C.cF_i)/1e3,C.sentence_i,B.SNR(C.snr_i),B.noiseTypes{C.noise_i},A.level(C.level_i) );
resultPostfix2 = sprintf(anal.resultFileName,AN.CF(C.cF_i)/1e3,C.sentence_i,B.SNR(C.snr_i), B.noiseTypes{C.noise_i}, A.level(C.level_i) );
resultPostfix{1} = resultPostfix1(~isspace(resultPostfix1));
resultPostfix{2} = resultPostfix2(~isspace(resultPostfix2));
plt=anal.plt;
verbose=anal.verbose;

%% MODEL params - general (for both conditions)
ANmodel_Fs_Hz=100000;
Nreps=100;

%% Condition parameters
%%%%%%%%%%%%%
% Condition A
% stim_FileName_A=stim_S;
CF_A_kHz=CF_kHz;
SRtype_A=fiberType; % fiberType: 1 - LSR, 2 - MSR, 3 - HSR
Cohc_A=Cohc;
Cihc_A=Cihc;

%%%%%%%%%%%%%
% Condition B
% stim_FileName_B=stim_N;
CF_B_kHz=CF_kHz;
SRtype_B=fiberType;
Cohc_B=Cohc;
Cihc_B=Cihc;

%%%%%%%%%%%%%
% Condition C
% stim_FileName_C=stim_SN;
CF_C_kHz=CF_kHz;
SRtype_C=fiberType;
Cohc_C=Cohc;
Cihc_C=Cihc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stimuli generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp(stim_FileName_A)
% [stim_A_model FsA_Hz] = audioread(stim_FileName_A);
% [stim_B_model FsB_Hz] = audioread(stim_FileName_B);
% [stim_C_model FsC_Hz] = audioread(stim_FileName_C);
stim_A_model = stim_S;   
stim_B_model= stim_N;
stim_C_model = stim_SN;
% if (FsA_Hz~=ANmodel_Fs_Hz)|(FsB_Hz~=ANmodel_Fs_Hz)|(FsC_Hz~=ANmodel_Fs_Hz)
%     error('WAV file (S, N, or SN) does not have correct ANmodel Sampling Rate Fs')  % Code assumes AN model Sampling Rate is already set in WAV files
% end
% clear FsA_Hz FsB_Hz FsC_Hz
dur_sec=min([length(stim_A_model)/ANmodel_Fs_Hz length(stim_B_model)/ANmodel_Fs_Hz length(stim_C_model)/ANmodel_Fs_Hz]);  % Set duration to minimum of all stimuli

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OALevel_dBSPL = BML_dBSPL;  % set OAL to BML for Speech IQ for this fiber/condition


dBSPL_A_after=20*log10(sqrt(mean(stim_A_model.^2))/(20e-6));
%  	dBSPL_B_after=20*log10(sqrt(mean(stim_B_model.^2))/(20e-6));
%  	dBSPL_C_after=20*log10(sqrt(mean(stim_C_model.^2))/(20e-6));

%% Scale stimuli to correct OALevel_dBSPL
% Note: Stim scaling for SNRenv study is done to set the SNR by
% adjusting the Noise level, so SIGNAL is same level in all WAV files,
% and SIGNAL should be scaled to provide desired OALlevel, and then all
% three WAV files scaled by the same scale factor.
stim_A_model=stim_A_model*10^((OALevel_dBSPL-dBSPL_A_after)/20);
stim_B_model=stim_B_model*10^((OALevel_dBSPL-dBSPL_A_after)/20);
stim_C_model=stim_C_model*10^((OALevel_dBSPL-dBSPL_A_after)/20);
%% REFIT and WINDOWwavefile at ANmodel_Fs
% Repeat or truncate waveform to fit requested stimulus duration:
stim_A_model = Library.refit_waveform(stim_A_model,ANmodel_Fs_Hz,dur_sec*1000);
stim_B_model = Library.refit_waveform(stim_B_model,ANmodel_Fs_Hz,dur_sec*1000);
stim_C_model = Library.refit_waveform(stim_C_model,ANmodel_Fs_Hz,dur_sec*1000);
% Window waveform using linear rise/fall:
stim_A_model = Library.window_waveform(stim_A_model,ANmodel_Fs_Hz,dur_sec*1000);
stim_B_model = Library.window_waveform(stim_B_model,ANmodel_Fs_Hz,dur_sec*1000);
stim_C_model = Library.window_waveform(stim_C_model,ANmodel_Fs_Hz,dur_sec*1000);

%% Listen to Stimuli at the Levels and SNR played to the model - listen at least once for each new condition to verify
PLAYstim = 0;  % Set to 0 to not hear stimuli
if PLAYstim
    PAUSEtime_sec=3;
    sound(stim_A_model, ANmodel_Fs_Hz)
    disp('... Playing 100kHz SIGNAL')
    pause(PAUSEtime_sec)
    sound(stim_B_model, ANmodel_Fs_Hz)
    disp('... Playing 100kHz NOISE')
    pause(PAUSEtime_sec)
    sound(stim_C_model, ANmodel_Fs_Hz)
    disp('... Playing 100kHz SIGNAL+NOISE')
    pause(PAUSEtime_sec)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Spike generation (from AN model for all 6 conditions: A+,A-,B+,B-,C+,C-)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp(sprintf('... Processing CONDITION:A \n    ... CF=%.3f kHz; SRtype=%d (1: LSR; 2: MSR; 3: HSR); Cohc=%.2f; Cihc=%.2f; Nreps=%d;\n    ... Stim_A="%s"\n    ... Overall (Signal) Level: %.1f dB SPL; SNR = %.1f dB; Duration = %.2f sec', ...
%    CF_A_kHz,SRtype_A,Cohc_A,Cihc_A,Nreps,stim_FileName_A,OALevel_dBSPL,SNR2use_dB,dur_sec))


%% % 2014 model
%%%%%%%%%%%%%%%%
%% STIM_A = Signal
%%%%%%%%%%%%%%%%
%% stim_A_plus
vIHC = ANModel.model_IHC(stim_A_model.',CF_A_kHz*1000,1,1/ANmodel_Fs_Hz,dur_sec+0.05,Cohc_A,Cihc_A,2); % 2013 model; "2" for human with BM tuning from Shera et al. (PNAS 2002)
if verbose
    [sout_Ap,~, ~] = ANModel.model_Synapse(vIHC,CF_A_kHz*1000,1,1/ANmodel_Fs_Hz,SRtype_A,1,1);
else
    evalc('[sout_Ap,~, ~] = ANModel.model_Synapse(vIHC,CF_A_kHz*1000,1,1/ANmodel_Fs_Hz,SRtype_A,1,1);');
end
[sptimes, ~]= Library.SGfast([1/ANmodel_Fs_Hz, Nreps],sout_Ap);
% Format spikes into NEL spikes format then cell array
NELspikes=Library.ANmodelSTs2nel(sptimes,Nreps);
SpikeTrainsA_plus=Library.nelSTs2cell(NELspikes);

%% stim_A_minus
vIHC = ANModel.model_IHC(-stim_A_model.',CF_A_kHz*1000,1,1/ANmodel_Fs_Hz,dur_sec+0.05,Cohc_A,Cihc_A,2); % 2013 model
if verbose
    [sout_An,~, ~] = ANModel.model_Synapse(vIHC,CF_A_kHz*1000,1,1/ANmodel_Fs_Hz,SRtype_A,0,0);
else
    evalc('[sout_An,~, ~] = ANModel.model_Synapse(vIHC,CF_A_kHz*1000,1,1/ANmodel_Fs_Hz,SRtype_A,0,0);');
end
[sptimes , ~]= Library.SGfast([1/ANmodel_Fs_Hz, Nreps],sout_An);
% Format spikes into NEL spikes format then cell array
NELspikes=Library.ANmodelSTs2nel(sptimes,Nreps);
SpikeTrainsA_minus=Library.nelSTs2cell(NELspikes);

%%%%%%%%%%%%%%%%
%% STIM_B = Noise
%%%%%%%%%%%%%%%%
%disp(sprintf('... Processing CONDITION:B \n    ... CF=%.3f kHz; SRtype=%d (1: LSR; 2: MSR; 3: HSR); Cohc=%.2f; Cihc=%.2f; Nreps=%d;\n    ... Stim_B="%s"\n    ... Overall (Signal) Level: %.1f dB SPL; SNR = %.1f dB; Duration = %.2f sec', ...
%    CF_B_kHz,SRtype_B,Cohc_B,Cihc_B,Nreps,stim_FileName_B,OALevel_dBSPL,SNR2use_dB,dur_sec))
%% stim_B_plus
vIHC = ANModel.model_IHC(stim_B_model.',CF_B_kHz*1000,1,1/ANmodel_Fs_Hz,dur_sec+0.05,Cohc_B,Cihc_B,2); % 2013 model
if verbose
    [sout_Bp,~, ~] = ANModel.model_Synapse(vIHC,CF_B_kHz*1000,1,1/ANmodel_Fs_Hz,SRtype_B,0,0);
else
    evalc('[sout_Bp,~, ~] = ANModel.model_Synapse(vIHC,CF_B_kHz*1000,1,1/ANmodel_Fs_Hz,SRtype_B,0,0);');
end

[sptimes , ~]= Library.SGfast([1/ANmodel_Fs_Hz, Nreps],sout_Bp);
% Format spikes into NEL spikes format then cell array
NELspikes=Library.ANmodelSTs2nel(sptimes,Nreps);
SpikeTrainsB_plus=Library.nelSTs2cell(NELspikes);

%% stim_B_minus
vIHC = ANModel.model_IHC(-stim_B_model.',CF_B_kHz*1000,1,1/ANmodel_Fs_Hz,dur_sec+0.05,Cohc_B,Cihc_B,2); % 2013 model
if verbose
    [sout_Bn,~, ~] = ANModel.model_Synapse(vIHC,CF_B_kHz*1000,1,1/ANmodel_Fs_Hz,SRtype_B,0,0);
else
    evalc('[sout_Bn,~, ~] = ANModel.model_Synapse(vIHC,CF_B_kHz*1000,1,1/ANmodel_Fs_Hz,SRtype_B,0,0);');
end


[sptimes , ~]= Library.SGfast([1/ANmodel_Fs_Hz, Nreps],sout_Bn);
% Format spikes into NEL spikes format then cell array
NELspikes=Library.ANmodelSTs2nel(sptimes,Nreps);
SpikeTrainsB_minus=Library.nelSTs2cell(NELspikes);

%%%%%%%%%%%%%%%%
%% STIM_C = Signal+Noise
%%%%%%%%%%%%%%%%
%disp(sprintf('... Processing CONDITION:C \n    ... CF=%.3f kHz; SRtype=%d (1: LSR; 2: MSR; 3: HSR); Cohc=%.2f; Cihc=%.2f; Nreps=%d;\n    ... Stim_C="%s"\n    ... Overall (Signal) Level: %.1f dB SPL; SNR = %.1f dB; Duration = %.2f sec', ...
%    CF_C_kHz,SRtype_C,Cohc_C,Cihc_C,Nreps,stim_FileName_C,OALevel_dBSPL,SNR2use_dB,dur_sec))
%% stim_C_plus
vIHC = ANModel.model_IHC(stim_C_model.',CF_C_kHz*1000,1,1/ANmodel_Fs_Hz,dur_sec+0.05,Cohc_C,Cihc_C,2); % 2013 model
if verbose
    [sout_Cp,~, ~] = ANModel.model_Synapse(vIHC,CF_C_kHz*1000,1,1/ANmodel_Fs_Hz,SRtype_C,0,0);
else
    evalc('[sout_Cp,~, ~] = ANModel.model_Synapse(vIHC,CF_C_kHz*1000,1,1/ANmodel_Fs_Hz,SRtype_C,0,0);');
end

[sptimes , ~]= Library.SGfast([1/ANmodel_Fs_Hz, Nreps],sout_Cp);
% Format spikes into NEL spikes format then cell array
NELspikes=Library.ANmodelSTs2nel(sptimes,Nreps);
SpikeTrainsC_plus=Library.nelSTs2cell(NELspikes);

%% stim_C_minus
vIHC = ANModel.model_IHC(-stim_C_model.',CF_C_kHz*1000,1,1/ANmodel_Fs_Hz,dur_sec+0.05,Cohc_C,Cihc_C,2); % 2013 model
if verbose
    [sout_Cn,~, ~] = ANModel.model_Synapse(vIHC,CF_C_kHz*1000,1,1/ANmodel_Fs_Hz,SRtype_C,0,0);
else
    evalc('[sout_Cn,~, ~] = ANModel.model_Synapse(vIHC,CF_C_kHz*1000,1,1/ANmodel_Fs_Hz,SRtype_C,0,0);');
end

[sptimes, ~]= Library.SGfast([1/ANmodel_Fs_Hz, Nreps],sout_Cn);
% Format spikes into NEL spikes format then cell array
NELspikes=Library.ANmodelSTs2nel(sptimes,Nreps);


SpikeTrainsC_minus=Library.nelSTs2cell(NELspikes);

if plt
    figure; clf; ymaxTEMP=-Inf;
    timeVEC1_sec=((1:length(stim_A_model))-1)/ANmodel_Fs_Hz;
    subplot(411)
    plot(timeVEC1_sec,stim_C_model,'b'); hold on; plot(timeVEC1_sec,stim_B_model,'r'); plot(timeVEC1_sec,stim_A_model,'g');
    legend('SN','N','S')
    title(sprintf('CF=%.3f kHz; SRtype=%d (1: LSR; 2: MSR; 3: HSR); Cohc=%.2f; Cihc=%.2f; Overall (Signal) Level: %.1f dB SPL; SNR = %.1f dB', ...
        CF_C_kHz,SRtype_C,Cohc_C,Cihc_C,OALevel_dBSPL,SNR2use_dB))
    ylabel(sprintf('Stimulus Amplitude\n(Pascals)'))
    subplot(412)
    timeVEC_sec=((1:length(sout_Ap))-1)/ANmodel_Fs_Hz;
    plot(timeVEC_sec,sout_Ap,'g'); hold on; plot(timeVEC_sec,sout_An,'g:')
    ylabel(sprintf('Synapse Output\n(spikes/sec)'))
    title('Signal')
    legend('POS','NEG')
    ymaxTEMP=max([max(ylim) ymaxTEMP]);
    subplot(413)
    plot(timeVEC_sec,sout_Bp,'r'); hold on; plot(timeVEC_sec,sout_Bn,'r:')
    ylabel(sprintf('Synapse Output\n(spikes/sec)'))
    title('Noise')
    legend('POS','NEG')
    ymaxTEMP=max([max(ylim) ymaxTEMP]);
    subplot(414)
    plot(timeVEC_sec,sout_Cp,'b'); hold on; plot(timeVEC_sec,sout_Cn,'b:')
    ylabel(sprintf('Synapse Output\n(spikes/sec)'))
    title('Signal + Noise')
    xlabel('Time (sec)')
    legend('POS','NEG')
    ymaxTEMP=max([max(ylim) ymaxTEMP]);
    subplot(412); ylim([0 ymaxTEMP]); subplot(413); ylim([0 ymaxTEMP]); subplot(414); ylim([0 ymaxTEMP]);
end
clear timeout meout c1filterout c2filterout c1vihc c2vihc vihc psth NELspikes sptimes
clear soutAp soutAn soutBp soutBn soutCp soutCn
clear stim_A stim_B timeVEC1 timeVEC


%% Organize variables for SACSCCanal
% SpikeTrains=cell(3,2); % {condition (1,2,3), polarity (plus,minus)}
SpikeTrains={SpikeTrainsA_plus,SpikeTrainsA_minus;SpikeTrainsB_plus,SpikeTrainsB_minus;SpikeTrainsC_plus,SpikeTrainsC_minus};

% specify params to be used
clear paramsIN
paramsIN.durA_msec=dur_sec*1000;
paramsIN.durB_msec=dur_sec*1000;
paramsIN.durC_msec=dur_sec*1000;
paramsIN.CF_A_Hz=CF_A_kHz*1000;
paramsIN.CF_B_Hz=CF_B_kHz*1000;
paramsIN.CF_C_Hz=CF_C_kHz*1000;
% Need to include CF_A, CF_B, CF_C for more generality
paramsIN.MAXspikes=2500;
paramsIN.SNR2use_dB=SNR2use_dB;
%% *MH July 7 2015 - change to longer delays to allow for lower modulations
paramsIN.MAXdelay_sec=1;  % sentence duration is 2.7 sec, so 1 sec delays are about 1/3 of stim duration, which is pushing what we can estimate
%% *MH: HOW TO OPTIMIZE??

[SACSCCfunctions,SACSCCmetrics,paramsOUT] = Library.SACSCCanal_SNRenv(SpikeTrains,paramsIN,1,resultDir,resultPostfix);
%Library.parsave([resultDir 'sac/sac' resultPostfix{2} '.mat'], SACSCCfunctions, SACSCCmetrics, paramsOUT )
Library.parsave([resultDir 'sacMet/sacMet' resultPostfix{2} '.mat'],  SACSCCmetrics )

%%% PSDenv computed in SACSCCanal_SNRenv and passed within SACSCCfunctions
% use BOOTSTRAP REP 6, which is the AVG of 5 bootstrap REPS
% env
PSDenv_S=SACSCCfunctions{6}.PSDenv_A;
PSDenv_N=SACSCCfunctions{6}.PSDenv_B;
PSDenv_SN=SACSCCfunctions{6}.PSDenv_C;
PSDenv_S_noisefloor=SACSCCfunctions{6}.rand.PSDenv_A;
PSDenv_N_noisefloor=SACSCCfunctions{6}.rand.PSDenv_B;
PSDenv_SN_noisefloor=SACSCCfunctions{6}.rand.PSDenv_C;
PSDfreqVEC_Hz=SACSCCfunctions{6}.PSD_freqVEC;

% tfs
tfsCutOff = 1.5;
if CF_A_kHz < tfsCutOff && CF_B_kHz < tfsCutOff && CF_C_kHz < tfsCutOff
    PSDtfs_S=SACSCCfunctions{6}.PSDtfs_A;
    PSDtfs_N=SACSCCfunctions{6}.PSDtfs_B;
    PSDtfs_SN=SACSCCfunctions{6}.PSDtfs_C;
    PSDtfs_S_noisefloor=SACSCCfunctions{6}.rand.PSDtfs_A;
    PSDtfs_N_noisefloor=SACSCCfunctions{6}.rand.PSDtfs_B;
    PSDtfs_SN_noisefloor=SACSCCfunctions{6}.rand.PSDtfs_C;
else  
    PSDtfs_S=0;
    PSDtfs_N=0;
    PSDtfs_SN=0;
    PSDtfs_S_noisefloor=0;
    PSDtfs_N_noisefloor=0;
    PSDtfs_SN_noisefloor=0;
end

%% *MH - look at spectra
if plt
    figure;
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
        CF_C_kHz,SRtype_C,Cohc_C,Cihc_C,OALevel_dBSPL,SNR2use_dB),'FontSize',12)
    Library.saveFigureAs([resultDir 'doc/psd' resultPostfix{1} '.eps'])
    Library.saveFigureAs([resultDir 'doc/png/psd' resultPostfix{1} '.png'])
end
% Divide the FFT output into modulation filters
ModFreq = [1,2,4,8,16,32,64];
PowerModS = zeros(length(ModFreq),1);
PowerModSN = zeros(length(ModFreq),1);
PowerModN = zeros(length(ModFreq),1);
PowerModS_noisefloor = zeros(length(ModFreq),1);
PowerModSN_noisefloor = zeros(length(ModFreq),1);
PowerModN_noisefloor = zeros(length(ModFreq),1);

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



% figure
% scatter(ModFreq, PowerModSN)
% hold on
% scatter(ModFreq, PowerModN)


end


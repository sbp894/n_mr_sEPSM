function [SpikeTrains,paramsIN, mrWindows]=get_model_spiketimes(stim_S, stim_N, stim_SN, A, resultPostfix, C, AN, anal, resultsDir, FixSNlevel, mrWindows) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SET PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CF_kHz = C.cF_i/1e3;
fiberType = C.ftype_i;
Cohc = AN.Cohc(AN.CF==C.cF_i);
Cihc = AN.Cihc(AN.CF==C.cF_i);
OALevel_dBSPL = C.level;
SNR2use_dB = C.snr_i;

% if isfield(ExpControlParams, 'fixSPL')
%     FixSNlevel=ExpControlParams.fixSPL;
% else
%     FixSNlevel=1;
% end
% verbose=anal.verbose;

%% MODEL params - general (for both conditions)
ANmodel_Fs_Hz=100e3;
Nreps=AN.nTrials;

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
%% Stimuli

stim_A_model = stim_S;
stim_B_model= stim_N;
stim_C_model = stim_SN;

dur_sec=min([length(stim_A_model)/ANmodel_Fs_Hz length(stim_B_model)/ANmodel_Fs_Hz length(stim_C_model)/ANmodel_Fs_Hz]);  % Set duration to minimum of all stimuli

%% Scale stimuli to correct OALevel_dBSPL
% Note: Stim scaling for SNRenv study is done to set the SNR by
% adjusting the Noise level, so SIGNAL is same level in all WAV files,
% and SIGNAL should be scaled to provide desired OALlevel, and then all
% three WAV files scaled by the same scale factor.

dBSPL_A_after=20*log10(sqrt(mean(stim_A_model.^2))/(20e-6));

scaleFactor=10^((OALevel_dBSPL-dBSPL_A_after)/20);
stim_A_model=stim_A_model*scaleFactor;
stim_B_model=stim_B_model*scaleFactor;
stim_C_model=stim_A_model+stim_B_model;

if FixSNlevel
    dBSPL_A_after=20*log10(sqrt(mean(stim_A_model.^2))/(20e-6));
    stim_A_model=stim_A_model*10^((OALevel_dBSPL-dBSPL_A_after)/20);
    
    dBSPL_B_after=20*log10(sqrt(mean(stim_B_model.^2))/(20e-6));
    stim_B_model=stim_B_model*10^((OALevel_dBSPL-dBSPL_B_after)/20);
    
    dBSPL_C_after=20*log10(sqrt(mean(stim_C_model.^2))/(20e-6));
    stim_C_model=stim_C_model*10^((OALevel_dBSPL-dBSPL_C_after)/20);
end

%% REFIT and WINDOWwavefile at ANmodel_Fs
% Repeat or truncate waveform to fit requested stimulus duration:
stim_A_model = Library.refit_waveform(stim_A_model,ANmodel_Fs_Hz,dur_sec*1000);
stim_B_model = Library.refit_waveform(stim_B_model,ANmodel_Fs_Hz,dur_sec*1000);
stim_C_model = Library.refit_waveform(stim_C_model,ANmodel_Fs_Hz,dur_sec*1000);

% Window waveform using linear rise/fall:
stim_A_model = Library.window_waveform(stim_A_model,ANmodel_Fs_Hz,dur_sec*1000);
stim_B_model = Library.window_waveform(stim_B_model,ANmodel_Fs_Hz,dur_sec*1000);
stim_C_model = Library.window_waveform(stim_C_model,ANmodel_Fs_Hz,dur_sec*1000);

%% Debugging
% stim_C_model=stim_B_model;

%% Listen to Stimuli at the Levels and SNR played to the model - listen at least once for each new condition to verify
PLAYstim = 0;  % Set to 0 to not hear stimuli
if PLAYstim
    %     PAUSEtime_sec=3; %#ok<*UNRCH>
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


%% Spike generation (from AN model for all 6 conditions: A+,A-,B+,B-,C+,C-)

%% STIM_A = Signal
implnt = 0;     % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse
noiseType = 0;  % 1 for variable fGn (0 for fixed fGn)

% stim_A_plus
vIHC = ANModel.model_IHC(stim_A_model.',CF_A_kHz*1000,1,1/ANmodel_Fs_Hz,dur_sec+0.05,Cohc_A,Cihc_A,AN.species); % 2013 model; "2" for human with BM tuning from Shera et al. (PNAS 2002)
[sout_Ap,meanrate_unad_Ap, ~] = ANModel.model_Synapse(vIHC,CF_A_kHz*1000,1,1/ANmodel_Fs_Hz,SRtype_A,noiseType ,implnt);
% SpikeTrainsA_plus=ensure_no_empty_spikes(meanrate_unad,dur_sec,ANmodel_Fs_Hz,Nreps, anal); %=get_sptimes(meanrate_unad,ANmodel_Fs_Hz,Nreps);
SpikeTrainsA_plus=get_sptimes(meanrate_unad_Ap,ANmodel_Fs_Hz,Nreps);

% stim_A_minus
vIHC = ANModel.model_IHC(-stim_A_model.',CF_A_kHz*1000,1,1/ANmodel_Fs_Hz,dur_sec+0.05,Cohc_A,Cihc_A,AN.species); % 2013 model
[sout_An,meanrate_unad_An, ~] = ANModel.model_Synapse(vIHC,CF_A_kHz*1000,1,1/ANmodel_Fs_Hz,SRtype_A,noiseType ,implnt);
% SpikeTrainsA_minus=ensure_no_empty_spikes(meanrate_unad,dur_sec,ANmodel_Fs_Hz,Nreps, anal); %get_sptimes(meanrate_unad,ANmodel_Fs_Hz,Nreps);
SpikeTrainsA_minus=get_sptimes(meanrate_unad_An,ANmodel_Fs_Hz,Nreps);

%% STIM_B = Noise
% stim_B_plus
vIHC = ANModel.model_IHC(stim_B_model.',CF_B_kHz*1000,1,1/ANmodel_Fs_Hz,dur_sec+0.05,Cohc_B,Cihc_B,AN.species); % 2013 model
[sout_Bp,meanrate_unad_Bp, ~] = ANModel.model_Synapse(vIHC,CF_B_kHz*1000,1,1/ANmodel_Fs_Hz,SRtype_B,noiseType ,implnt);
% SpikeTrainsB_plus=ensure_no_empty_spikes(meanrate_unad,dur_sec,ANmodel_Fs_Hz,Nreps, anal); %get_sptimes(meanrate_unad,ANmodel_Fs_Hz,Nreps);
SpikeTrainsB_plus=get_sptimes(meanrate_unad_Bp,ANmodel_Fs_Hz,Nreps);

% stim_B_minus
vIHC = ANModel.model_IHC(-stim_B_model.',CF_B_kHz*1000,1,1/ANmodel_Fs_Hz,dur_sec+0.05,Cohc_B,Cihc_B,AN.species); % 2013 model
[sout_Bn,meanrate_unad_Bn, ~] = ANModel.model_Synapse(vIHC,CF_B_kHz*1000,1,1/ANmodel_Fs_Hz,SRtype_B,noiseType ,implnt);
% SpikeTrainsB_minus=ensure_no_empty_spikes(meanrate_unad,dur_sec,ANmodel_Fs_Hz,Nreps, anal); %get_sptimes(meanrate_unad,ANmodel_Fs_Hz,Nreps);
SpikeTrainsB_minus=get_sptimes(meanrate_unad_Bn,ANmodel_Fs_Hz,Nreps);

%%%%%%%%%%%%%%%%
%% STIM_C = Signal+Noise
% stim_C_plus
vIHC = ANModel.model_IHC(stim_C_model.',CF_C_kHz*1000,1,1/ANmodel_Fs_Hz,dur_sec+0.05,Cohc_C,Cihc_C,AN.species); % 2013 model
[sout_Cp,meanrate_unad_Cp, ~] = ANModel.model_Synapse(vIHC,CF_C_kHz*1000,1,1/ANmodel_Fs_Hz,SRtype_C,noiseType ,implnt);
% SpikeTrainsC_plus=ensure_no_empty_spikes(meanrate_unad,dur_sec,ANmodel_Fs_Hz,Nreps, anal); %get_sptimes(meanrate_unad,ANmodel_Fs_Hz,Nreps);
SpikeTrainsC_plus=get_sptimes(meanrate_unad_Cp,ANmodel_Fs_Hz,Nreps);

% stim_C_minus
vIHC = ANModel.model_IHC(-stim_C_model.',CF_C_kHz*1000,1,1/ANmodel_Fs_Hz,dur_sec+0.05,Cohc_C,Cihc_C,AN.species); % 2013 model
[sout_Cn,meanrate_unad_Cn, ~] = ANModel.model_Synapse(vIHC,CF_C_kHz*1000,1,1/ANmodel_Fs_Hz,SRtype_C,noiseType ,implnt);
% SpikeTrainsC_minus=ensure_no_empty_spikes(meanrate_unad,dur_sec,ANmodel_Fs_Hz,Nreps, anal); %get_sptimes(meanrate_unad,ANmodel_Fs_Hz,Nreps);
SpikeTrainsC_minus=get_sptimes(meanrate_unad_Cn,ANmodel_Fs_Hz,Nreps);


if anal.plot_meanrate
%     figure(figHandles.meanRates);
    clf;
    ymaxTEMP=-Inf;
    timeVEC1_sec=((1:length(stim_A_model))-1)/ANmodel_Fs_Hz;
    
    subplot(411)
    plot(timeVEC1_sec,stim_C_model,'b'); hold on; plot(timeVEC1_sec,stim_B_model,'r'); plot(timeVEC1_sec,stim_A_model,'g');
    legend('SN','N','S')
    title(sprintf('CF=%.3f kHz; SRtype=%d (1: LSR; 2: MSR; 3: HSR); Cohc=%.2f; Cihc=%.2f; Overall (Signal) Level: %.1f dB SPL; SNR = %.1f dB', CF_C_kHz,SRtype_C,Cohc_C,Cihc_C,OALevel_dBSPL,SNR2use_dB));
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
    
    subplot(412); ylim([0 ymaxTEMP]);
    subplot(413); ylim([0 ymaxTEMP]);
    subplot(414); ylim([0 ymaxTEMP]);
    
    saveas(gcf,[resultsDir 'SpeechSpikePlots_bmp' filesep 'SSPlots' resultPostfix '.bmp']);
    Library.saveFigureAs([resultsDir 'SpeechSpikePlots_eps' filesep 'SSPlots' resultPostfix '.eps']);
    
end

% Library.saveFigureAs([resultsDir 'SpeechSpikePlots_bmp' filesep 'SSPlots' resultPostfix '.bmp']);

SpikeTrains={SpikeTrainsA_plus,SpikeTrainsA_minus;SpikeTrainsB_plus,SpikeTrainsB_minus;SpikeTrainsC_plus,SpikeTrainsC_minus};


%% specify params to be used
clear paramsIN
paramsIN.durA_msec=dur_sec*1000;
paramsIN.durB_msec=dur_sec*1000;
paramsIN.durC_msec=dur_sec*1000;
paramsIN.CF_A_Hz=CF_A_kHz*1000;
paramsIN.CF_B_Hz=CF_B_kHz*1000;
paramsIN.CF_C_Hz=CF_C_kHz*1000;
paramsIN.MAXspikes=2500;
paramsIN.SNR2use_dB=SNR2use_dB;
paramsIN.SCC_onsetIGNORE_sec=anal.onsetIgnore;
%% *MH July 7 2015 - change to longer delays to allow for lower modulations
paramsIN.MAXdelay_sec=1;  % sentence duration is 2.7 sec, so 1 sec delays are about 1/3 of stim duration, which is pushing what we can estimate
paramsIN.FiberType=AN.fiberType.sr;
paramsIN.plt=anal.plot_PSD;
paramsIN.level=1;
paramsIN.spls_s_n_sn=[20*log10(sqrt(mean(stim_A_model.^2))/(20e-6)), 20*log10(sqrt(mean(stim_B_model.^2))/(20e-6)), 20*log10(sqrt(mean(stim_C_model.^2))/(20e-6))];
paramsIN.MINspikes=anal.MINspikes;
mrWindows(mrWindows>(dur_sec-anal.onsetIgnore))=dur_sec-anal.onsetIgnore;
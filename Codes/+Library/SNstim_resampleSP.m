function [stimVEC_S_model,stimVEC_N_model,stimVEC_SN_model]=SNstim_resampleSP(A,B,C,anal)
% This function generates (saves) 3 wave files sampled at 100 kHz so that
% rest of code doesn't have to resample each time.  Also, this allows exact
% SN-N=S, which is fouled up after resampling.

%for Hint sentences:
stim_S_Fname = [A.StimDir 'Stim_S_P.wav'];
% stim_N_Fname = ['stimuli' filesep B.noisePrefix{strcmp(B.noiseTypes,C.noise_i)} B.fileExtension];
stim_N_Fname = sprintf('%sStim%1.0fdB_N_P.wav', A.StimDir ,-6); %SP Hardcoded. should change this.

SNR2use = C.snr_i;
verbose = anal.verbose;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AN-model Stimuli generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[stimVEC_S, FsS_Hz] = audioread(stim_S_Fname);
[stimVEC_N, FsN_Hz] = audioread(stim_N_Fname);

%%%%%%%%%%%
%% RESAMPLE StimVEC_S to ANmodel_Fs_Hz, with original as Fs_Hz
ANmodel_Fs_Hz=A.Fs;
if verbose
    fprintf('... resampling stimVEC_S (%s) from %.f Hz to %.f Hz',stim_S_Fname,FsS_Hz,ANmodel_Fs_Hz);
end

dBSPL_S_before=20*log10(sqrt(mean(stimVEC_S.^2))/(20e-6));
sfreq=FsS_Hz;	   
sfreqNEW=ANmodel_Fs_Hz;

P=round(sfreqNEW/10); 
Q=round(sfreq/10);  %Integers used to up sample
if(P/Q*sfreq~=sfreqNEW && verbose) 
    disp('Integer sfreq conversion NOT exact'); 
end

Nfir=30;  % proportional to FIR filter length used for resampling: higher Nfir, better accuracy & longer comp time
stimVEC_S_model=resample(stimVEC_S,P,Q,Nfir);
dBSPL_S_after=20*log10(sqrt(mean(stimVEC_S_model.^2))/(20e-6));
if abs(dBSPL_S_before-dBSPL_S_after)>0.2;
    error('RESAMPLING CHANGED stim_A dBSPL by %f dB',dBSPL_A_after-dBSPL_A_before);
end

%%%%%%%%%%%
%% RESAMPLE StimVEC_N to ANmodel_Fs_Hz, with original as Fs_Hz
if verbose
    fprintf('... resampling stimVEC_N (%s) from %.f Hz to %.f Hz',stim_N_Fname,FsN_Hz,ANmodel_Fs_Hz);
end
dBSPL_N_before=20*log10(sqrt(mean(stimVEC_N.^2))/(20e-6));

sfreq=FsN_Hz;	  
sfreqNEW=ANmodel_Fs_Hz;
P=round(sfreqNEW/10); 
Q=round(sfreq/10);  %Integers used to up sample
if(P/Q*sfreq~=sfreqNEW && verbose) 
    disp('Integer sfreq conversion NOT exact'); 
end

Nfir=30;  % proportional to FIR filter length used for resampling: higher Nfir, better accuracy & longer comp time
stimVEC_N_model=resample(stimVEC_N,P,Q,Nfir);
dBSPL_N_after=20*log10(sqrt(mean(stimVEC_N_model.^2))/(20e-6));
if abs(dBSPL_N_before-dBSPL_N_after)>0.2;
    error('RESAMPLING CHANGED stim_A dBSPL by %f dB',dBSPL_A_after-dBSPL_A_before);
end

% Adjust the length of the noise file. It should start 1 sec. before the speech signal
% and end 0.6 sec. after; this equals 1.6 x Fs samples. Start with
% a random sample:
startSampleMax =1;  length(stimVEC_N_model) - length(stimVEC_S_model);
rng('default');
startSample = ceil(rand(1,1)*startSampleMax);
if startSample==0
    startSample=1;
end
stimVEC_N_model = stimVEC_N_model(startSample:startSample+length(stimVEC_S_model)-1);

%% Set SNR:
RMS_S = rms(stimVEC_S_model);
RMS_N = rms(stimVEC_N_model);
SNRorig = 20*log10(RMS_S/RMS_N);
if verbose
    fprintf('..SNR2use = %.2f dB',SNR2use);
    sprintf('..SNR before = %.2f dB',SNRorig);
end

%% %Adjust the Noise to achieve desired SNR
% Note: Stim scaling is done to set the SNR by adjusting the Noise level,
% so SIGNAL is same level in all WAV files, and SIGNAL should be scaled to
% provide deisred OALlevel, and then all three WAV files scaled by the same
% scale factor.
stimVEC_N_model=stimVEC_N_model*10^((SNRorig-SNR2use)/20);

%% Save new WAV files
if length(stimVEC_S_model) ~= length(stimVEC_N_model)
    numZeros = length(stimVEC_N_model)-length(stimVEC_S_model);
    stimVEC_S_model = stimVEC_S_model(1:end + numZeros);
end
stimVEC_SN_model = stimVEC_S_model+stimVEC_N_model;  % Create SN condition after resampling and after SNR is set.

%% Compute Acoustic SNR:
RMS_S = rms(stimVEC_S_model);
RMS_N = rms(stimVEC_N_model);

stimSNR = 20*log10(RMS_S/RMS_N);
if verbose
    fprintf('***Acoustic SNR = %.2f dB',stimSNR);
end
if abs(SNR2use-stimSNR)>0.1
    error('Desired SNR not generated')
end
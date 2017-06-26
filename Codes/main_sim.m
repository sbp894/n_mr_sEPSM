%% 
clear;
close all;
clc;

%% Define model parameters and generate stimuli
verbose = 0;
Fs = 100000; % Model sampling frequency
dt = 1/Fs;

%% create progress dir: each itrn of the || loop creats a file in this directory. The files serve both as progress indicator and as a way start the simulation again at the point where it stopped in the event of an crash
progressDir = ['progress' filesep];
if ~exist(progressDir,'dir')
    mkdir(progressDir)
end
resultsDir=Simulation.create_result_dir; % different prefix/postfix for data analysis or simulation

%% Get Parameters
[A,B]=Simulation.get_speech_params(Fs);
AN=Simulation.get_AN_params(verbose);
anal=Simulation.get_anal_params(Fs,AN,resultsDir);
cs=Simulation.get_conditions(A,B,AN,resultsDir);

for c_i = 1 %: size(cs,1)
    C = struct(); % necessary to create structure in parfor
    C.cF_i = cs(c_i,1);
    C.sentence_i = cs(c_i,2);
    C.snr_i = cs(c_i,3);
    C.noise_i = cs(c_i,4);
    C.level_i = cs(c_i,5);
    
    resultPostfix2 = sprintf(anal.resultFileName,C.cF_i,C.sentence_i,C.snr_i, C.noise_i, C.level_i);
    if ~exist([progressDir filesep resultPostfix2 '.mat'],'file')
        % This function generates 3 wave files sampled at 100 kHz so that rest of code doesn't have to resample each time. Also, this allows exact relation of S=SN-N, which is fouled up is SN=S+N is done before resampling.
        [stim_S,stim_N,stim_SN] = Library.SNstim_resample(A,B,C,anal);

        [SpikeTrains,paramsIN]=Simulation.get_model_spiketimes(stim_S, stim_N, stim_SN, A, B, C, AN, anal);
        
        [PSDenv_STRUCT,PSDtfs_STRUCT,PowerMod_STRUCT,PowerTfs_STRUCT] = Library.sumcors_sp(A, B, C, AN, anal ,SpikeTrains,paramsIN);
        
        Simulation.save_analysis_results(PSDenv_STRUCT,PSDtfs_STRUCT,PowerMod_STRUCT,PowerTfs_STRUCT,resultsDir,resultPostfix2);
        
        Simulation.update_progress(progressDir,resultPostfix2,A,B,C,AN,anal,cs);    
    end
end
time_taken=toc;
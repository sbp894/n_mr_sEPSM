%%
clear;
close all;
% clc;

tic;
%%
Simulation1DataAnal0=0;
DataDir='D:\Study Stuff\Matlab\NelData\SP-2016_10_28-Q285_AN_Normal_forsnrenvcodes';
verbose = 0;

%% create progress dir: each itrn of the || loop creats a file in this directory. The files serve both as progress indicator and as a way start the simulation again at the point where it stopped in the event of an crash
% progressDir = ['progress' filesep];
% if ~exist(progressDir,'dir')
%     mkdir(progressDir)
% end


%% Get Parameters
if Simulation1DataAnal0
    resultsDir=Library.create_output_dir(Simulation1DataAnal0,datestr(now,'yyyymmdd')); %#ok<*UNRCH> % Create directories
    Fs = 100000; % Model sampling frequency
    dt = 1/Fs;
    [A,B]=Simulation.get_speech_params(Fs);
    AN=Simulation.get_AN_params(verbose);
    anal=Simulation.get_anal_params(Fs,AN,resultsDir);
    cs=Simulation.get_conditions(A,B,AN,resultsDir);
    MaxIter=size(cs,1);
    CFs=AN.CF;

else
    resultsDir=Library.create_output_dir(Simulation1DataAnal0,DataDir(strfind(DataDir,fileparts(DataDir))+length(fileparts(DataDir))+1:end)); % Create directories
    spike_data=load_data(DataDir,resultsDir);
    [A,B]=DataAnal.get_speech_params;
    anal=DataAnal.get_anal_params(resultsDir);
    MaxIter=length(spike_data);
end


parfor condition_var = 1 : MaxIter
    % parfor condition_var = 1 : MaxIter
    if Simulation1DataAnal0
        C = struct(); % necessary to create structure in parfor
        C.cF_i = CFs(cs(condition_var,1));
        C.sentence_i = 1;%cs(c_i,2);
        C.snr_i = B.SNR(cs(condition_var,3));
        C.noise_i = 'SSN';%cs(c_i,4);
        
    else
        C = struct(); % necessary to create structure in parfor
        C.cF_i=spike_data(condition_var).CF;
        C.sentence_i = 1;
        C.snr_i = spike_data(condition_var).SNR;
        C.noise_i = 'SSN';
    end
    
    resultPostfix = sprintf(resultTxt,       C.cF_i/1e3,  C.sentence_i,    C.snr_i,  C.noise_i);
    
    if ~exist([resultsDir 'progress' filesep resultPostfix '.mat'],'file')
        
        % This function generates 3 wave files sampled at 100 kHz so that rest of code doesn't have to resample each time. Also, this allows exact relation of S=SN-N, which is fouled up is SN=S+N is done before resampling.
        if Simulation1DataAnal0
            [stim_S,stim_N,stim_SN] = Library.SNstim_resample(A,B,C,anal);
            [SpikeTrains,paramsIN]=Simulation.get_model_spiketimes(stim_S, stim_N, stim_SN, A, B, C, AN, anal);
        else
            SpikeTrains=spike_data(condition_var).SpikeTrains;
            paramsIN=spike_data(condition_var).paramsIN;
        end
        
        [PSDenv_STRUCT,PSDtfs_STRUCT,PowerMod_STRUCT,PowerTfs_STRUCT] = Library.sumcors_sp(SpikeTrains,paramsIN, resultsDir,resultPostfix);
        
        Simulation.save_analysis_results(PSDenv_STRUCT,PSDtfs_STRUCT,PowerMod_STRUCT,PowerTfs_STRUCT,resultsDir,resultPostfix,paramsIN);
        Simulation.update_progress(resultsDir,resultPostfix,MaxIter);
    end
end

parse_saved_data_for_SNRenv(resultsDir);
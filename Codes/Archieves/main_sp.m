%%
% function main_sp(DataDir):
%       Data Analysis: Input DataDir name, the function looks for the directory under NELData and does the analysis
% function main_sp(Simulation1DataAnal0):
%       if 1, Simulation.
%       if 0, Data Analysis, User would be asked to input directory.
% function main_sp():
%       Default: Simulation
%
% Created by SP [3/1/16]

%% Set up Conditions
function main_sp(varargin)

if nargin==0
    Simulation1DataAnal0=1;
elseif ischar(varargin{1})
    Simulation1DataAnal0=0;
    DataDir=strcat(pwd, '\NELData\',varargin{1});
elseif varargin{1}==0 || varargin{1}==1
    Simulation1DataAnal0=varargin{1};
    if ~Simulation1DataAnal0
        DataDir=uigetdir([pwd '\NELData']);
    end
else
    error('Type help main_sp to see usage');
end
verbose=0;

%% Get Parameters
if Simulation1DataAnal0
    resultsDir=Library.create_output_dir(Simulation1DataAnal0,datestr(now,'yyyymmdd')); % Create directories
    Fs = 100000; % Model sampling frequency
    [A,B]=Simulation.get_speech_params(Fs);
    AN=Simulation.get_AN_params;
    anal=Simulation.get_anal_params(Fs,AN,resultsDir);
    cndts=Simulation.get_conditions(A,B,AN,resultsDir);
    MaxIter=size(cndts,1);
    CFs=AN.CF;
    spike_data={};
    StimData=[];
    
else
    % Create directories
    resultsDir=Library.create_output_dir(Simulation1DataAnal0,DataDir(strfind(DataDir,fileparts(DataDir))+length(fileparts(DataDir))+1:end));
    [spike_data,StimData]=load_data(DataDir,resultsDir);
    [A,B]=DataAnal.get_speech_params;
    anal=DataAnal.get_anal_params(resultsDir);
    MaxIter=length(spike_data);
    CFs={};
    cndts={};
    AN={};
end

resultTxt=anal.resultTxt;

% for condition_var = 6 : 8
for condition_var = 1 : MaxIter
    if Simulation1DataAnal0
        C=getCsim(cndts,condition_var,B,CFs);
    else
        C=getCdata(spike_data,condition_var);
    end
    
    resultPostfix = sprintf(resultTxt,       C.cF_i/1e3,  C.sentence_i,    C.snr_i,  C.noise_i);
    if ~exist([resultsDir 'progress' filesep resultPostfix '.mat'],'file')
        
        
        if Simulation1DataAnal0
            % This function generates 3 wave files sampled at 100 kHz so that rest of code doesn't have to resample each time....
            % Also, this allows exact relation of S=SN-N, which is fouled up is SN=S+N is done before resampling.
            [stim_S,stim_N,stim_SN] = Library.SNstim_resample(A,B,C,anal);
            [SpikeTrains,paramsIN]=Simulation.get_model_spiketimes(stim_S, stim_N, stim_SN, A, resultPostfix, C, AN, anal, resultsDir);
        else
            SpikeTrains=spike_data(condition_var).SpikeTrains;
            plot_assist=data_plot_assist(spike_data,StimData,condition_var);
            plot_stimulus_spike_SNSN(SpikeTrains,plot_assist,resultsDir,resultPostfix);
            paramsIN=spike_data(condition_var).paramsIN;
        end
        
        if ~isempty(SpikeTrains)
            [PSDenv_STRUCT,PSDtfs_STRUCT,PowerMod_STRUCT,PowerTfs_STRUCT] = Library.sumcors_sp(SpikeTrains,paramsIN, resultsDir,resultPostfix);
            %             PowerMod_STRUCT.PowerModS
            Simulation.save_analysis_results(PSDenv_STRUCT,PSDtfs_STRUCT,PowerMod_STRUCT,PowerTfs_STRUCT,resultsDir,resultPostfix,paramsIN);
            Simulation.update_progress(resultsDir,resultPostfix,MaxIter);
        else
            disp(['Whoa!' num2str(condition_var)]);
        end
    end
end

parse_saved_data_for_SNRenv2(resultsDir);
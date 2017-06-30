function [resultsDir,resultTxt]=Analyze_Data(DataDir)

resultsDir=Library.create_output_dir(0,DataDir(strfind(DataDir,fileparts(DataDir))+length(fileparts(DataDir))+1:end)); % Create directories
[spike_data,StimData,StimsFNames]=DataAnal.load_data(DataDir,resultsDir); %#ok<ASGLU>
save([resultsDir 'SpikeStimulusData.mat'],'spike_data','StimsFNames');

anal=DataAnal.get_anal_params(resultsDir);
MaxIter=length(spike_data);
resultTxt=anal.resultTxt;

% paramsIN=spike_data{end}.paramsIN;
for condition_var = 1 : MaxIter
    
    C=DataAnal.getCdata(spike_data,condition_var);
    resultPostfix = sprintf(resultTxt,       C.cF_i/1e3,  C.sentence_i,  C.noise_i, C.level,    C.snr_i);
    if spike_data(condition_var).nReps>=14
        
        
        if ~exist([resultsDir 'progress' filesep resultPostfix '.mat'],'file')
            
            SpikeTrains=spike_data(condition_var).SpikeTrains;
            plot_assist=DataAnal.data_plot_assist(spike_data,StimData,condition_var);
            DataAnal.plot_stimulus_spike_SNSN(SpikeTrains,plot_assist,resultsDir,resultPostfix);
            
            if ~isempty(SpikeTrains)
                paramsIN=spike_data(condition_var).paramsIN;
                
                [PSDenv_STRUCT,PSDtfs_STRUCT,PowerMod_STRUCT,PowerTfs_STRUCT] = ...
                    Library.sumcors_bootstrap(SpikeTrains,paramsIN, resultsDir,resultPostfix);
                
                Library.save_analysis_results(PSDenv_STRUCT,PSDtfs_STRUCT,PowerMod_STRUCT,PowerTfs_STRUCT,resultsDir,resultPostfix,paramsIN);
                Library.update_progress(resultsDir,resultPostfix,MaxIter);
            else
                disp(['Whoa!' num2str(condition_var)]);
            end
        end
    else
        fprintf('File: %s skipped\n',resultPostfix);
    end
end
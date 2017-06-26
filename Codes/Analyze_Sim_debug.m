function [resultsDir,resultTxt]=Analyze_Sim_debug(ExpControlParams, RootOUTPUTDir)

resultsDir=Library.create_output_dir(1,'SimOut', RootOUTPUTDir); % Create directories
Fs = 100e3; % Model sampling frequency
[A,B]=Simulation.get_speech_params(Fs, ExpControlParams);
AN=Simulation.get_AN_params(ExpControlParams);
anal=Simulation.get_anal_params(Fs,AN,resultsDir);
cndts=Simulation.get_conditions(A,B,AN,resultsDir);
MaxIter=size(cndts,1);
CFs=AN.CF;
resultTxt=anal.resultTxt;

FixSNlevel=ExpControlParams.fixSPL;
mrWindows=ExpControlParams.mrWindows;
BootstrapLoopMax=ExpControlParams.BootstrapLoopMax;%300;
BootstrapLoopReport=ExpControlParams.BootstrapLoopReport;


for condition_var = 1:MaxIter
    
    C=getCsim(cndts,condition_var,A,B,AN, CFs);
    resultPostfix = sprintf(resultTxt,       C.cF_i/1e3,  C.sentence_i,  C.noise_i, C.level, C.snr_i, C.ftype_i);
%     show_current_parloop(resultsDir ,resultPostfix);

    if ~exist([resultsDir 'progress' filesep resultPostfix '.mat'],'file') && exist([resultsDir 'current' filesep resultPostfix '.mat'],'file')
        
        [stim_S,stim_N,stim_SN] = Library.SNstim_resample(A,B,C,anal);
        %         [stim_S,stim_N,stim_SN] = Library.SNstim_resampleSP(A,B,C,anal);
        [SpikeTrains,paramsIN, mrWindowsUpd]=Simulation.get_model_spiketimes(stim_S, stim_N, stim_SN, A, resultPostfix, C, AN, anal, resultsDir, FixSNlevel, mrWindows);
        
        [mrSPikeTrains, paramsIN, PowerModCell]=get_mrSPikeTrains(SpikeTrains, mrWindowsUpd, anal.onsetIgnore, paramsIN);
        
        
        if ~isempty(SpikeTrains)
            for windowVar=1:length(mrSPikeTrains)
%                 paramsIN.SCCdur_sec=ExpControlParams.mrWindows(windowVar);
                for loopVar=1:length(mrSPikeTrains{windowVar})
                    paramsIN.curBoundaries=paramsIN.mrTimeBoundaries{windowVar}{loopVar};
                    paramsIN.SCCdur_sec=diff(paramsIN.curBoundaries);
                    paramsIN.MAXdelay_sec=.8*paramsIN.SCCdur_sec;
                    [~,~,PowerMod_STRUCT,~] = Library.sumcors_bootstrap(mrSPikeTrains{windowVar}{loopVar},paramsIN, resultsDir,resultPostfix, BootstrapLoopMax, BootstrapLoopReport);
                    PowerModCell{windowVar}{loopVar}=PowerMod_STRUCT;
                    fprintf('%i/%i :: %i/%i\n',windowVar, length(mrSPikeTrains), loopVar, length(mrSPikeTrains{windowVar}));
                end
            end
            %             Simulation.save_analysis_results(PSDenv_STRUCT,PSDtfs_STRUCT,PowerMod_STRUCT,PowerTfs_STRUCT,resultsDir,resultPostfix,paramsIN);
            mr_save_analysis_results(PowerModCell,resultsDir,resultPostfix,paramsIN);
            Simulation.update_progress(resultsDir,resultPostfix,MaxIter);
        else
            disp(['Whoa!' num2str(condition_var)]);
        end
    end
end
function [resultsDir,resultTxt]=Analyze_Sim(ExpControlParams, RootOUTPUTDir, figHandles)

resultsDir=Library.create_output_dir(1,'SimOut', RootOUTPUTDir); % Create directories
[A,B]=Simulation.get_speech_params(ExpControlParams);
AN=Simulation.get_AN_params(ExpControlParams);
anal=Simulation.get_anal_params(ExpControlParams,AN,resultsDir,figHandles);
cndts=Simulation.get_conditions(A,B,AN,resultsDir);
MaxIter=size(cndts,1);
CFs=AN.CF;
resultTxt=anal.resultTxt;


parfor condition_var = 1:MaxIter
    
    C=Simulation.getCsim(cndts,condition_var,A,B,AN, CFs);
    resultPostfix = sprintf(resultTxt,       C.cF_i/1e3,  C.sentence_i,  C.noise_i, C.level, C.snr_i, C.ftype_i);
    
    if ~exist([resultsDir 'progress' filesep resultPostfix '.mat'],'file')  %&& exist([resultsDir 'current' filesep resultPostfix '.mat'],'file')
        
        [stim_S,stim_N,stim_SN] = Library.SNstim_resample(A,B,C,anal);
        [SpikeTrains,paramsIN, mrWindowsUpd]=Simulation.get_model_spiketimes(stim_S, stim_N, stim_SN, A, resultPostfix, C, AN, anal, resultsDir, ExpControlParams.mrWindows);
        
        [mrSPikeTrains, paramsIN, PowerModCell]=Simulation.get_mrSPikeTrains(SpikeTrains, mrWindowsUpd, anal.onsetIgnore, paramsIN);

        for windowVar=1:length(mrSPikeTrains)
            for loopVar=1:length(mrSPikeTrains{windowVar})
               
                paramsIN.curBoundaries=paramsIN.mrTimeBoundaries{windowVar}{loopVar};
                paramsIN.SCCdur_sec=diff(paramsIN.curBoundaries);
                paramsIN.MAXdelay_sec=.8*paramsIN.SCCdur_sec;
                
                [~,~,PowerMod_STRUCT,~] = Library.sumcors_bootstrap(mrSPikeTrains{windowVar}{loopVar},paramsIN, resultsDir,resultPostfix, anal);
                PowerModCell{windowVar}{loopVar}=PowerMod_STRUCT;
                fprintf('%i/%i :: %i/%i \n',windowVar, length(mrSPikeTrains), loopVar, length(mrSPikeTrains{windowVar}));
            end
        end
        Library.mr_save_analysis_results(PowerModCell,resultsDir,resultPostfix,paramsIN);
        Library.update_progress(resultsDir,resultPostfix,MaxIter);
    end
end
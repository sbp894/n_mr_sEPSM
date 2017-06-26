function [resultsDir,resultTxt]=Analyze_Sim(ExpControlParams, RootOUTPUTDir)

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
winCorr0Add1Mul=ExpControlParams.winCorr0Add1Mul;


parfor condition_var = 1:MaxIter
    
    C=getCsim(cndts,condition_var,A,B,AN, CFs);
    resultPostfix = sprintf(resultTxt,       C.cF_i/1e3,  C.sentence_i,  C.noise_i, C.level, C.snr_i, C.ftype_i);
    
    if ~exist([resultsDir 'progress' filesep resultPostfix '.mat'],'file')  %&& exist([resultsDir 'current' filesep resultPostfix '.mat'],'file')
        
        curFile=show_current_parloop(resultsDir ,resultPostfix);
        [stim_S,stim_N,stim_SN] = Library.SNstim_resample(A,B,C,anal);
        [SpikeTrains,paramsIN, mrWindowsUpd]=Simulation.get_model_spiketimes(stim_S, stim_N, stim_SN, A, resultPostfix, C, AN, anal, resultsDir, FixSNlevel, mrWindows);
        
        [mrSPikeTrains, paramsIN, PowerModCell]=get_mrSPikeTrains(SpikeTrains, mrWindowsUpd, anal.onsetIgnore, paramsIN);
        Library.parsave([resultsDir 'current' filesep 'spkData_' resultPostfix '.mat'], SpikeTrains, mrSPikeTrains);

        for windowVar=1:length(mrSPikeTrains)
            for loopVar=1:length(mrSPikeTrains{windowVar})
               
                paramsIN.curBoundaries=paramsIN.mrTimeBoundaries{windowVar}{loopVar};
                paramsIN.SCCdur_sec=diff(paramsIN.curBoundaries);
                paramsIN.MAXdelay_sec=.8*paramsIN.SCCdur_sec;
                
                [~,~,PowerMod_STRUCT,~] = Library.sumcors_bootstrap...
                    (mrSPikeTrains{windowVar}{loopVar},paramsIN, resultsDir,resultPostfix, BootstrapLoopMax, BootstrapLoopReport, winCorr0Add1Mul);
                PowerModCell{windowVar}{loopVar}=PowerMod_STRUCT;
%                 s=whos('PowerModCell');
                fprintf('%i/%i :: %i/%i \n',windowVar, length(mrSPikeTrains), loopVar, length(mrSPikeTrains{windowVar}));
            end
        end
        mr_save_analysis_results(PowerModCell,resultsDir,resultPostfix,paramsIN);
        Simulation.update_progress(resultsDir,resultPostfix,MaxIter);
        delete(curFile);
    end
end
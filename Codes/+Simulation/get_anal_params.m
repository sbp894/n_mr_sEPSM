function anal=get_anal_params(ExpControlParams,AN,resultsDir,figHandles)
%% %%%%%%%%%%%%%%%%%%% SAC/SCC Analysis Parameters %%%%%%%%%%%%%%%%%%%%
anal.binWidth=1/50000;
anal.onsetIgnore=550e-3;
anal.maxDelaySec=1; 
anal.maxSpikes=2500;
anal.fs=ExpControlParams.fs; % unnecessary? 
anal.nTrials=AN.nTrials;
anal.plt=1;
anal.ModFreqs=[1,2,4,8,16,32,64];
anal.resultsDir=resultsDir;
anal.resultTxt='CF_%1.4fk_Sent_%i_%s_level_%1.2f_SNR%i_fType%1.0f';
anal.verbose=0;
anal.plot=0;
anal.FsPSD=10e3;

% anal.plot_meanrate=0;
% anal.plot_modPower=0;
% anal.plot_PSD=0;
anal.figHandles=figHandles;
anal.MINspikes=15;

anal.FixSNlevel=ExpControlParams.fixSPL;
anal.BootstrapLoopMax=ExpControlParams.BootstrapLoopMax;%300;
anal.BootstrapLoopReport=ExpControlParams.BootstrapLoopReport;
anal.nPSDs2Avg=ExpControlParams.nPSDs2Avg;
anal.winCorr0Add1Mul=ExpControlParams.winCorr0Add1Mul;

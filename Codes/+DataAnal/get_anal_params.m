function anal=get_anal_params(resultsDir)

%% %%%%%%%%%%%%%%%%%%% SAC/SCC Analysis Parameters %%%%%%%%%%%%%%%%%%%%
anal.binWidth = 1/50000;
anal.onsetIgnore = 50e-3;
anal.maxLag = 1.5;
anal.maxSpikes = 25000;
% anal.fs = Fs;
anal.nTrials = 1;
anal.plt = 1;
anal.ModFreq = [1,2,4,8,16,32,64];
anal.resultsDir = resultsDir;
anal.resultTxt = 'CF_%1.4fk_Sent_%i_%s_level_%1.2f_SNR%i';
% anal.resultTxt = 'CF_%1.0fk_Sent_%i_level_%1.0f_SNR%i';
anal.verbose = 0;
anal.plot = 1 ;

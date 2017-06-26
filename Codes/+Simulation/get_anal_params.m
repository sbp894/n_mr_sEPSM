function anal=get_anal_params(Fs,AN,resultsDir)

%% %%%%%%%%%%%%%%%%%%% SAC/SCC Analysis Parameters %%%%%%%%%%%%%%%%%%%%
anal.binWidth=1/50000;
anal.onsetIgnore=550e-3;
anal.maxLag=1.2;
anal.maxSpikes=2500;
anal.fs=Fs;
anal.nTrials=AN.nTrials;
anal.plt=1;
anal.ModFreq=[1,2,4,8,16,32,64];
anal.resultsDir=resultsDir;
anal.resultTxt='CF_%1.4fk_Sent_%i_%s_level_%1.2f_SNR%i_fType%1.0f';
anal.verbose=0;
anal.plot=0;
anal.FsPSD=10e3;
anal.MINspikes=15;

anal.plot_meanrate=0;
anal.plot_modPower=0;
anal.plot_PSD=0;

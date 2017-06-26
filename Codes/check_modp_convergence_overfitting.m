%%
% function master_main(DataDir):
%       Data Analysis: Input DataDir name, the function looks for the directory under NELData and does the analysis
% function master_main(Simulation1DataAnal0):
%       if 1, Simulation.
%       if 0, Data Analysis, User will be asked for input directory.
% function master_main():
%       Default: Simulation
%
% Created by SP [5/18/16]

%% Set up Conditions
% function check_modp_convergence_overfitting(varargin)
clear all; %#ok<CLALL>

global RootOUTPUTDir RootCodesDir ExpControlParams

RootCodesDir= [pwd filesep];
RootOUTPUTDir=[fileparts(pwd) filesep 'OUTPUT' filesep];

global figHandles
figHandles.PSDplot=11;
figHandles.SACSCC=2;
figHandles.meanRates=1;
figHandles.modPPlots=13;


cd(RootCodesDir);


%%

ExpControlParams.SNR=6;
ExpControlParams.level=65;

ExpControlParams.noiseTypes={'SSN'};
ExpControlParams.noisePrefix = {['noise' filesep 'ssn_simulation_dtu']};    %     B.noiseTypes = {'SAM'};

% ExpControlParams.noiseTypes={'SAM', 'SSN'};
% ExpControlParams.noisePrefix = {['noise' filesep 'ssn_simulation_dtu'], ['noise' filesep 'sam_simulation_dtu']};    %     B.noiseTypes = {'SAM'};

% ExpControlParams.noiseTypes={'SAM'};
% ExpControlParams.noisePrefix = {['noise' filesep 'sam_simulation_dtu']};

ExpControlParams.fiberType=2; %1:3; % L/M/H <--> 1/2/3
ExpControlParams.CF=[1]*1e3; %logspace(log10(125), log10(8e3), 21);
ExpControlParams.species=2;    % 1 for cat (2 for human with Shera et al. tuning; 3 for human with Glasberg & Moore tuning)
ExpControlParams.sentences=3;

ExpControlParams.nRep=30;
ExpControlParams.BootstrapLoopMax=24;
ExpControlParams.BootstrapLoopReport=60;
ExpControlParams.nPSDs2Avg=12;


resultsDir=Library.create_output_dir(-1,datestr(now,'yyyymmdd')); % Create directories
Fs = 100e3; % Model sampling frequency
[A,B]=Simulation.get_speech_params(Fs);
AN=Simulation.get_AN_params;
anal=Simulation.get_anal_params(Fs,AN,resultsDir);
cndts=Simulation.get_conditions(A,B,AN,resultsDir);
MaxIter=size(cndts,1);
CFs=AN.CF;
resultTxt=anal.resultTxt;

for condition_var = 1:MaxIter
    C=getCsim(cndts,condition_var,A,B,AN, CFs);
    resultPostfix = sprintf(resultTxt,       C.cF_i/1e3,  C.sentence_i,  C.noise_i, C.level, C.snr_i, C.ftype_i);
    
    if ~exist([resultsDir 'progress' filesep resultPostfix '.mat'],'file')
        
        [stim_S,stim_N,stim_SN] = Library.SNstim_resample(A,B,C,anal);
        %         [stim_S,stim_N,stim_SN] = Library.SNstim_resampleSP(A,B,C,anal);
        [SpikeTrains,paramsIN]=Simulation.get_model_spiketimes(stim_S, stim_N, stim_SN, A, resultPostfix, C, AN, anal, resultsDir);
        
        if ~isempty(SpikeTrains)
            [PSDenv_STRUCT,PSDtfs_STRUCT,PowerMod_STRUCT,PowerTfs_STRUCT] = ...
                Library.sumcors_bootstrap(SpikeTrains,paramsIN, resultsDir,resultPostfix);
            
            Simulation.save_analysis_results(PSDenv_STRUCT,PSDtfs_STRUCT,PowerMod_STRUCT,PowerTfs_STRUCT,resultsDir,resultPostfix,paramsIN);
            Simulation.update_progress(resultsDir,resultPostfix,MaxIter);
        else
            disp(['Whoa!' num2str(condition_var)]);
        end
    end
end

figure(1);
clf; hold on; 
co=get(gca, 'colororder');
h1=errorbar(mean(PowerMod_STRUCT.PowerModS), std(PowerMod_STRUCT.PowerModS));
set(h1, 'color', co(1,:));
h2=errorbar(mean(PowerMod_STRUCT.PowerModS_noisefloor), std(PowerMod_STRUCT.PowerModS_noisefloor));
set(h2, 'color', co(1,:), 'LineStyle', '--');
h3=errorbar(mean(PowerMod_STRUCT.PowerModSN), std(PowerMod_STRUCT.PowerModSN));
set(h3, 'color', co(2,:));
h4=errorbar(mean(PowerMod_STRUCT.PowerModSN_noisefloor), std(PowerMod_STRUCT.PowerModSN_noisefloor));
set(h4, 'color', co(2,:), 'LineStyle', '--');
h5=errorbar(mean(PowerMod_STRUCT.PowerModN), std(PowerMod_STRUCT.PowerModN));
set(h5, 'color', co(3,:));
h6=errorbar(mean(PowerMod_STRUCT.PowerModN_noisefloor), std(PowerMod_STRUCT.PowerModN_noisefloor));
set(h6, 'color', co(3,:), 'LineStyle', '--');

legend('PowerModS','PowerModS-NF','PowerModSN','PowerModSN-NF','PowerModN','PowerModN-NF');

% parse_saved_data_for_SNRenvDTU;
% parse_saved_data_for_SNRenvDTU(resultsDir,resultTxt);
parse_saved_data_for_SNRenv3(resultsDir,resultTxt);
% parse_saved_data_for_SNRenv3;

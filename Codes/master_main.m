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
function master_main(varargin)
clear all; %#ok<CLALL>

% global RootDataDir RootOUTPUTDir RootCodesDir ExpControlParams paramsIfCallingFromOutside

RootCodesDir= [pwd filesep];
RootDataDir=[fileparts(pwd) filesep 'MATData' filesep];
RootOUTPUTDir=[pwd filesep 'OUTPUT' filesep];

% global figHandles
% figHandles.PSDplot=11;
% figHandles.SACSCC=2;
% figHandles.meanRates=1;
% figHandles.modPPlots=13;


cd(RootCodesDir);

if nargin==0
    Simulation1DataAnal0=1;
elseif ischar(varargin{1})
    Simulation1DataAnal0=0;
    DataDir=strcat(RootDataDir ,varargin{1});
elseif varargin{1}==0 || varargin{1}==1
    Simulation1DataAnal0=varargin{1};
    if ~Simulation1DataAnal0
        DataDir=uigetdir(RootDataDir);
    end
elseif isnumeric(varargin{1})
    Simulation1DataAnal0=0;
    ChinID=varargin{1};
    allfiles=dir(sprintf('%s*%d*',RootDataDir,ChinID));
    if ~isempty(allfiles)
        DataDir=[RootDataDir allfiles.name];
    end
else
    error('Type help master_main to see usage');
end

%%
if Simulation1DataAnal0
        ExpControlParams.SNR=0; %-9:3:0;
        ExpControlParams.level=65;
        
        ExpControlParams.noiseTypes={'SAM'};
%         ExpControlParams.noiseTypes={'SSN'};        
%         ExpControlParams.noiseTypes={'SAM', 'SSN'};
        
        
        ExpControlParams.noisePrefix = cell(1, length(ExpControlParams.noiseTypes));
        for i=1:length(ExpControlParams.noiseTypes)
            ExpControlParams.noisePrefix{i} = ['noise' filesep lower(ExpControlParams.noiseTypes{i}) '_simulation_dtu'];
        end
        
        ExpControlParams.fiberType=2; %1:3; % L/M/H <--> 1/2/3
        ExpControlParams.CF=1e3; %logspace(log10(125), log10(8e3), 21);
        ExpControlParams.species=2;    % 1 for cat (2 for human with Shera et al. tuning; 3 for human with Glasberg & Moore tuning)
        ExpControlParams.sentences=6; %1:10;
       ExpControlParams.ohcLoss_dB=0;
       ExpControlParams.ihcLoss_dB=0;
     
        ExpControlParams.nRep=60;
        ExpControlParams.BootstrapLoopMax=24;
        ExpControlParams.BootstrapLoopReport=60;
        ExpControlParams.nPSDs2Avg=12;
        ExpControlParams.fixSPL=0;
        ExpControlParams.winCorr0Add1Mul=1;
        ExpControlParams.modFreqs=[1 2 4 8 16 32 64];
        ExpControlParams.mrWindows=2.^10./ExpControlParams.modFreqs*1e-3;
else
    clear global ExpControlParams;
end

%%
if Simulation1DataAnal0 % Simulate
    [resultsDir,resultTxt]=Analyze_Sim(ExpControlParams, RootOUTPUTDir);
    save([resultsDir 'ExpControlParams.mat'],'ExpControlParams');
else
    [resultsDir,resultTxt]=Analyze_Data(DataDir);
end

% parse_saved_data_for_SNRenv2(resultsDir,resultTxt);
if ~strcmp(resultsDir(end),filesep)
    resultsDir=[resultsDir filesep];
end

% parse_saved_data_for_SNRenvDTU;
% parse_saved_data_for_SNRenvDTU(resultsDir,resultTxt);
% parse_saved_data_for_SNRenv_mr(resultsDir,resultTxt);
if exist('ChinID','var')
    %     create_summary(ChinID);
    % %     plot_summary(ChinID);
end
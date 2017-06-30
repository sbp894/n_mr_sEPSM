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

RootDataDir=[fileparts(pwd) filesep 'MATData' filesep];

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
RootOUTPUTDir=[pwd filesep 'OUTPUT' filesep];
figHandles=Library.get_figHandles;

if Simulation1DataAnal0 % Simulate
    ExpControlParams=Simulation.get_ExpControlParams;
    [resultsDir,resultTxt]=Simulation.Analyze_Sim(ExpControlParams, RootOUTPUTDir, figHandles);
    save([resultsDir 'ExpControlParams.mat'],'ExpControlParams');
else
    [resultsDir,resultTxt]=DataAnal.Analyze_Data(DataDir);
end




parse_saved_data_for_SNRenv_mr(resultsDir,resultTxt);
if exist('ChinID','var')
    %     create_summary(ChinID);
    % %     plot_summary(ChinID);
end
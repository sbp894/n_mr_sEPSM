function save_analysis_results(PSDenv_STRUCT,PSDtfs_STRUCT,PowerMod_STRUCT,PowerTfs_STRUCT,resultsDir,resultPostfix,paramsIN)

PSDenv_S=PSDenv_STRUCT.PSDenv_S;
PSDenv_N=PSDenv_STRUCT.PSDenv_N;
PSDenv_SN=PSDenv_STRUCT.PSDenv_SN;
PSDenv_S_noisefloor=PSDenv_STRUCT.PSDenv_S_noisefloor;
PSDenv_N_noisefloor=PSDenv_STRUCT.PSDenv_N_noisefloor;
PSDenv_SN_noisefloor=PSDenv_STRUCT.PSDenv_SN_noisefloor;
PSDfreqVEC_Hz=PSDenv_STRUCT.PSDfreqVEC_Hz;

%tfs
% PSDtfs_S=PSDtfs_STRUCT.PSDtfs_S;
% PSDtfs_N=PSDtfs_STRUCT.PSDtfs_N;
% PSDtfs_SN=PSDtfs_STRUCT.PSDtfs_SN;
% PSDtfs_S_noisefloor=PSDtfs_STRUCT.PSDtfs_S_noisefloor;
% PSDtfs_N_noisefloor=PSDtfs_STRUCT.PSDtfs_N_noisefloor;
% PSDtfs_SN_noisefloor=PSDtfs_STRUCT.PSDtfs_SN_noisefloor;



PowerMod_S=PowerMod_STRUCT.PowerModS;
PowerMod_N=PowerMod_STRUCT.PowerModN;
PowerMod_SN=PowerMod_STRUCT.PowerModSN;
PowerMod_S_noisefloor=PowerMod_STRUCT.PowerModS_noisefloor;
PowerMod_N_noisefloor=PowerMod_STRUCT.PowerModN_noisefloor;
PowerMod_SN_noisefloor=PowerMod_STRUCT.PowerModSN_noisefloor;
paramsFile=PowerMod_STRUCT.params';

% PowerTfs_S=PowerTfs_STRUCT.PowerTfsS;
% PowerTfs_N=PowerTfs_STRUCT.PowerTfsN;
% PowerTfs_SN=PowerTfs_STRUCT.PowerTfsSN;
% PowerTfs_S_noisefloor=PowerTfs_STRUCT.PowerTfsS_noisefloor;
% PowerTfs_N_noisefloor=PowerTfs_STRUCT.PowerTfsN_noisefloor;
% PowerTfs_SN_noisefloor=PowerTfs_STRUCT.PowerTfsSN_noisefloor;

% Library.parsave([resultsDir 'psd' filesep 'psd' resultPostfix '.mat'],...
%     PSDenv_S,PSDenv_N,PSDenv_SN,PSDenv_S_noisefloor,PSDenv_N_noisefloor,...
%     PSDenv_SN_noisefloor,PSDtfs_S,PSDtfs_N,PSDtfs_SN,PSDtfs_S_noisefloor,...
%     PSDtfs_N_noisefloor,PSDtfs_SN_noisefloor,PSDfreqVEC_Hz, PowerMod_S,...
%     PowerMod_N, PowerMod_SN, PowerMod_S_noisefloor, PowerMod_N_noisefloor,...
%     PowerMod_SN_noisefloor, PowerTfs_S, PowerTfs_N, PowerTfs_SN, ...
%     PowerTfs_S_noisefloor, PowerTfs_N_noisefloor, PowerTfs_SN_noisefloor, ...
%     paramsFile);

% Library.parsave([resultsDir 'psd' filesep 'psd' resultPostfix '.mat'],...
%     PSDenv_S,PSDenv_N,PSDenv_SN,PSDenv_S_noisefloor,PSDenv_N_noisefloor,PSDenv_SN_noisefloor,PSDfreqVEC_Hz, ...
%     PowerMod_S, PowerMod_N, PowerMod_SN, PowerMod_S_noisefloor, PowerMod_N_noisefloor, PowerMod_SN_noisefloor,paramsFile);

Library.parsave([resultsDir 'psd' filesep 'psd' resultPostfix '.mat'],...
    PowerMod_S, PowerMod_N, PowerMod_SN, PowerMod_S_noisefloor, PowerMod_N_noisefloor, PowerMod_SN_noisefloor,paramsFile);
    
Library.parsave([resultsDir 'paramsIN' filesep 'paramsIN' resultPostfix '.mat'], paramsIN);
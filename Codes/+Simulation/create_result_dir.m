function resultsDir =create_result_dir()

%% restults path: path where the result files are stored 
resultsDir = [pwd filesep 'SCRATCH' filesep 'csche' filesep datestr(now,'yyyymmdd') filesep];
if ~exist(resultsDir,'dir')
    mkdir(resultsDir)
    mkdir([resultsDir filesep 'sac'])
    mkdir([resultsDir filesep 'psd'])
    mkdir([resultsDir filesep 'spTrains'])
    mkdir([resultsDir filesep 'sacMet'])
    mkdir([resultsDir filesep 'doc'])
    mkdir([resultsDir filesep 'doc/png'])
    mkdir([resultsDir filesep 'sepsm'])
    mkdir([resultsDir filesep 'sepsm/doc'])
    mkdir([resultsDir filesep 'sepsm/doc/png'])
end
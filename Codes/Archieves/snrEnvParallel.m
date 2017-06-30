clear;
close all
clc;
tic;

%% Define model parameters and generate stimuli
verbose = 0;
fs = 100000; % Model sampling frequency
dt = 1/fs;
% create progress dir: each iteration of the parallel loop creats a file in
% this directory. The files serve both as progress indicator and as a way
% start the simulation again at the point where it stopped in the event of
% an crash
progressDir = ['progress' filesep];
if ~exist(progressDir,'dir')
    mkdir(progressDir)
end
% restults path: path where the result files are stored
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
%%%%%%%%%%%%%%%%%%%%%%%%% Stimulus A (speech) %%%%%%%%%%%%%%%%%%%%%
A.path = ['stimuli' filesep];
A.fileExtension = '.wav';
A.prefix = ['quiet' filesep 'quiet'];
A.numberOfSentences = 10;
A.level = 50; %dB SPL according to email from Varsha
A.fs = fs;

%%%%%%%%%%%%%%%%%%%%%%%%% Stimulus B (noise) %%%%%%%%%%%%%%%%%%%%%
B.noiseTypes = {'SSN','SAM'};
B.path = ['stimuli' filesep];
B.fileExtension = '.wav';
B.noisePrefix = {['noise' filesep 'SN_Varsha_Sentences'];...
    ['noise' filesep 'SN_Varsha_Sentences_SAM']};
B.prefix = ['speechInNoise' filesep 'speechInNoise'];
B.numberOfSentences = A.numberOfSentences;
B.level = A.level; %dB SPL according to email from Varsha
B.fs = fs;
B.SNR = [-6 -3 0 3 6];

%%%%%%%%%%%%%%%%%% Auditory Nerve Model params %%%%%%%%%%%%%%%
AN.CF = [125 160 200 250 315 400 500 630 800 1000 1250 1600 2000 ...
    2500 3150 4000 5000 6300 8000]; % CF in Hz;
AN.CF = AN.CF(13); % only select every third fiber
%AN.bml = [25 25 22 22 25]; 
AN.species = 2;    % 1 for cat (2 for human with Shera et al. tuning; 3 for human with Glasberg & Moore tuning)
AN.totalLoss_dB=0;  AN.ohcLoss_dB=0;  AN.ihcLoss_dB=AN.totalLoss_dB-AN.ohcLoss_dB;
% Set OHC/IHC params for desired loss
for i = 1:length(AN.CF)
    [AN.Cohc(i),AN.Cihc(i),AN.ohcLoss(i)]=Library.fitaudiogram2(AN.CF(i),AN.totalLoss_dB,AN.species,AN.ohcLoss_dB);
end

if verbose
    cfStr=sprintf('%.2f, ',AN.CF); %#ok<*UNRCH>
    ohcStr=sprintf('%.2f, ',AN.ohcLoss);
    disp(sprintf(['IMPAIRMENT: CFs=[%s]; \nHLTOT=%.2f dB; \nHLOHC=%.2f dB' ...
        '[Actual: %s dB]; \nHLIHC=%.2f dB'],cfStr(1:end-2),...
        AN.totalLoss_dB,AN.ohcLoss_dB,ohcStr,AN.ihcLoss_dB))
end

AN.noiseType = 1;  % 1 for variable fGn (0 for fixed fGn)
AN.fiberType.sr = 2;  % spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects; "1" = Low; "2" = Medium; "3" = High
AN.fiberType.names = {'low','medium','high'}; % names of the fiber types
AN.fiberType.color = {'r','b','g'};
AN.implnt = 0;     % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse
AN.nTrials = 100;
% AN.fiberType.names{AN.fiberType.sr};

%%%%%%%%%%%%%%%%%%% SAC/SCC Parameters %%%%%%%%%%%%%%%%%%%%
anal.binWidth = 1/50000;
anal.onsetIgnore = 50e-3;
anal.maxLag = 1.5;
anal.maxSpikes = 2500;
anal.fs = fs;
anal.nTrials = AN.nTrials;
anal.plt = 1;
anal.ModFreq = [1,2,4,8,16,32,64];
anal.resultsDir = resultsDir;
anal.resultTxt = 'CF_%1.2f kHz_ Sentence_ %i_ SNR_ %i dB_ Noise_ %s_ Level_ %i dB_SPL';
anal.resultTxt1 = 'CF_%1.2f kHz_Sentence_%i_SNR_dB_Noise_%s_Level_%i_dB_SPL';
anal.resultFileName = '_%i_%i_%i_%i_%i';
anal.verbose = 0;
anal.plot = 0;

% create condition table
nCF = length( AN.CF );
nSentences = A.numberOfSentences;
nSnr = length(B.SNR);
nNoise = length(B.noiseTypes);
nLevel = length(A.level);

cs = Library.cartesianProduct(nCF,nSentences,nSnr,nNoise,nLevel);
csCell = cell(size(cs));
for i = 1 : size(cs,1)
    csCell{i,1} =  AN.CF(cs(i,1));
    csCell{i,2} =  cs(i,2);
    csCell{i,3} =  B.SNR(cs(i,3));
    csCell{i,4} =  B.noiseTypes{cs(i,4)};
    csCell{i,5} =  A.level(cs(i,5));
end

csLabels = {'CF','sentences','snr','noise','level'};
nConditions = [nCF, nSentences, nSnr, nNoise, nLevel];
save([resultsDir 'conditions.mat'],'cs','csCell','csLabels','nConditions');

%% Process stimuli with AN model
% parfor c_i = 1 : size(cs,1)
for c_i = 1 : size(cs,1)
    C = struct(); % necessary to create structure in parfor
    C.cF_i = cs(c_i,1);
    C.sentence_i = cs(c_i,2);
    C.snr_i = cs(c_i,3);
    C.noise_i = cs(c_i,4);
    C.level_i = cs(c_i,5);
    
    resultPostfix2 = sprintf(anal.resultFileName,C.cF_i,C.sentence_i,...
        C.snr_i, C.noise_i, C.level_i);
    if ~exist([progressDir filesep resultPostfix2 '.mat'],'file')
        % This function generates 3 wave files sampled at 100 kHz so that
        % rest of code doesn't have to resample each time. Also, this allows exact
        % relation of S=SN-N, which is fouled up is SN=S+N is done before resampling.
        [stim_S,stim_N,stim_SN] = Library.SNstim_resample(A,B,C,anal);
        
        [PSDenv_STRUCT,PSDtfs_STRUCT,PowerMod_STRUCT,PowerTfs_STRUCT] = ...
            Library.sumcors_mod(stim_S, stim_N, stim_SN, A, B, C, AN, anal );
        
        PSDenv_S=PSDenv_STRUCT.PSDenv_S;
        PSDenv_N=PSDenv_STRUCT.PSDenv_N;
        PSDenv_SN=PSDenv_STRUCT.PSDenv_SN;
        PSDenv_S_noisefloor=PSDenv_STRUCT.PSDenv_S_noisefloor;
        PSDenv_N_noisefloor=PSDenv_STRUCT.PSDenv_N_noisefloor;
        PSDenv_SN_noisefloor=PSDenv_STRUCT.PSDenv_SN_noisefloor;
        PSDfreqVEC_Hz=PSDenv_STRUCT.PSDfreqVEC_Hz;
        
        %tfs
        PSDtfs_S=PSDtfs_STRUCT.PSDtfs_S;
        PSDtfs_N=PSDtfs_STRUCT.PSDtfs_N;
        PSDtfs_SN=PSDtfs_STRUCT.PSDtfs_SN;
        PSDtfs_S_noisefloor=PSDtfs_STRUCT.PSDtfs_S_noisefloor;
        PSDtfs_N_noisefloor=PSDtfs_STRUCT.PSDtfs_N_noisefloor;
        PSDtfs_SN_noisefloor=PSDtfs_STRUCT.PSDtfs_SN_noisefloor;
        
        
        PowerMod_S=PowerMod_STRUCT.PowerModS;
        PowerMod_N=PowerMod_STRUCT.PowerModN;
        PowerMod_SN=PowerMod_STRUCT.PowerModSN;
        PowerMod_S_noisefloor=PowerMod_STRUCT.PowerModS_noisefloor;
        PowerMod_N_noisefloor=PowerMod_STRUCT.PowerModN_noisefloor;
        PowerMod_SN_noisefloor=PowerMod_STRUCT.PowerModSN_noisefloor;
        
        PowerTfs_S=PowerTfs_STRUCT.PowerTfsS;
        PowerTfs_N=PowerTfs_STRUCT.PowerTfsN;
        PowerTfs_SN=PowerTfs_STRUCT.PowerTfsSN;
        PowerTfs_S_noisefloor=PowerTfs_STRUCT.PowerTfsS_noisefloor;
        PowerTfs_N_noisefloor=PowerTfs_STRUCT.PowerTfsN_noisefloor;
        PowerTfs_SN_noisefloor=PowerTfs_STRUCT.PowerTfsSN_noisefloor;
        
        Library.parsave([resultsDir 'psd' filesep 'psd' resultPostfix2 '.mat'],...
            PSDenv_S,PSDenv_N,PSDenv_SN,PSDenv_S_noisefloor,PSDenv_N_noisefloor,...
            PSDenv_SN_noisefloor,PSDtfs_S,PSDtfs_N,PSDtfs_SN,PSDtfs_S_noisefloor,...
            PSDtfs_N_noisefloor,PSDtfs_SN_noisefloor,PSDfreqVEC_Hz, PowerMod_S,...
            PowerMod_N, PowerMod_SN, PowerMod_S_noisefloor, PowerMod_N_noisefloor,...
            PowerMod_SN_noisefloor, PowerTfs_S, PowerTfs_N, PowerTfs_SN, ...
            PowerTfs_S_noisefloor, PowerTfs_N_noisefloor, PowerTfs_SN_noisefloor)
        
        % update progress dir
        a = 0==0; % dummy variable to save in progress mat file
        Library.parsave([progressDir filesep resultPostfix2 '.mat'],a);
        D = dir([progressDir filesep '*.mat']);
        progress = length(D(not([D.isdir])));
        resultPostfix1 = sprintf(anal.resultTxt,AN.CF(C.cF_i)/1e3,C.sentence_i,...
            B.SNR(C.snr_i), B.noiseTypes{C.noise_i}, A.level(C.level_i) );
        disp([datestr(now,'HH:MM:SS') ' Progress: ' num2str(progress) '/'...
            num2str(size(cs,1)) ' (' resultPostfix1 ')'])
    end
end
time_taken=toc;
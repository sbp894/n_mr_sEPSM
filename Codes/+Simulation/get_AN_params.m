%%%%%%%%%%%%%%%%%% Auditory Nerve Model params %%%%%%%%%%%%%%%
function AN=get_AN_params(ExpControlParams)

verbose=0;
%%
if isempty(ExpControlParams)
    AN.CF = [1e3 ] ;%[125 160 200 250 315 400 500 630 800 1000 1250 1600 2000 2500 3150 4000 5000 6300 8000]; % CF in Hz;
    AN.fiberType.sr = 2;  % spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects; "1" = Low; "2" = Medium; "3" = High
    AN.nTrials = 25;
    AN.species = 1;    % 1 for cat (2 for human with Shera et al. tuning; 3 for human with Glasberg & Moore tuning)
    AN.ohcLoss_dB=0;
    AN.ihcLoss_dB=0;
    AN.totalLoss_dB=AN.ohcLoss_dB+ AN.ihcLoss_dB;
else
    AN.CF=ExpControlParams.CF;
    AN.fiberType.sr =ExpControlParams.fiberType;
    AN.nTrials=ExpControlParams.nRep;
    AN.species = ExpControlParams.species;    % 1 for cat (2 for human with Shera et al. tuning; 3 for human with Glasberg & Moore tuning)   
    
    AN.ohcLoss_dB=ExpControlParams.ohcLoss_dB;
    AN.ihcLoss_dB=ExpControlParams.ihcLoss_dB;
    AN.totalLoss_dB=AN.ohcLoss_dB+ AN.ihcLoss_dB;
end

%%
% AN.CF = AN.CF(1); % only select every third fiber
%AN.bml = [25 25 22 22 25];


% Set OHC/IHC params for desired loss
for i = 1:length(AN.CF)
    [AN.Cohc(i),AN.Cihc(i),AN.ohcLoss(i)]=Library.fitaudiogram2(AN.CF(i),AN.totalLoss_dB,AN.species,AN.ohcLoss_dB);
end
AN.Cohc(AN.Cohc~=1)=1;
AN.Cihc(AN.Cihc~=1)=1;

if verbose
    cfStr=sprintf('%.2f, ',AN.CF); %#ok<*UNRCH>
    ohcStr=sprintf('%.2f, ',AN.ohcLoss);
    fprintf(['IMPAIRMENT: CFs=[%s]; \nHLTOT=%.2f dB; \nHLOHC=%.2f dB' '[Actual: %s dB]; \nHLIHC=%.2f dB'],cfStr(1:end-2), AN.totalLoss_dB,AN.ohcLoss_dB,ohcStr,AN.ihcLoss_dB);
end

AN.noiseType = 0;  % 1 for variable fGn (0 for fixed fGn)
AN.fiberType.names = {'low','medium','high'}; % names of the fiber types
AN.fiberType.color = {'r','b','g'};
AN.implnt = 0;     % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse
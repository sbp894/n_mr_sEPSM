function AN=get_AN_params(verbose) %#ok<*INUSD>
%% %%%%%%%%%%%%%%%%%% Auditory Nerve Model params %%%%%%%%%%%%%%%
% AN.CF = [125 160 200 250 315 400 500 630 800 1000 1250 1600 2000 ...
%     2500 3150 4000 5000 6300 8000]; % CF in Hz;
% AN.CF = AN.CF(13); % only select every third fiber
% %AN.bml = [25 25 22 22 25]; 
% AN.species = 2;    % 1 for cat (2 for human with Shera et al. tuning; 3 for human with Glasberg & Moore tuning)
% 
% AN.totalLoss_dB=0;  
% AN.ohcLoss_dB=0;  
% AN.ihcLoss_dB=AN.totalLoss_dB-AN.ohcLoss_dB;
% 
% % Set OHC/IHC params for desired loss
% for i = 1:length(AN.CF)
%     [AN.Cohc(i),AN.Cihc(i),AN.ohcLoss(i)]=Library.fitaudiogram2(AN.CF(i),AN.totalLoss_dB,AN.species,AN.ohcLoss_dB);
% end
% 
% if verbose
%     cfStr=sprintf('%.2f, ',AN.CF); %#ok<*UNRCH>
%     ohcStr=sprintf('%.2f, ',AN.ohcLoss);
%     fprintf(['IMPAIRMENT: CFs=[%s]; \nHLTOT=%.2f dB; \nHLOHC=%.2f dB' '[Actual: %s dB]; \nHLIHC=%.2f dB'],cfStr(1:end-2), AN.totalLoss_dB,AN.ohcLoss_dB,ohcStr,AN.ihcLoss_dB);
% end
AN.Fs=100e3;
AN.noiseType = 1;  % 1 for variable fGn (0 for fixed fGn)
AN.fiberType.sr = 2;  % spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects; "1" = Low; "2" = Medium; "3" = High
AN.fiberType.names = {'low','medium','high'}; % names of the fiber types
AN.fiberType.color = {'r','b','g'};
AN.implnt = 0;     % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse
AN.nTrials = 1;

function [SACSCCfunctions,SACSCCmetrics,paramsOUT] = SACSCCanal_SNRenv(SpikeTrains,paramsIN,PLOT_15panel,resultDir,resultPostfix)
% File: CCCanal.m
% M. Heinz/ J. Swaminathan
% May 20, 2008
%
% Modified: July 29, 2008 to generalize SCadjpeak and CCCenv comps based on
% different frequency regions.
% *** MODIFIED CCCanal_3 ONLY, 0, 1, 2 have not been updated with:
%           1) max([0 NFremoval]) to imaginary numbers from sqrt(<0)
%           2) generalization to compute several frequency ranges
%
% CCCanal_0: 1 basic RUN using all spike REPs
% CCCanal_1: Use all spike REPs, run 5 RS reps, avg PSDs, compute 1 CCC
% CCCanal_2: Use all spike REPs, run 5 indep RUNS, avg 5 CCCs
% CCCanal_3: BOOTSTRAP 5 RUNS using 80% reps each, run RS rep for each,
%        3A: compute individual CCCs for each RUN, then AVG to one CCC
%        3B: AVG 5 PSDs (AN and RS), the compute 1 CCC (store as RUN: 6)
%
% General function for neural cross-correlation analyses/metrics.
%
% this program computes the envelope and fine structure coding and
% correlations from SAC/SCC/XpAC/XpCC analyses of neurophysiological spike
% train data.
%
% SpikeTrains: cell array with 4 spiketrains {condition (1,2), polarity (plus,minus)}
% paramsIN(OUT): defaults used if not included
%         .SACbinwidth_msec
%         .ignoreONSET_msec
%         .dur1_msec
%         .dur2_msec
% PLOT_15panel: 1=yes, plot fill 15-panel figure; 0=don't plot
%
% SACSCCfunctions: structure with all functions returned
%                .delays_usec
%                .SAC_A_plus   : SAC for condition A+
%                .SAC_A_minus  : SAC for condition A-
%                .SAC_A_avg    : SAC for AVG(A+,A-)
%                .SAC_B_plus   : SAC for condition B+
%                .SAC_B_minus  : SAC for condition B-
%                .SAC_B_avg    : SAC for AVG(B+,B-)
%                .SCC_C_plus  : SCC for condition C+
%                .SCC_C_minus : SCC for condition C-
%                .SCC_C_avg   : SCC for AVG(C+/C-)
%                .XpAC_A_plus  : XpAC for condition A+/A-
%                .XpAC_A_minus : XpAC for condition A-/A+
%                .XpAC_A_avg   : XpAC for AVG(A+/A-,A-/A+)
%                .XpAC_B_plus  : XpAC for condition B+/B-
%                .XpAC_B_minus : XpAC for condition B-/B+
%                .XpAC_B_avg   : XpAC for AVG(B+/B-,B-/B+)
%                .XpCC_C_plus : XpCC for condition C+/C-
%                .XpCC_C_minus: XpCC for condition C-/C+
%                .XpCC_C_avg  : XpCC for AVG(C+/C-,C-/C+)
%                .SUMCOR_A     : SUMCOR for condition A = AVG(SAC_A_avg,XpAC_A_avg)
%                .SUMCOR_B     : SUMCOR for condition B = AVG(SAC_B_avg,XpAC_B_avg)
%                .SUMCOR_C    : SUMCOR for condition 3 = AVG(SAC_C_avg,XpAC_C_avg)
%                .DIFCOR_A     : DIFCOR for condition A = SAC_A_avg-XpAC_A_avg
%                .DIFCOR_B     : DIFCOR for condition B = SAC_B_avg-XpAC_B_avg
%                .DIFCOR_C    : DIFCOR for condition C = SAC_C_avg-XpAC_C_avg
%                .PSDenv_A     : FFT of SUMCOR (PSD) for condition A = AVG(SAC_A_avg,XpAC_A_avg)
%                .PSDenv_B     : FFT of SUMCOR (PSD) for condition B = AVG(SAC_B_avg,XpAC_B_avg)
%                .PSDenv_C    : FFT of SUMCOR (PSD) for condition C = AVG(SAC_C_avg,XpAC_C_avg)
%                .PSDtfs_A     : FFT of DIFCOR (PSD) for condition A = SAC_A_avg-XpAC_A_avg
%                .PSDtfs_B     : FFT of DIFCOR (PSD) for condition B = SAC_B_avg-XpAC_B_avg
%                .PSDtfs_C    : FFT of DIFCOR (PSD) for condition C = SAC_C_avg-XpAC_C_avg
%
% SACSCCmetrics: structure with all summary metrics returned
%                .CCCenv1
%                .CCCtfs
%                .CDenv_usec
%                .CDtfs_usec
%                .CDscc_usec
%                .DCpeak
%                .SCpeak_1
%                .NumDrivenSpikes(3,2)
%                .AvgRate_sps(3,2)

if ~exist('PLOT_15panel','var'),    PLOT_15panel=0;   end
verbose = 0;
%% set the parameters for SAC and SCC %%%%%%%
paramsOUT=paramsIN;
% ignore ONSET - set default if not specified
if ~isfield(paramsOUT,'SCC_onsetIGNORE_sec')
    paramsOUT.SCC_onsetIGNORE_sec=0.05;
end
% Default SACpeak, DCpeak, SCpeak and CCCenvs to use
if ~isfield(paramsOUT,'SACpeak_TOUSE')
    paramsOUT.SACpeak_TOUSE='SACpeak_0';
end
if ~isfield(paramsOUT,'DCpeak_TOUSE')
    paramsOUT.DCpeak_TOUSE='DCpeak_0';
end
if ~isfield(paramsOUT,'CCCtfs_TOUSE')
    paramsOUT.CCCtfs_TOUSE='DCpeak_0';
end
if ~isfield(paramsOUT,'SCpeak_TOUSE')
    paramsOUT.SCpeak_TOUSE='IFFTraw_0';
end
if ~isfield(paramsOUT,'CCCenv_TOUSE')
    % 	paramsOUT.CCCenv_TOUSE='0-300, subBIAS';
    paramsOUT.CCCenv_TOUSE='IFFTrawSC_0';
end

% smoothing for functions
TriFiltWidthSAC=5; TriFiltWidthDC=1; TriFiltWidthSC=1;   % CHANGED 9/27/08
paramsOUT.TriFiltWidthSAC=TriFiltWidthSAC;
paramsOUT.TriFiltWidthDC=TriFiltWidthDC;
paramsOUT.TriFiltWidthSC=TriFiltWidthSC;
% SCC binwidth (Joris, 2003)
paramsOUT.DELAYbinwidth_sec=50e-6;
% FFT window length to get desired (was 1-Hz) sampling in frequency domain
% Nfft_psd=round(1/paramsOUT.DELAYbinwidth_sec);  % 20000
%% *MH July 7 2015 - change to get <0.2 Hz resolution
Nfft_psd=2^nextpow2(5*round(1/paramsOUT.DELAYbinwidth_sec));  % 5*20000 - Next power of 2
paramsOUT.Nfft_psd=Nfft_psd;
% SAC/SCC window length
if ~isfield(paramsOUT,'MAXdelay_sec')
    paramsOUT.MAXdelay_sec=1;
end
% stim duration for SAC analysis = min of all conditions
stimdur_sec=min([paramsIN.durA_msec paramsIN.durB_msec paramsIN.durC_msec])/1000;
paramsOUT.stimdur_sec=stimdur_sec;
% Limit the number of spikes in the SAC/SCC analysis
if ~isfield(paramsOUT,'MAXspikes')
    paramsOUT.MAXspikes=3700;  % If used, Library.windowSTs will cut out extra REPs beyond MAXspikes - used to avoid memory limits
end
% minimum DCpeak for A and B to compute CCCtfs
paramsOUT.minDCpeak_CCCtfs=0.1;
% CF: use minimum of all CFs to be most conservative for adjSCpeak
CF_Hz=min([paramsOUT.CF_A_Hz paramsOUT.CF_B_Hz paramsOUT.CF_C_Hz]);
paramsOUT.SACSCC_CF_Hz=CF_Hz;
%% BOOTSTRAPPING params
paramsOUT.BOOTSTRAP_percentdata=0.80;  % Percent of AN reps to include in each RUN

%% BOOTSTRAPPING - SETUP spike REPS to use for each RUN
clear Nreps
Nreps=zeros(3,2);
Ninclude=zeros(3,2);
Nexclude=zeros(3,2);
Navgs=zeros(3,2);

for i=1:3
    for j=1:2
        Nreps(i,j)=length(SpikeTrains{i,j});
        Ninclude(i,j)=ceil(paramsOUT.BOOTSTRAP_percentdata*Nreps(i,j)); % some rounding issues with ceil/floor - have to do it this way
        Nexclude(i,j)=max([1 Nreps(i,j)-Ninclude(i,j)]);  % at least 1 rep excluded to make sure runs are different
        Navgs(i,j)=floor(Nreps(i,j)/Nexclude(i,j));  % Number of RUNS for Averaging CCCs, or PSDs
    end
end

paramsOUT.BOOTSTRAP_Navgs=min(min(Navgs));  % Number of RUNS for Averaging CCCs, or PSDs
BOOTinds=cell(3,2);
for i=1:3
    for j=1:2
        BOOTinds{i,j}=cell(1,paramsOUT.BOOTSTRAP_Navgs);
    end
end

for BOOTrep1=1:paramsOUT.BOOTSTRAP_Navgs
    for i=1:3
        for j=1:2
            TEMPinds=setdiff((1:Nreps(i,j)),(Nreps(i,j)-Nexclude(i,j)+1:Nreps(i,j))-(BOOTrep1-1)*Nexclude(i,j));
            %% NEED TO RANDOMIZE SPIKES to avoid replication if > 5000 spikes
            %% in 50 reps (e.g., if only 20 reps are needed to get 5000 spikes, then
            %% 1st 3 bootstrap samples will all use 1:20
            randINDs=randperm(length(TEMPinds));
            BOOTinds{i,j}{BOOTrep1}=TEMPinds(randINDs);
        end
    end
end

%% First paramsOUT.BOOTSTRAP_Navgs reps are boot strap reps using 80% of
%% reps, last BOOTrep uses the AVG of all reps for all computations
SACSCCfunctions=cell(paramsOUT.BOOTSTRAP_Navgs+1,1);
SACSCCmetrics=cell(paramsOUT.BOOTSTRAP_Navgs+1,1);

for BOOTrep=1:paramsOUT.BOOTSTRAP_Navgs+1
    % Run BOOTSTRAP_Navgs different random samples
    if BOOTrep~=paramsOUT.BOOTSTRAP_Navgs+1
        %% Window spike times (ignore 1st 50 ms, and cut to shorter of 2 stim)
        [ST_A_plus,NumDrivenSpikes(1,1)]=Library.windowSTs(SpikeTrains{1,1}(BOOTinds{1,1}{BOOTrep}), ...
            paramsOUT.SCC_onsetIGNORE_sec,stimdur_sec,paramsOUT.MAXspikes);
        [ST_A_minus,NumDrivenSpikes(1,2)]=Library.windowSTs(SpikeTrains{1,2}(BOOTinds{1,2}{BOOTrep}), ...
            paramsOUT.SCC_onsetIGNORE_sec,stimdur_sec,paramsOUT.MAXspikes);
        [ST_B_plus,NumDrivenSpikes(2,1)]=Library.windowSTs(SpikeTrains{2,1}(BOOTinds{2,1}{BOOTrep}), ...
            paramsOUT.SCC_onsetIGNORE_sec,stimdur_sec,paramsOUT.MAXspikes);
        [ST_B_minus,NumDrivenSpikes(2,2)]=Library.windowSTs(SpikeTrains{2,2}(BOOTinds{2,2}{BOOTrep}), ...
            paramsOUT.SCC_onsetIGNORE_sec,stimdur_sec,paramsOUT.MAXspikes);
        [ST_C_plus,NumDrivenSpikes(3,1)]=Library.windowSTs(SpikeTrains{3,1}(BOOTinds{3,1}{BOOTrep}), ...
            paramsOUT.SCC_onsetIGNORE_sec,stimdur_sec,paramsOUT.MAXspikes);
        [ST_C_minus,NumDrivenSpikes(3,2)]=Library.windowSTs(SpikeTrains{3,2}(BOOTinds{3,2}{BOOTrep}), ...
            paramsOUT.SCC_onsetIGNORE_sec,stimdur_sec,paramsOUT.MAXspikes);
        paramsOUT.SCCdur_sec=stimdur_sec-paramsOUT.SCC_onsetIGNORE_sec;
        
        NumReps(1,1)=length(ST_A_plus); NumReps(1,2)=length(ST_A_minus);
        NumReps(2,1)=length(ST_B_plus); NumReps(2,2)=length(ST_B_minus);
        NumReps(3,1)=length(ST_C_plus); NumReps(3,2)=length(ST_C_minus);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if verbose
            fprintf('NUMBER OF SPIKES (A_plus, A_minus, B_plus, B_minus, C_plus, C_minus) = %d, %d, %d, %d, %d, %d',NumDrivenSpikes(1,1),NumDrivenSpikes(1,2), ...
                NumDrivenSpikes(2,1),NumDrivenSpikes(2,2),NumDrivenSpikes(3,1),NumDrivenSpikes(3,2));
            fprintf('NUMBER OF REPS (A_plus, A_minus, B_plus, B_minus, C_plus, C_minus) = %d, %d, %d, %d, %d, %d',NumReps(1,1),NumReps(1,2),NumReps(2,1),NumReps(2,2),NumReps(3,1),NumReps(3,2));
            fprintf('   ... Computing SACs/SCCs for AN spikes: REP %d of %d ...',BOOTrep,paramsOUT.BOOTSTRAP_Navgs);
        end
        [SACSCCs,AvgRate_sps]=SACSCCanal(ST_A_plus,ST_A_minus,ST_B_plus,ST_B_minus,ST_C_plus,ST_C_minus,paramsOUT);
        %%% OK to HERE - going NOW to SACSCCanal
        
        % 		if ~isnan(NumDrivenSpikes(2,1))
        %% RANDOMIZE spike times within same window used above (ignore 1st 50 ms, and cut to shorter of 2 stim)
        STrand_A_plus=Library.randomizeSTs(ST_A_plus,paramsOUT.SCC_onsetIGNORE_sec,stimdur_sec);
        STrand_A_minus=Library.randomizeSTs(ST_A_minus,paramsOUT.SCC_onsetIGNORE_sec,stimdur_sec);
        STrand_B_plus=Library.randomizeSTs(ST_B_plus,paramsOUT.SCC_onsetIGNORE_sec,stimdur_sec);
        STrand_B_minus=Library.randomizeSTs(ST_B_minus,paramsOUT.SCC_onsetIGNORE_sec,stimdur_sec);
        STrand_C_plus=Library.randomizeSTs(ST_C_plus,paramsOUT.SCC_onsetIGNORE_sec,stimdur_sec);
        STrand_C_minus=Library.randomizeSTs(ST_C_minus,paramsOUT.SCC_onsetIGNORE_sec,stimdur_sec);
        if verbose
            fprintf('   ... Computing SACs/SCCs for RANDOM spikes (for Noise Floor) ...'); %#ok<*UNRCH>
        end
        [SACSCCs_rand,~]=SACSCCanal(STrand_A_plus,STrand_A_minus,STrand_B_plus,STrand_B_minus,STrand_C_plus,STrand_C_minus,paramsOUT);
        % 		else
        % 			SACSCCs_rand=[];
        % 		end
        
    else	% For LAST loop index, take AVG of all functions for computations
        % AVG everything! Then run rest of code ASIS!!
        % currently it has the last REP already
        for BOOTrep2=1:paramsOUT.BOOTSTRAP_Navgs-1
            SACSCCs.SAC_A_avg=SACSCCs.SAC_A_avg+SACSCCfunctions{BOOTrep2}.SAC_A_avg;
            SACSCCs.SAC_B_avg=SACSCCs.SAC_B_avg+SACSCCfunctions{BOOTrep2}.SAC_B_avg;
            SACSCCs.SAC_C_avg=SACSCCs.SAC_C_avg+SACSCCfunctions{BOOTrep2}.SAC_C_avg;
            SACSCCs.SCC_AC_avg=SACSCCs.SCC_AC_avg+SACSCCfunctions{BOOTrep2}.SCC_AC_avg;
            SACSCCs.XpAC_A_avg=SACSCCs.XpAC_A_avg+SACSCCfunctions{BOOTrep2}.XpAC_A_avg;
            SACSCCs.XpAC_B_avg=SACSCCs.XpAC_B_avg+SACSCCfunctions{BOOTrep2}.XpAC_B_avg;
            SACSCCs.XpAC_C_avg=SACSCCs.XpAC_C_avg+SACSCCfunctions{BOOTrep2}.XpAC_C_avg;
            SACSCCs.XpCC_AC_avg=SACSCCs.XpCC_AC_avg+SACSCCfunctions{BOOTrep2}.XpCC_AC_avg;
            SACSCCs.SUMCOR_A=SACSCCs.SUMCOR_A+SACSCCfunctions{BOOTrep2}.SUMCOR_A;
            SACSCCs.SUMCOR_B=SACSCCs.SUMCOR_B+SACSCCfunctions{BOOTrep2}.SUMCOR_B;
            SACSCCs.SUMCOR_C=SACSCCs.SUMCOR_C+SACSCCfunctions{BOOTrep2}.SUMCOR_C;
            SACSCCs.SUMCOR_AC=SACSCCs.SUMCOR_AC+SACSCCfunctions{BOOTrep2}.SUMCOR_AC;
            SACSCCs.SUMCORadj_A=SACSCCs.SUMCORadj_A+SACSCCfunctions{BOOTrep2}.SUMCORadj_A;
            SACSCCs.SUMCORadj_B=SACSCCs.SUMCORadj_B+SACSCCfunctions{BOOTrep2}.SUMCORadj_B;
            SACSCCs.SUMCORadj_C=SACSCCs.SUMCORadj_C+SACSCCfunctions{BOOTrep2}.SUMCORadj_C;
            SACSCCs.SUMCORadj_AC=SACSCCs.SUMCORadj_AC+SACSCCfunctions{BOOTrep2}.SUMCORadj_AC;
            SACSCCs.DIFCOR_A=SACSCCs.DIFCOR_A+SACSCCfunctions{BOOTrep2}.DIFCOR_A;
            SACSCCs.DIFCOR_B=SACSCCs.DIFCOR_B+SACSCCfunctions{BOOTrep2}.DIFCOR_B;
            SACSCCs.DIFCOR_C=SACSCCs.DIFCOR_C+SACSCCfunctions{BOOTrep2}.DIFCOR_C;
            SACSCCs.DIFCOR_AC=SACSCCs.DIFCOR_AC+SACSCCfunctions{BOOTrep2}.DIFCOR_AC;
            %PSDenv
            SACSCCs.PSDenv_A=SACSCCs.PSDenv_A+SACSCCfunctions{BOOTrep2}.PSDenv_A;
            SACSCCs.PSDenv_B=SACSCCs.PSDenv_B+SACSCCfunctions{BOOTrep2}.PSDenv_B;
            SACSCCs.PSDenv_C=SACSCCs.PSDenv_C+SACSCCfunctions{BOOTrep2}.PSDenv_C;
            SACSCCs.PSDenv_AC=SACSCCs.PSDenv_AC+SACSCCfunctions{BOOTrep2}.PSDenv_AC;
            % RAND
            SACSCCs_rand.PSDenv_A=SACSCCs_rand.PSDenv_A+SACSCCfunctions{BOOTrep2}.rand.PSDenv_A;
            SACSCCs_rand.PSDenv_B=SACSCCs_rand.PSDenv_B+SACSCCfunctions{BOOTrep2}.rand.PSDenv_B;
            SACSCCs_rand.PSDenv_C=SACSCCs_rand.PSDenv_C+SACSCCfunctions{BOOTrep2}.rand.PSDenv_C;
            SACSCCs_rand.PSDenv_AC=SACSCCs_rand.PSDenv_AC+SACSCCfunctions{BOOTrep2}.rand.PSDenv_AC;
            %PSDtfs
            SACSCCs.PSDtfs_A=SACSCCs.PSDtfs_A+SACSCCfunctions{BOOTrep2}.PSDtfs_A;
            SACSCCs.PSDtfs_B=SACSCCs.PSDtfs_B+SACSCCfunctions{BOOTrep2}.PSDtfs_B;
            SACSCCs.PSDtfs_C=SACSCCs.PSDtfs_C+SACSCCfunctions{BOOTrep2}.PSDtfs_C;
            SACSCCs.PSDtfs_AC=SACSCCs.PSDtfs_AC+SACSCCfunctions{BOOTrep2}.PSDtfs_AC;
            % RAND
            SACSCCs_rand.PSDtfs_A=SACSCCs_rand.PSDtfs_A+SACSCCfunctions{BOOTrep2}.rand.PSDtfs_A;
            SACSCCs_rand.PSDtfs_B=SACSCCs_rand.PSDtfs_B+SACSCCfunctions{BOOTrep2}.rand.PSDtfs_B;
            SACSCCs_rand.PSDtfs_C=SACSCCs_rand.PSDtfs_C+SACSCCfunctions{BOOTrep2}.rand.PSDtfs_C;
            SACSCCs_rand.PSDtfs_AC=SACSCCs_rand.PSDtfs_AC+SACSCCfunctions{BOOTrep2}.rand.PSDtfs_AC;
            
        end
        SACSCCs.SAC_A_avg=SACSCCs.SAC_A_avg/paramsOUT.BOOTSTRAP_Navgs;
        SACSCCs.SAC_B_avg=SACSCCs.SAC_B_avg/paramsOUT.BOOTSTRAP_Navgs;
        SACSCCs.SAC_C_avg=SACSCCs.SAC_C_avg/paramsOUT.BOOTSTRAP_Navgs;
        SACSCCs.SCC_AC_avg=SACSCCs.SCC_AC_avg/paramsOUT.BOOTSTRAP_Navgs;
        SACSCCs.XpAC_A_avg=SACSCCs.XpAC_A_avg/paramsOUT.BOOTSTRAP_Navgs;
        SACSCCs.XpAC_B_avg=SACSCCs.XpAC_B_avg/paramsOUT.BOOTSTRAP_Navgs;
        SACSCCs.XpAC_C_avg=SACSCCs.XpAC_C_avg/paramsOUT.BOOTSTRAP_Navgs;
        SACSCCs.XpCC_AC_avg=SACSCCs.XpCC_AC_avg/paramsOUT.BOOTSTRAP_Navgs;
        SACSCCs.SUMCOR_A=SACSCCs.SUMCOR_A/paramsOUT.BOOTSTRAP_Navgs;
        SACSCCs.SUMCOR_B=SACSCCs.SUMCOR_B/paramsOUT.BOOTSTRAP_Navgs;
        SACSCCs.SUMCOR_C=SACSCCs.SUMCOR_C/paramsOUT.BOOTSTRAP_Navgs;
        SACSCCs.SUMCOR_AC=SACSCCs.SUMCOR_AC/paramsOUT.BOOTSTRAP_Navgs;
        SACSCCs.SUMCORadj_A=SACSCCs.SUMCORadj_A/paramsOUT.BOOTSTRAP_Navgs;
        SACSCCs.SUMCORadj_B=SACSCCs.SUMCORadj_B/paramsOUT.BOOTSTRAP_Navgs;
        SACSCCs.SUMCORadj_C=SACSCCs.SUMCORadj_C/paramsOUT.BOOTSTRAP_Navgs;
        SACSCCs.SUMCORadj_AC=SACSCCs.SUMCORadj_AC/paramsOUT.BOOTSTRAP_Navgs;
        SACSCCs.DIFCOR_A=SACSCCs.DIFCOR_A/paramsOUT.BOOTSTRAP_Navgs;
        SACSCCs.DIFCOR_B=SACSCCs.DIFCOR_B/paramsOUT.BOOTSTRAP_Navgs;
        SACSCCs.DIFCOR_C=SACSCCs.DIFCOR_C/paramsOUT.BOOTSTRAP_Navgs;
        SACSCCs.DIFCOR_AC=SACSCCs.DIFCOR_AC/paramsOUT.BOOTSTRAP_Navgs;
        %%%%%%%%%%%%%%%%%%%%%%%
        % PSDs - Best to AVG in PSD domain, not in time domain to reduce variance (e.g., Welch technique)
        SACSCCs.PSDenv_A=SACSCCs.PSDenv_A/paramsOUT.BOOTSTRAP_Navgs;
        SACSCCs.PSDenv_B=SACSCCs.PSDenv_B/paramsOUT.BOOTSTRAP_Navgs;
        SACSCCs.PSDenv_C=SACSCCs.PSDenv_C/paramsOUT.BOOTSTRAP_Navgs;
        SACSCCs.PSDenv_AC=SACSCCs.PSDenv_AC/paramsOUT.BOOTSTRAP_Navgs;
        SACSCCs_rand.PSDenv_A=SACSCCs_rand.PSDenv_A/paramsOUT.BOOTSTRAP_Navgs;
        SACSCCs_rand.PSDenv_B=SACSCCs_rand.PSDenv_B/paramsOUT.BOOTSTRAP_Navgs;
        SACSCCs_rand.PSDenv_C=SACSCCs_rand.PSDenv_C/paramsOUT.BOOTSTRAP_Navgs;
        SACSCCs_rand.PSDenv_AC=SACSCCs_rand.PSDenv_AC/paramsOUT.BOOTSTRAP_Navgs;
        SACSCCs.PSDtfs_A=SACSCCs.PSDtfs_A/paramsOUT.BOOTSTRAP_Navgs;
        SACSCCs.PSDtfs_B=SACSCCs.PSDtfs_B/paramsOUT.BOOTSTRAP_Navgs;
        SACSCCs.PSDtfs_C=SACSCCs.PSDtfs_C/paramsOUT.BOOTSTRAP_Navgs;
        SACSCCs.PSDtfs_AC=SACSCCs.PSDtfs_AC/paramsOUT.BOOTSTRAP_Navgs;
        SACSCCs_rand.PSDtfs_A=SACSCCs_rand.PSDtfs_A/paramsOUT.BOOTSTRAP_Navgs;
        SACSCCs_rand.PSDtfs_B=SACSCCs_rand.PSDtfs_B/paramsOUT.BOOTSTRAP_Navgs;
        SACSCCs_rand.PSDtfs_C=SACSCCs_rand.PSDtfs_C/paramsOUT.BOOTSTRAP_Navgs;
        SACSCCs_rand.PSDtfs_AC=SACSCCs_rand.PSDtfs_AC/paramsOUT.BOOTSTRAP_Navgs;
        % 		%%%%%%%%%%%%%%%%%%%%%%%
        % 		%%% Recompute PSDs/CSD based on AVG SCs
        % 		%%%%%%%%%%%%%%%%%%%%%%%
        % 		SACSCCs.PSDenv_A=abs(fft((SACSCCs.SUMCOR_A-1),paramsOUT.Nfft_psd));
        % 		SACSCCs.PSDenv_B=abs(fft((SACSCCs.SUMCOR_B-1),paramsOUT.Nfft_psd));
        % 		SACSCCs.PSDenv_C=abs(fft((SACSCCs.SUMCOR_AC-1),paramsOUT.Nfft_psd));
        
        NumDrivenSpikes=zeros(3,2);
        AvgRate_sps=zeros(3,2);
        for BOOTrep3=1:paramsOUT.BOOTSTRAP_Navgs
            NumDrivenSpikes=NumDrivenSpikes+SACSCCmetrics{BOOTrep3}.NumDrivenSpikes/paramsOUT.BOOTSTRAP_Navgs;
            AvgRate_sps=AvgRate_sps+SACSCCmetrics{BOOTrep3}.AvgRate_sps/paramsOUT.BOOTSTRAP_Navgs;
        end
    end
    
    
    %% PSD/CSD summations - Compute summed energy in ENVELOPE spectral (densities over various frequency ranges)
    
    % create complete list of LH frequencies over which sums are needed for
    % SCpeak_adjusted and/or CCCenv
    NYQ_Hz=0.5*(1/paramsOUT.DELAYbinwidth_sec);
    % default list
    default_PSD_LHfreqs_Hz=[0 NYQ_Hz; 0 CF_Hz; 0 300; 10 300; 10 CF_Hz; 0 50; 0 100; 0 64; 0 150];
    % add user-passed list, and create unique list
    if isfield(paramsOUT,'PSD_LHfreqs_Hz')
        paramsOUT.PSD_LHfreqs_Hz=unique([default_PSD_LHfreqs_Hz; paramsOUT.PSD_LHfreqs_Hz],'rows');
    else
        paramsOUT.PSD_LHfreqs_Hz=default_PSD_LHfreqs_Hz;
    end
    
    % find INDs in freqVEC_Hz for relevant cutoffs
    PSD_LHfreqs_inds=zeros(size(paramsOUT.PSD_LHfreqs_Hz,1),2);
    for i=1:size(paramsOUT.PSD_LHfreqs_Hz,1)
        [~,PSD_LHfreqs_inds(i,1)]=min(abs(SACSCCs.freqVEC-paramsOUT.PSD_LHfreqs_Hz(i,1)));
        [~,PSD_LHfreqs_inds(i,2)]=min(abs(SACSCCs.freqVEC-paramsOUT.PSD_LHfreqs_Hz(i,2)));
    end
    
    % Compute all sums for all spectral densities
    for i=1:size(paramsOUT.PSD_LHfreqs_Hz,1)
        SACSCCmetrics{BOOTrep}.sums.sumPSD_A(i) = sum(SACSCCs.PSDenv_A(PSD_LHfreqs_inds(i,1):PSD_LHfreqs_inds(i,2)));
        if ~isnan(NumDrivenSpikes(2,1))
            SACSCCmetrics{BOOTrep}.sums.sumPSD_B(i) = sum(SACSCCs.PSDenv_B(PSD_LHfreqs_inds(i,1):PSD_LHfreqs_inds(i,2)));
            SACSCCmetrics{BOOTrep}.sums.sumCSD_AC(i) = sum(SACSCCs.PSDenv_C(PSD_LHfreqs_inds(i,1):PSD_LHfreqs_inds(i,2)));
            SACSCCmetrics{BOOTrep}.sums.sumPSDrand_A(i) = sum(SACSCCs_rand.PSDenv_A(PSD_LHfreqs_inds(i,1):PSD_LHfreqs_inds(i,2)));
            SACSCCmetrics{BOOTrep}.sums.sumPSDrand_B(i) = sum(SACSCCs_rand.PSDenv_B(PSD_LHfreqs_inds(i,1):PSD_LHfreqs_inds(i,2)));
            SACSCCmetrics{BOOTrep}.sums.sumCSDrand_AC(i) = sum(SACSCCs_rand.PSDenv_C(PSD_LHfreqs_inds(i,1):PSD_LHfreqs_inds(i,2)));
        end
    end
    SACSCCmetrics{BOOTrep}.sums.PSD_LHfreqs_Hz=paramsOUT.PSD_LHfreqs_Hz;
    
    % TFS
    NYQ_Hz=0.5*(1/paramsOUT.DELAYbinwidth_sec);
    % default list
%     default_PSD_LHfreqs_Hz=[0 NYQ_Hz];
    % find INDs in freqVEC_Hz for relevant cutoffs
    for i=1:size(paramsOUT.PSD_LHfreqs_Hz,1)
        [~,PSD_LHfreqs_inds(i,1)]=min(abs(SACSCCs.freqVEC-paramsOUT.PSD_LHfreqs_Hz(i,1)));
        [~,PSD_LHfreqs_inds(i,2)]=min(abs(SACSCCs.freqVEC-paramsOUT.PSD_LHfreqs_Hz(i,2)));
    end
    % Compute all sums for all spectral densities
    for i=1:size(paramsOUT.PSD_LHfreqs_Hz,1)
        SACSCCmetrics{BOOTrep}.sums.sumPSDtfs_A(i) = sum(SACSCCs.PSDtfs_A(PSD_LHfreqs_inds(i,1):PSD_LHfreqs_inds(i,2)));
        if ~isnan(NumDrivenSpikes(2,1))
            SACSCCmetrics{BOOTrep}.sums.sumPSDtfs_B(i) = sum(SACSCCs.PSDtfs_B(PSD_LHfreqs_inds(i,1):PSD_LHfreqs_inds(i,2)));
            SACSCCmetrics{BOOTrep}.sums.sumCSDtfs_C(i) = sum(SACSCCs.PSDtfs_C(PSD_LHfreqs_inds(i,1):PSD_LHfreqs_inds(i,2)));
            SACSCCmetrics{BOOTrep}.sums.sumPSDrandTfs_A(i) = sum(SACSCCs_rand.PSDtfs_A(PSD_LHfreqs_inds(i,1):PSD_LHfreqs_inds(i,2)));
            SACSCCmetrics{BOOTrep}.sums.sumPSDrandTfs_B(i) = sum(SACSCCs_rand.PSDtfs_B(PSD_LHfreqs_inds(i,1):PSD_LHfreqs_inds(i,2)));
            SACSCCmetrics{BOOTrep}.sums.sumCSDrandTfs_AC(i) = sum(SACSCCs_rand.PSDtfs_C(PSD_LHfreqs_inds(i,1):PSD_LHfreqs_inds(i,2)));
        end
    end
    SACSCCmetrics{BOOTrep}.sums.PSD_LHfreqsTfs_Hz=paramsOUT.PSD_LHfreqs_Hz;
    
    
    %%%%%%%%%%%%%%%%%%%%
    %% Compute metrics
    %%%%%%%%%%%%%%%%%%%%
    
    %% Characteristic Delays
    %% NEED TO IMPROVE - find largest, unless there is a "close-to-largest" that
    %% is closer to zero!
    % probably use some criterion for 2nd largest peak (if with 5%)???
    SACSCCmetrics{BOOTrep}.CDscc_usec=Library.findCD_SCC(SACSCCs.SCC_AC_avg,SACSCCs.delays_usec);
    SACSCCmetrics{BOOTrep}.CDenv_usec=Library.findCD_SCC(SACSCCs.SUMCORadj_AC,SACSCCs.delays_usec);  % use IFFTadjusted SUMCOR!
    SACSCCmetrics{BOOTrep}.CDtfs_usec=Library.findCD_SCC(SACSCCs.DIFCOR_AC,SACSCCs.delays_usec);
    
    
    %% SAC/DC/SC Peak Heights
    %%%%%%%%%
    %% SAC peaks - this is CI from Joris et al 2006 (HR) [don't subtract 1]
    % 1) SACpeak (pure MAX)
    SACSCCmetrics{BOOTrep}.SACpeaks_legend{1}='SACpeak_max';
    SACSCCmetrics{BOOTrep}.SACpeak_A(1)=max(SACSCCs.SAC_A_avg);
    SACSCCmetrics{BOOTrep}.SACpeak_B(1)=max(SACSCCs.SAC_B_avg);
    SACSCCmetrics{BOOTrep}.SACpeak_C(1)=max(SACSCCs.SAC_C_avg);
    % 2) SACpeak (0 delay)
    SACSCCmetrics{BOOTrep}.SACpeaks_legend{2}='SACpeak_0';
    INDEX_0=find(SACSCCs.delays_usec==0);
    SACSCCmetrics{BOOTrep}.SACpeak_A(2)=SACSCCs.SAC_A_avg(INDEX_0);
    if ~isnan(NumDrivenSpikes(2,1))
        SACSCCmetrics{BOOTrep}.SACpeak_B(2)=SACSCCs.SAC_B_avg(INDEX_0);
        SACSCCmetrics{BOOTrep}.SACpeak_C(2)=SACSCCs.SAC_C_avg(INDEX_0);
    end
    % 3) SACpeak (CD delay)
    SACSCCmetrics{BOOTrep}.SACpeaks_legend{3}='SACpeak_CD';
    INDEX_CD=find(SACSCCs.delays_usec==SACSCCmetrics{BOOTrep}.CDscc_usec);
    SACSCCmetrics{BOOTrep}.SACpeak_A(3)=SACSCCs.SAC_A_avg(INDEX_CD);
    if ~isnan(NumDrivenSpikes(2,1))
        SACSCCmetrics{BOOTrep}.SACpeak_B(3)=SACSCCs.SAC_B_avg(INDEX_CD);
        SACSCCmetrics{BOOTrep}.SCCpeak_AC(3)=SACSCCs.SCC_AC_avg(INDEX_CD);
    end
    %% User-passed list of delays to compute
    if isfield(paramsOUT,'UserDelays_usec')
        for i=1:length(paramsOUT.UserDelays_usec)
            % 4+) SACpeak (user delays)
            SACSCCmetrics{BOOTrep}.SACpeaks_legend{end+1}=sprintf('SACpeak_%d',round(paramsOUT.UserDelays_usec(i)));
            [~,INDEX_user]=min(abs(SACSCCs.delays_usec-paramsOUT.UserDelays_usec(i))); % find closest delay
            SACSCCmetrics{BOOTrep}.SACpeak_A(end+1)=SACSCCs.SAC_A_avg(INDEX_user);
            if ~isnan(NumDrivenSpikes(2,1))
                SACSCCmetrics{BOOTrep}.SACpeak_B(end+1)=SACSCCs.SAC_B_avg(INDEX_user);
                SACSCCmetrics{BOOTrep}.SCCpeak_AC(end+1)=SACSCCs.SCC_AC_avg(INDEX_user);
            end
        end
    end
    
    %%%%%%%%%
    %% DC peaks -
    % 1) DCpeak (pure MAX)
    SACSCCmetrics{BOOTrep}.DCpeaks_legend{1}='DCpeak_max';
    SACSCCmetrics{BOOTrep}.DCpeak_A(1)=max(SACSCCs.DIFCOR_A);
    if ~isnan(NumDrivenSpikes(2,1))
        SACSCCmetrics{BOOTrep}.DCpeak_B(1)=max(SACSCCs.DIFCOR_B);
        SACSCCmetrics{BOOTrep}.DCpeak_C(1)=max(SACSCCs.DIFCOR_C);
    end
    % 2) DCpeak (0 delay)
    SACSCCmetrics{BOOTrep}.DCpeaks_legend{2}='DCpeak_0';
    SACSCCmetrics{BOOTrep}.DCpeak_A(2)=SACSCCs.DIFCOR_A(INDEX_0);
    if ~isnan(NumDrivenSpikes(2,1))
        SACSCCmetrics{BOOTrep}.DCpeak_B(2)=SACSCCs.DIFCOR_B(INDEX_0);
        SACSCCmetrics{BOOTrep}.DCpeak_C(2)=SACSCCs.DIFCOR_C(INDEX_0);
        SACSCCmetrics{BOOTrep}.DCpeak_AC(2)=SACSCCs.DIFCOR_AC(INDEX_0);
    end
    % 3) DCpeak (CD delay)
    SACSCCmetrics{BOOTrep}.DCpeaks_legend{3}='DCpeak_CD';
    INDEX_CD=find(SACSCCs.delays_usec==SACSCCmetrics{BOOTrep}.CDtfs_usec);
    SACSCCmetrics{BOOTrep}.DCpeak_A(3)=SACSCCs.DIFCOR_A(INDEX_CD);
    if ~isnan(NumDrivenSpikes(2,1))
        SACSCCmetrics{BOOTrep}.DCpeak_B(3)=SACSCCs.DIFCOR_B(INDEX_CD);
        SACSCCmetrics{BOOTrep}.DCpeak_AC(3)=SACSCCs.DIFCOR_AC(INDEX_CD);
    end
    %% User-passed list of delays to compute
    if isfield(paramsOUT,'UserDelays_usec')
        for i=1:length(paramsOUT.UserDelays_usec)
            % 4+) DCpeak (user delays)
            SACSCCmetrics{BOOTrep}.DCpeaks_legend{end+1}=sprintf('DCpeak_%d',round(paramsOUT.UserDelays_usec(i)));
            [~,INDEX_user]=min(abs(SACSCCs.delays_usec-paramsOUT.UserDelays_usec(i))); % find closest delay
            SACSCCmetrics{BOOTrep}.DCpeak_A(end+1)=SACSCCs.DIFCOR_A(INDEX_user);
            if ~isnan(NumDrivenSpikes(2,1))
                SACSCCmetrics{BOOTrep}.DCpeak_B(end+1)=SACSCCs.DIFCOR_B(INDEX_user);
                SACSCCmetrics{BOOTrep}.DCpeak_AC(end+1)=SACSCCs.DIFCOR_AC(INDEX_user);
            end
        end
    end
    
    %%%%%%%%%
    %% SUMCOR peaks (don't subtract 1; Louage et al 2004)
    % 1) raw peaks
    SACSCCmetrics{BOOTrep}.SCpeaks_legend{1}='raw';
    SACSCCmetrics{BOOTrep}.SCpeaks_A(1)=max(SACSCCs.SUMCOR_A);
    if ~isnan(NumDrivenSpikes(2,1))
        SACSCCmetrics{BOOTrep}.SCpeaks_B(1)=max(SACSCCs.SUMCOR_B);
        SACSCCmetrics{BOOTrep}.SCpeaks_C(1)=max(SACSCCs.SUMCOR_C);
    end
    % 2) Adjusted SCpeak (Eq. 2): Adj SCpeak = rawSCpeak*sum[0,CF]/sum[0,Fs/2]
    SACSCCmetrics{BOOTrep}.SCpeaks_legend{2}='adj: 0-CF';
    CFsum_index=find((paramsOUT.PSD_LHfreqs_Hz(:,1)==0)&(paramsOUT.PSD_LHfreqs_Hz(:,2)==CF_Hz));
    NYQsum_index=find((paramsOUT.PSD_LHfreqs_Hz(:,1)==0)&(paramsOUT.PSD_LHfreqs_Hz(:,2)==NYQ_Hz));
    SACSCCmetrics{BOOTrep}.SCpeaks_A(2)=1+(max(SACSCCs.SUMCOR_A)-1)* ...
        SACSCCmetrics{BOOTrep}.sums.sumPSD_A(CFsum_index)/SACSCCmetrics{BOOTrep}.sums.sumPSD_A(NYQsum_index);
    if ~isnan(NumDrivenSpikes(2,1))
        SACSCCmetrics{BOOTrep}.SCpeaks_B(2)=1+(max(SACSCCs.SUMCOR_B)-1)* ...
            SACSCCmetrics{BOOTrep}.sums.sumPSD_B(CFsum_index)/SACSCCmetrics{BOOTrep}.sums.sumPSD_B(NYQsum_index);
        SACSCCmetrics{BOOTrep}.SCpeaks_AC(2)=1+(max(SACSCCs.SUMCOR_AC)-1)* ...
            SACSCCmetrics{BOOTrep}.sums.sumCSD_AC(CFsum_index)/SACSCCmetrics{BOOTrep}.sums.sumCSD_AC(NYQsum_index);
    end
    % 3) alternative adjusted SC peak height [0,300]/[0,Fs/2]
    % just for FIG 2 in NM paper
    SACSCCmetrics{BOOTrep}.SCpeaks_legend{3}='adj: 0-300';
    ENVsum_index=find((paramsOUT.PSD_LHfreqs_Hz(:,1)==0)&(paramsOUT.PSD_LHfreqs_Hz(:,2)==300));
    SACSCCmetrics{BOOTrep}.SCpeaks_A(3)=1+(max(SACSCCs.SUMCOR_A)-1)* ...
        SACSCCmetrics{BOOTrep}.sums.sumPSD_A(ENVsum_index)/SACSCCmetrics{BOOTrep}.sums.sumPSD_A(NYQsum_index);
    if ~isnan(NumDrivenSpikes(2,1))
        SACSCCmetrics{BOOTrep}.SCpeaks_B(3)=1+(max(SACSCCs.SUMCOR_B)-1)* ...
            SACSCCmetrics{BOOTrep}.sums.sumPSD_B(ENVsum_index)/SACSCCmetrics{BOOTrep}.sums.sumPSD_B(NYQsum_index);
        SACSCCmetrics{BOOTrep}.SCpeaks_AC(3)=1+(max(SACSCCs.SUMCOR_AC)-1)* ...
            SACSCCmetrics{BOOTrep}.sums.sumCSD_AC(ENVsum_index)/SACSCCmetrics{BOOTrep}.sums.sumCSD_AC(NYQsum_index);
    end
    % 4) raw peaks of IFFTraw
    SACSCCmetrics{BOOTrep}.SCpeaks_legend{4}='IFFTraw';
    SACSCCmetrics{BOOTrep}.SCpeaks_A(4)=max(SACSCCs.SUMCORadj_A);
    if ~isnan(NumDrivenSpikes(2,1))
        SACSCCmetrics{BOOTrep}.SCpeaks_B(4)=max(SACSCCs.SUMCORadj_B);
        SACSCCmetrics{BOOTrep}.SCpeaks_C(4)=max(SACSCCs.SUMCORadj_C);
        SACSCCmetrics{BOOTrep}.SCpeaks_AC(4)=max(SACSCCs.SUMCORadj_AC);
    end
    % 5) IFFTraw (0 delay)
    SACSCCmetrics{BOOTrep}.SCpeaks_legend{5}='IFFTraw_0';
    INDEX_0=find(SACSCCs.delays_usec==0);
    SACSCCmetrics{BOOTrep}.SCpeaks_A(5)=SACSCCs.SUMCORadj_A(INDEX_0);
    if ~isnan(NumDrivenSpikes(2,1))
        SACSCCmetrics{BOOTrep}.SCpeaks_B(5)=SACSCCs.SUMCORadj_B(INDEX_0);
        SACSCCmetrics{BOOTrep}.SCpeaks_C(5)=SACSCCs.SUMCORadj_C(INDEX_0);
        SACSCCmetrics{BOOTrep}.SCpeaks_AC(5)=SACSCCs.SUMCORadj_AC(INDEX_0);
    end
    %%%Add SUMCOR peaks from 0-64, 0-5, 5-64, 64-300 and 0-300 here%%%%
    
    SACSCCmetrics{BOOTrep}.SCpeaks_legend{end+1}='IFFTraw_5@0';
    INDEX_0=find(SACSCCs.delays_usec==0);
    SACSCCmetrics{BOOTrep}.SCpeaks_A(end+1)=SACSCCs.SUMCORadj_A_5(INDEX_0);
    if ~isnan(NumDrivenSpikes(2,1))
        SACSCCmetrics{BOOTrep}.SCpeaks_B(end+1)=SACSCCs.SUMCORadj_B_5(INDEX_0);
        SACSCCmetrics{BOOTrep}.SCpeaks_C(end+1)=SACSCCs.SUMCORadj_C_5(INDEX_0);
        SACSCCmetrics{BOOTrep}.SCpeaks_AC(end+1)=SACSCCs.SUMCORadj_AC_5(INDEX_0);
    end
    
    SACSCCmetrics{BOOTrep}.SCpeaks_legend{end+1}='IFFTraw_64@0';
    INDEX_0=find(SACSCCs.delays_usec==0);
    SACSCCmetrics{BOOTrep}.SCpeaks_A(end+1)=SACSCCs.SUMCORadj_A_64(INDEX_0);
    if ~isnan(NumDrivenSpikes(2,1))
        SACSCCmetrics{BOOTrep}.SCpeaks_B(end+1)=SACSCCs.SUMCORadj_B_64(INDEX_0);
        SACSCCmetrics{BOOTrep}.SCpeaks_C(end+1)=SACSCCs.SUMCORadj_C_64(INDEX_0);
        SACSCCmetrics{BOOTrep}.SCpeaks_AC(end+1)=SACSCCs.SUMCORadj_AC_64(INDEX_0);
    end
    
    SACSCCmetrics{BOOTrep}.SCpeaks_legend{end+1}='IFFTraw_300@0';
    INDEX_0=find(SACSCCs.delays_usec==0);
    SACSCCmetrics{BOOTrep}.SCpeaks_A(end+1)=SACSCCs.SUMCORadj_A_300(INDEX_0);
    if ~isnan(NumDrivenSpikes(2,1))
        SACSCCmetrics{BOOTrep}.SCpeaks_B(end+1)=SACSCCs.SUMCORadj_B_300(INDEX_0);
        SACSCCmetrics{BOOTrep}.SCpeaks_C(end+1)=SACSCCs.SUMCORadj_C_300(INDEX_0);
        SACSCCmetrics{BOOTrep}.SCpeaks_AC(end+1)=SACSCCs.SUMCORadj_AC_300(INDEX_0);
    end
    
    SACSCCmetrics{BOOTrep}.SCpeaks_legend{end+1}='IFFTraw_5_64@0';
    INDEX_0=find(SACSCCs.delays_usec==0);
    SACSCCmetrics{BOOTrep}.SCpeaks_A(end+1)=SACSCCs.SUMCORadj_A_5_64(INDEX_0);
    if ~isnan(NumDrivenSpikes(2,1))
        SACSCCmetrics{BOOTrep}.SCpeaks_B(end+1)=SACSCCs.SUMCORadj_B_5_64(INDEX_0);
        SACSCCmetrics{BOOTrep}.SCpeaks_C(end+1)=SACSCCs.SUMCORadj_C_5_64(INDEX_0);
        SACSCCmetrics{BOOTrep}.SCpeaks_AC(end+1)=SACSCCs.SUMCORadj_AC_5_64(INDEX_0);
    end
    
    SACSCCmetrics{BOOTrep}.SCpeaks_legend{end+1}='IFFTraw_64_300@0';
    INDEX_0=find(SACSCCs.delays_usec==0);
    SACSCCmetrics{BOOTrep}.SCpeaks_A(end+1)=SACSCCs.SUMCORadj_A_64_300(INDEX_0);
    if ~isnan(NumDrivenSpikes(2,1))
        SACSCCmetrics{BOOTrep}.SCpeaks_B(end+1)=SACSCCs.SUMCORadj_B_64_300(INDEX_0);
        SACSCCmetrics{BOOTrep}.SCpeaks_C(end+1)=SACSCCs.SUMCORadj_C_64_300(INDEX_0);
        SACSCCmetrics{BOOTrep}.SCpeaks_AC(end+1)=SACSCCs.SUMCORadj_AC_64_300(INDEX_0);
    end
    % 6) IFFTraw (CD delay)
    SACSCCmetrics{BOOTrep}.SCpeaks_legend{end+1}='IFFTraw_CD';
    INDEX_CD=find(SACSCCs.delays_usec==SACSCCmetrics{BOOTrep}.CDenv_usec);
    SACSCCmetrics{BOOTrep}.SCpeaks_A(end+1)=SACSCCs.SUMCORadj_A(INDEX_CD);
    if ~isnan(NumDrivenSpikes(2,1))
        SACSCCmetrics{BOOTrep}.SCpeaks_B(end+1)=SACSCCs.SUMCORadj_B(INDEX_CD);
        SACSCCmetrics{BOOTrep}.SCpeaks_C(end+1)=SACSCCs.SUMCORadj_C(INDEX_CD);
        SACSCCmetrics{BOOTrep}.SCpeaks_AC(end+1)=SACSCCs.SUMCORadj_AC(INDEX_CD);
    end
    %% User-passed list of delays to compute
    if isfield(paramsOUT,'UserDelays_usec')
        for i=1:length(paramsOUT.UserDelays_usec)
            % 7+) IFFTraw (user delays)
            SACSCCmetrics{BOOTrep}.SCpeaks_legend{end+1}=sprintf('IFFTraw_%d',round(paramsOUT.UserDelays_usec(i)));
            [~,INDEX_user]=min(abs(SACSCCs.delays_usec-paramsOUT.UserDelays_usec(i))); % find closest delay
            SACSCCmetrics{BOOTrep}.SCpeaks_A(end+1)=SACSCCs.SUMCORadj_A(INDEX_user);
            if ~isnan(NumDrivenSpikes(2,1))
                SACSCCmetrics{BOOTrep}.SCpeaks_B(end+1)=SACSCCs.SUMCORadj_B(INDEX_user);
                SACSCCmetrics{BOOTrep}.SCpeaks_AC(end+1)=SACSCCs.SUMCORadj_AC(INDEX_user);
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%
    %% Neural Cross Correlation Coefficients
    %%%%%%%%%%%%%%%%%
    if ~isnan(NumDrivenSpikes(2,1))
        %% CCCtfs
        DC_max_index=find(strcmp(SACSCCmetrics{BOOTrep}.DCpeaks_legend,'DCpeak_max'));
        if (SACSCCmetrics{BOOTrep}.DCpeak_A(DC_max_index)>=paramsOUT.minDCpeak_CCCtfs) ...
                && (SACSCCmetrics{BOOTrep}.DCpeak_B(DC_max_index)>=paramsOUT.minDCpeak_CCCtfs)
            % 1) using DCpeak_max
            SACSCCmetrics{BOOTrep}.CCCtfs_legend{1}='DCpeak_max';
            SACSCCmetrics{BOOTrep}.CCCtfs(1) = SACSCCmetrics{BOOTrep}.DCpeak_AC(DC_max_index) ...
                / (sqrt(SACSCCmetrics{BOOTrep}.DCpeak_A(DC_max_index)*SACSCCmetrics{BOOTrep}.DCpeak_B(DC_max_index)));
            % 2) using DCpeak_0
            DC_0_index=find(strcmp(SACSCCmetrics{BOOTrep}.DCpeaks_legend,'DCpeak_0'));
            SACSCCmetrics{BOOTrep}.CCCtfs_legend{2}='DCpeak_0';
            SACSCCmetrics{BOOTrep}.CCCtfs(2) = SACSCCmetrics{BOOTrep}.DCpeak_AC(DC_0_index) ...
                / (sqrt(SACSCCmetrics{BOOTrep}.DCpeak_A(DC_0_index)*SACSCCmetrics{BOOTrep}.DCpeak_B(DC_0_index)));
            % 3) using DCpeak_CD
            DC_CD_index= strcmp(SACSCCmetrics{BOOTrep}.DCpeaks_legend,'DCpeak_CD');
            SACSCCmetrics{BOOTrep}.CCCtfs_legend{3}='DCpeak_CD';
            % NOTE: peaks for ACFs are taken at 0 delay (by definition),
            % rather than at specified CCF delay
            SACSCCmetrics{BOOTrep}.CCCtfs(3) = SACSCCmetrics{BOOTrep}.DCpeak_AC(DC_CD_index) ...
                / (sqrt(SACSCCmetrics{BOOTrep}.DCpeak_A(DC_0_index)*SACSCCmetrics{BOOTrep}.DCpeak_B(DC_0_index)));
            % 4+) User-passed list of delays to compute
            if isfield(paramsOUT,'UserDelays_usec')
                for i=1:length(paramsOUT.UserDelays_usec)
                    DC_user_index= strcmp(SACSCCmetrics{BOOTrep}.DCpeaks_legend,sprintf('DCpeak_%d',round(paramsOUT.UserDelays_usec(i))));
                    SACSCCmetrics{BOOTrep}.CCCtfs_legend{end+1}=sprintf('DCpeak_%d',round(paramsOUT.UserDelays_usec(i)));
                    % NOTE: peaks for ACFs are taken at 0 delay (by definition),
                    % rather than at specified CCF delay
                    SACSCCmetrics{BOOTrep}.CCCtfs(end+1) = SACSCCmetrics{BOOTrep}.DCpeak_AC(DC_user_index) ...
                        / (sqrt(SACSCCmetrics{BOOTrep}.DCpeak_A(DC_0_index)*SACSCCmetrics{BOOTrep}.DCpeak_B(DC_0_index)));
                end
            end
            
        else
            for i=length(SACSCCmetrics{BOOTrep}.DCpeaks_legend)
                SACSCCmetrics{BOOTrep}.CCCtfs(i) = NaN;
            end
        end
        
        %% CCCenv
        % first, go through whole LHfreq list, for subBIAS (subtract bias) and
        % withBIAS (Don't subtract bias)
        for i=1:size(paramsOUT.PSD_LHfreqs_Hz,1)
            if paramsOUT.PSD_LHfreqs_Hz(i,2)==CF_Hz
                SACSCCmetrics{BOOTrep}.CCCenvs_legend{(i-1)*2+1}=sprintf('%.f-CF, subBIAS',paramsOUT.PSD_LHfreqs_Hz(i,1));
                SACSCCmetrics{BOOTrep}.CCCenvs_legend{(i-1)*2+2}=sprintf('%.f-CF, withBIAS',paramsOUT.PSD_LHfreqs_Hz(i,1));
            elseif paramsOUT.PSD_LHfreqs_Hz(i,2)==NYQ_Hz
                SACSCCmetrics{BOOTrep}.CCCenvs_legend{(i-1)*2+1}=sprintf('%.f-NYQ, subBIAS',paramsOUT.PSD_LHfreqs_Hz(i,1));
                SACSCCmetrics{BOOTrep}.CCCenvs_legend{(i-1)*2+2}=sprintf('%.f-NYQ, withBIAS',paramsOUT.PSD_LHfreqs_Hz(i,1));
            else
                SACSCCmetrics{BOOTrep}.CCCenvs_legend{(i-1)*2+1}=sprintf('%.f-%.f, subBIAS',paramsOUT.PSD_LHfreqs_Hz(i,1),paramsOUT.PSD_LHfreqs_Hz(i,2));
                SACSCCmetrics{BOOTrep}.CCCenvs_legend{(i-1)*2+2}=sprintf('%.f-%.f, withBIAS',paramsOUT.PSD_LHfreqs_Hz(i,1),paramsOUT.PSD_LHfreqs_Hz(i,2));
            end
            % subBIAS
            SACSCCmetrics{BOOTrep}.CCCenvs((i-1)*2+1) = max([0 (SACSCCmetrics{BOOTrep}.sums.sumCSD_AC(i)-SACSCCmetrics{BOOTrep}.sums.sumCSDrand_AC(i)) ])/ ...
                sqrt(max([0 (SACSCCmetrics{BOOTrep}.sums.sumPSD_A(i)-SACSCCmetrics{BOOTrep}.sums.sumPSDrand_A(i)) ])* ...
                max([0 (SACSCCmetrics{BOOTrep}.sums.sumPSD_B(i)-SACSCCmetrics{BOOTrep}.sums.sumPSDrand_B(i)) ]));
            % withBIAS
            SACSCCmetrics{BOOTrep}.CCCenvs((i-1)*2+2) = SACSCCmetrics{BOOTrep}.sums.sumCSD_AC(i)/ ...
                sqrt(SACSCCmetrics{BOOTrep}.sums.sumPSD_A(i)*SACSCCmetrics{BOOTrep}.sums.sumPSD_B(i));
        end
        % using raw SCpeaks
        RAWSCindex=find(strcmp(SACSCCmetrics{BOOTrep}.SCpeaks_legend,'raw'));
        SACSCCmetrics{BOOTrep}.CCCenvs_legend{end+1}='rawSC';
        SACSCCmetrics{BOOTrep}.CCCenvs(end+1) = (SACSCCmetrics{BOOTrep}.SCpeaks_AC(RAWSCindex)-1)/ ...
            (sqrt((SACSCCmetrics{BOOTrep}.SCpeaks_A(RAWSCindex)-1)*(SACSCCmetrics{BOOTrep}.SCpeaks_B(RAWSCindex)-1)));
        % using adj SCpeaks
        ADJSCindex=find(strcmp(SACSCCmetrics{BOOTrep}.SCpeaks_legend,'adj: 0-CF'));
        SACSCCmetrics{BOOTrep}.CCCenvs_legend{end+1}='adjSC';
        SACSCCmetrics{BOOTrep}.CCCenvs(end+1) = (SACSCCmetrics{BOOTrep}.SCpeaks_AC(ADJSCindex)-1)/ ...
            (sqrt((SACSCCmetrics{BOOTrep}.SCpeaks_A(ADJSCindex)-1)*(SACSCCmetrics{BOOTrep}.SCpeaks_B(ADJSCindex)-1)));
        % using IFFTraw SCpeaks
        IFFTSCindex=find(strcmp(SACSCCmetrics{BOOTrep}.SCpeaks_legend,'IFFTraw'));
        SACSCCmetrics{BOOTrep}.CCCenvs_legend{end+1}='IFFTrawSC';
        SACSCCmetrics{BOOTrep}.CCCenvs(end+1) = (SACSCCmetrics{BOOTrep}.SCpeaks_AC(IFFTSCindex)-1)/ ...
            (sqrt((SACSCCmetrics{BOOTrep}.SCpeaks_A(IFFTSCindex)-1)*(SACSCCmetrics{BOOTrep}.SCpeaks_B(IFFTSCindex)-1)));
        % using IFFT(0 delay) SCpeaks
        IFFTSC_0_index=find(strcmp(SACSCCmetrics{BOOTrep}.SCpeaks_legend,'IFFTraw_0'));
        SACSCCmetrics{BOOTrep}.CCCenvs_legend{end+1}='IFFTrawSC_0';
        SACSCCmetrics{BOOTrep}.CCCenvs(end+1) = (SACSCCmetrics{BOOTrep}.SCpeaks_AC(IFFTSC_0_index)-1)/ ...
            (sqrt((SACSCCmetrics{BOOTrep}.SCpeaks_A(IFFTSC_0_index)-1)*(SACSCCmetrics{BOOTrep}.SCpeaks_B(IFFTSC_0_index)-1)));
        % using IFFT(CD delay) SCpeaks
        IFFTSC_CD_index= strcmp(SACSCCmetrics{BOOTrep}.SCpeaks_legend,'IFFTraw_CD');
        SACSCCmetrics{BOOTrep}.CCCenvs_legend{end+1}='IFFTrawSC_CD';
        % NOTE: peaks for ACFs are taken at 0 delay (by definition), rather
        % than at specified CCF delay
        SACSCCmetrics{BOOTrep}.CCCenvs(end+1) = (SACSCCmetrics{BOOTrep}.SCpeaks_AC(IFFTSC_CD_index)-1)/ ...
            (sqrt((SACSCCmetrics{BOOTrep}.SCpeaks_A(IFFTSC_0_index)-1)*(SACSCCmetrics{BOOTrep}.SCpeaks_B(IFFTSC_0_index)-1)));
        % User-passed list of delays to compute
        if isfield(paramsOUT,'UserDelays_usec')
            for i=1:length(paramsOUT.UserDelays_usec)
                IFFTSC_user_index= strcmp(SACSCCmetrics{BOOTrep}.SCpeaks_legend,sprintf('IFFTraw_%d',round(paramsOUT.UserDelays_usec(i))));
                SACSCCmetrics{BOOTrep}.CCCenvs_legend{end+1}=sprintf('IFFTrawSC_%d',round(paramsOUT.UserDelays_usec(i)));
                % NOTE: peaks for ACFs are taken at 0 delay (by definition),
                % rather than at specified CCF delay
                SACSCCmetrics{BOOTrep}.CCCenvs(end+1) = (SACSCCmetrics{BOOTrep}.SCpeaks_AC(IFFTSC_user_index)-1)/ ...
                    (sqrt((SACSCCmetrics{BOOTrep}.SCpeaks_A(IFFTSC_0_index)-1)*(SACSCCmetrics{BOOTrep}.SCpeaks_B(IFFTSC_0_index)-1)));
            end
        end
        
    end
    
    %% Book-keeping
    SACSCCfunctions{BOOTrep}=SACSCCs;
    if ~isnan(NumDrivenSpikes(2,1))
        SACSCCfunctions{BOOTrep}.rand.PSDenv_A=SACSCCs_rand.PSDenv_A;
        SACSCCfunctions{BOOTrep}.rand.PSDenv_B=SACSCCs_rand.PSDenv_B;
        SACSCCfunctions{BOOTrep}.rand.PSDenv_C=SACSCCs_rand.PSDenv_C;
        SACSCCfunctions{BOOTrep}.rand.PSDenv_AC=SACSCCs_rand.PSDenv_AC;
        SACSCCfunctions{BOOTrep}.rand.PSDtfs_A=SACSCCs_rand.PSDtfs_A;
        SACSCCfunctions{BOOTrep}.rand.PSDtfs_B=SACSCCs_rand.PSDtfs_B;
        SACSCCfunctions{BOOTrep}.rand.PSDtfs_C=SACSCCs_rand.PSDtfs_C;
        SACSCCfunctions{BOOTrep}.rand.PSDtfs_AC=SACSCCs_rand.PSDtfs_AC;
    end
    
    SACSCCmetrics{BOOTrep}.NumDrivenSpikes=NumDrivenSpikes;
    SACSCCmetrics{BOOTrep}.AvgRate_sps=AvgRate_sps;
    
end

if PLOT_15panel
    Library.plot_SACSCCanal_SNRenv(SACSCCfunctions,SACSCCmetrics,paramsOUT)
    Library.saveFigureAs([resultDir 'eps' filesep 'sac' resultPostfix '.eps']);
    Library.saveFigureAs([resultDir 'png' filesep 'sac' resultPostfix '.png']);
    Library.plot_CCCanal_NRSA(SACSCCfunctions,SACSCCmetrics,paramsOUT);
    Library.saveFigureAs([resultDir 'eps' filesep 'sacScc' resultPostfix '.eps']);
    Library.saveFigureAs([resultDir 'png' filesep 'sacScc' resultPostfix '.png']);
end

return;


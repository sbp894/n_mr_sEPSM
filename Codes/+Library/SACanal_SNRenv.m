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
	paramsOUT.MAXspikes=3700;  % If used, windowSTs will cut out extra REPs beyond MAXspikes - used to avoid memory limits
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
for i=1:3
	for j=1:2
		Nreps(i,j)=length(SpikeTrains{i,j});
		Ninclude(i,j)=ceil(paramsOUT.BOOTSTRAP_percentdata*Nreps(i,j)); % some rounding issues with ceil/floor - have to do it this way
		Nexclude(i,j)=max([1 Nreps(i,j)-Ninclude(i,j)]);  % at least 1 rep excluded to make sure runs are different
		Navgs(i,j)=floor(Nreps(i,j)/Nexclude(i,j));  % Number of RUNS for Averaging CCCs, or PSDs
	end
end
paramsOUT.BOOTSTRAP_Navgs=min(min(Navgs));  % Number of RUNS for Averaging CCCs, or PSDs
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
for BOOTrep=1:paramsOUT.BOOTSTRAP_Navgs+1
	% Run BOOTSTRAP_Navgs different random samples
	if BOOTrep~=paramsOUT.BOOTSTRAP_Navgs+1
		%% Window spike times (ignore 1st 50 ms, and cut to shorter of 2 stim)
		[ST_A_plus,NumDrivenSpikes(1,1)]=windowSTs(SpikeTrains{1,1}(BOOTinds{1,1}{BOOTrep}), ...
			paramsOUT.SCC_onsetIGNORE_sec,stimdur_sec,paramsOUT.MAXspikes);
		[ST_A_minus,NumDrivenSpikes(1,2)]=windowSTs(SpikeTrains{1,2}(BOOTinds{1,2}{BOOTrep}), ...
			paramsOUT.SCC_onsetIGNORE_sec,stimdur_sec,paramsOUT.MAXspikes);
		[ST_B_plus,NumDrivenSpikes(2,1)]=windowSTs(SpikeTrains{2,1}(BOOTinds{2,1}{BOOTrep}), ...
			paramsOUT.SCC_onsetIGNORE_sec,stimdur_sec,paramsOUT.MAXspikes);
		[ST_B_minus,NumDrivenSpikes(2,2)]=windowSTs(SpikeTrains{2,2}(BOOTinds{2,2}{BOOTrep}), ...
			paramsOUT.SCC_onsetIGNORE_sec,stimdur_sec,paramsOUT.MAXspikes);
		[ST_C_plus,NumDrivenSpikes(3,1)]=windowSTs(SpikeTrains{3,1}(BOOTinds{3,1}{BOOTrep}), ...
			paramsOUT.SCC_onsetIGNORE_sec,stimdur_sec,paramsOUT.MAXspikes);
		[ST_C_minus,NumDrivenSpikes(3,2)]=windowSTs(SpikeTrains{3,2}(BOOTinds{3,2}{BOOTrep}), ...
			paramsOUT.SCC_onsetIGNORE_sec,stimdur_sec,paramsOUT.MAXspikes);
		paramsOUT.SCCdur_sec=stimdur_sec-paramsOUT.SCC_onsetIGNORE_sec;

		NumReps(1,1)=length(ST_A_plus); NumReps(1,2)=length(ST_A_minus); 
        NumReps(2,1)=length(ST_B_plus); NumReps(2,2)=length(ST_B_minus);
        NumReps(3,1)=length(ST_C_plus); NumReps(3,2)=length(ST_C_minus);
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if verbose
        disp(sprintf('NUMBER OF SPIKES (A_plus, A_minus, B_plus, B_minus, C_plus, C_minus) = %d, %d, %d, %d, %d, %d',NumDrivenSpikes(1,1),NumDrivenSpikes(1,2), ...
			NumDrivenSpikes(2,1),NumDrivenSpikes(2,2),NumDrivenSpikes(3,1),NumDrivenSpikes(3,2)))
		disp(sprintf('NUMBER OF REPS (A_plus, A_minus, B_plus, B_minus, C_plus, C_minus) = %d, %d, %d, %d, %d, %d',NumReps(1,1),NumReps(1,2),NumReps(2,1),NumReps(2,2),NumReps(3,1),NumReps(3,2)))

		disp(sprintf('   ... Computing SACs/SCCs for AN spikes: REP %d of %d ...',BOOTrep,paramsOUT.BOOTSTRAP_Navgs))
        end
        [SACSCCs,AvgRate_sps]=SACSCCanal(ST_A_plus,ST_A_minus,ST_B_plus,ST_B_minus,ST_C_plus,ST_C_minus,paramsOUT);
        %%% OK to HERE - going NOW to SACSCCanal

        % 		if ~isnan(NumDrivenSpikes(2,1))
        %% RANDOMIZE spike times within same window used above (ignore 1st 50 ms, and cut to shorter of 2 stim)
        STrand_A_plus=randomizeSTs(ST_A_plus,paramsOUT.SCC_onsetIGNORE_sec,stimdur_sec);
        STrand_A_minus=randomizeSTs(ST_A_minus,paramsOUT.SCC_onsetIGNORE_sec,stimdur_sec);
        STrand_B_plus=randomizeSTs(ST_B_plus,paramsOUT.SCC_onsetIGNORE_sec,stimdur_sec);
        STrand_B_minus=randomizeSTs(ST_B_minus,paramsOUT.SCC_onsetIGNORE_sec,stimdur_sec);
        STrand_C_plus=randomizeSTs(ST_C_plus,paramsOUT.SCC_onsetIGNORE_sec,stimdur_sec);
        STrand_C_minus=randomizeSTs(ST_C_minus,paramsOUT.SCC_onsetIGNORE_sec,stimdur_sec);
        if verbose
        disp(sprintf('   ... Computing SACs/SCCs for RANDOM spikes (for Noise Floor) ...'))
        end
        [SACSCCs_rand,AvgRate_rand_sps]=SACSCCanal(STrand_A_plus,STrand_A_minus,STrand_B_plus,STrand_B_minus,STrand_C_plus,STrand_C_minus,paramsOUT);
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
            SACSCCs.SCC_AB_avg=SACSCCs.SCC_AB_avg+SACSCCfunctions{BOOTrep2}.SCC_AB_avg;
			SACSCCs.XpAC_A_avg=SACSCCs.XpAC_A_avg+SACSCCfunctions{BOOTrep2}.XpAC_A_avg;
			SACSCCs.XpAC_B_avg=SACSCCs.XpAC_B_avg+SACSCCfunctions{BOOTrep2}.XpAC_B_avg;
			SACSCCs.XpAC_C_avg=SACSCCs.XpAC_C_avg+SACSCCfunctions{BOOTrep2}.XpAC_C_avg;
            SACSCCs.XpCC_AB_avg=SACSCCs.XpCC_AB_avg+SACSCCfunctions{BOOTrep2}.XpCC_AB_avg;
			SACSCCs.SUMCOR_A=SACSCCs.SUMCOR_A+SACSCCfunctions{BOOTrep2}.SUMCOR_A;
			SACSCCs.SUMCOR_B=SACSCCs.SUMCOR_B+SACSCCfunctions{BOOTrep2}.SUMCOR_B;
			SACSCCs.SUMCOR_C=SACSCCs.SUMCOR_C+SACSCCfunctions{BOOTrep2}.SUMCOR_C;
            SACSCCs.SUMCOR_AB=SACSCCs.SUMCOR_AB+SACSCCfunctions{BOOTrep2}.SUMCOR_AB;
			SACSCCs.SUMCORadj_A=SACSCCs.SUMCORadj_A+SACSCCfunctions{BOOTrep2}.SUMCORadj_A;
			SACSCCs.SUMCORadj_B=SACSCCs.SUMCORadj_B+SACSCCfunctions{BOOTrep2}.SUMCORadj_B;
			SACSCCs.SUMCORadj_C=SACSCCs.SUMCORadj_C+SACSCCfunctions{BOOTrep2}.SUMCORadj_C;
            SACSCCs.SUMCORadj_AB=SACSCCs.SUMCORadj_AB+SACSCCfunctions{BOOTrep2}.SUMCORadj_AB;
			SACSCCs.DIFCOR_A=SACSCCs.DIFCOR_A+SACSCCfunctions{BOOTrep2}.DIFCOR_A;
			SACSCCs.DIFCOR_B=SACSCCs.DIFCOR_B+SACSCCfunctions{BOOTrep2}.DIFCOR_B;
			SACSCCs.DIFCOR_C=SACSCCs.DIFCOR_C+SACSCCfunctions{BOOTrep2}.DIFCOR_C;
            SACSCCs.DIFCOR_AB=SACSCCs.DIFCOR_AB+SACSCCfunctions{BOOTrep2}.DIFCOR_AB;
			SACSCCs.PSDenv_A=SACSCCs.PSDenv_A+SACSCCfunctions{BOOTrep2}.PSDenv_A;
			SACSCCs.PSDenv_B=SACSCCs.PSDenv_B+SACSCCfunctions{BOOTrep2}.PSDenv_B;
			SACSCCs.PSDenv_C=SACSCCs.PSDenv_C+SACSCCfunctions{BOOTrep2}.PSDenv_C;
            SACSCCs.PSDenv_AB=SACSCCs.PSDenv_AB+SACSCCfunctions{BOOTrep2}.PSDenv_AB;
			% RAND
			SACSCCs_rand.PSDenv_A=SACSCCs_rand.PSDenv_A+SACSCCfunctions{BOOTrep2}.rand.PSDenv_A;
			SACSCCs_rand.PSDenv_B=SACSCCs_rand.PSDenv_B+SACSCCfunctions{BOOTrep2}.rand.PSDenv_B;
			SACSCCs_rand.PSDenv_C=SACSCCs_rand.PSDenv_C+SACSCCfunctions{BOOTrep2}.rand.PSDenv_C;
            SACSCCs_rand.PSDenv_AB=SACSCCs_rand.PSDenv_AB+SACSCCfunctions{BOOTrep2}.rand.PSDenv_AB;
		end
		SACSCCs.SAC_A_avg=SACSCCs.SAC_A_avg/paramsOUT.BOOTSTRAP_Navgs;
		SACSCCs.SAC_B_avg=SACSCCs.SAC_B_avg/paramsOUT.BOOTSTRAP_Navgs;
		SACSCCs.SAC_C_avg=SACSCCs.SAC_C_avg/paramsOUT.BOOTSTRAP_Navgs;
        SACSCCs.SCC_AB_avg=SACSCCs.SCC_AB_avg/paramsOUT.BOOTSTRAP_Navgs;
		SACSCCs.XpAC_A_avg=SACSCCs.XpAC_A_avg/paramsOUT.BOOTSTRAP_Navgs;
		SACSCCs.XpAC_B_avg=SACSCCs.XpAC_B_avg/paramsOUT.BOOTSTRAP_Navgs;
		SACSCCs.XpAC_C_avg=SACSCCs.XpAC_C_avg/paramsOUT.BOOTSTRAP_Navgs;
        SACSCCs.XpCC_AB_avg=SACSCCs.XpCC_AB_avg/paramsOUT.BOOTSTRAP_Navgs;
		SACSCCs.SUMCOR_A=SACSCCs.SUMCOR_A/paramsOUT.BOOTSTRAP_Navgs;
		SACSCCs.SUMCOR_B=SACSCCs.SUMCOR_B/paramsOUT.BOOTSTRAP_Navgs;
		SACSCCs.SUMCOR_C=SACSCCs.SUMCOR_C/paramsOUT.BOOTSTRAP_Navgs;
        SACSCCs.SUMCOR_AB=SACSCCs.SUMCOR_AB/paramsOUT.BOOTSTRAP_Navgs;
		SACSCCs.SUMCORadj_A=SACSCCs.SUMCORadj_A/paramsOUT.BOOTSTRAP_Navgs;
		SACSCCs.SUMCORadj_B=SACSCCs.SUMCORadj_B/paramsOUT.BOOTSTRAP_Navgs;
		SACSCCs.SUMCORadj_C=SACSCCs.SUMCORadj_C/paramsOUT.BOOTSTRAP_Navgs;
        SACSCCs.SUMCORadj_AB=SACSCCs.SUMCORadj_AB/paramsOUT.BOOTSTRAP_Navgs;
		SACSCCs.DIFCOR_A=SACSCCs.DIFCOR_A/paramsOUT.BOOTSTRAP_Navgs;
		SACSCCs.DIFCOR_B=SACSCCs.DIFCOR_B/paramsOUT.BOOTSTRAP_Navgs;
		SACSCCs.DIFCOR_C=SACSCCs.DIFCOR_C/paramsOUT.BOOTSTRAP_Navgs;
        SACSCCs.DIFCOR_AB=SACSCCs.DIFCOR_AB/paramsOUT.BOOTSTRAP_Navgs;
		%%%%%%%%%%%%%%%%%%%%%%%
		% PSDs - Best to AVG in PSD domain, not in time domain to reduce variance (e.g., Welch technique)
		SACSCCs.PSDenv_A=SACSCCs.PSDenv_A/paramsOUT.BOOTSTRAP_Navgs;
		SACSCCs.PSDenv_B=SACSCCs.PSDenv_B/paramsOUT.BOOTSTRAP_Navgs;
		SACSCCs.PSDenv_C=SACSCCs.PSDenv_C/paramsOUT.BOOTSTRAP_Navgs;
        SACSCCs.PSDenv_AB=SACSCCs.PSDenv_AB/paramsOUT.BOOTSTRAP_Navgs;
		SACSCCs_rand.PSDenv_A=SACSCCs_rand.PSDenv_A/paramsOUT.BOOTSTRAP_Navgs;
		SACSCCs_rand.PSDenv_B=SACSCCs_rand.PSDenv_B/paramsOUT.BOOTSTRAP_Navgs;
		SACSCCs_rand.PSDenv_C=SACSCCs_rand.PSDenv_C/paramsOUT.BOOTSTRAP_Navgs;
        SACSCCs_rand.PSDenv_AB=SACSCCs_rand.PSDenv_AB/paramsOUT.BOOTSTRAP_Navgs;
        % 		%%%%%%%%%%%%%%%%%%%%%%%
        % 		%%% Recompute PSDs/CSD based on AVG SCs
        % 		%%%%%%%%%%%%%%%%%%%%%%%
        % 		SACSCCs.PSDenv_A=abs(fft((SACSCCs.SUMCOR_A-1),paramsOUT.Nfft_psd));
        % 		SACSCCs.PSDenv_B=abs(fft((SACSCCs.SUMCOR_B-1),paramsOUT.Nfft_psd));
        % 		SACSCCs.PSDenv_C=abs(fft((SACSCCs.SUMCOR_AB-1),paramsOUT.Nfft_psd));
		
		NumDrivenSpikes=zeros(3,2);
		AvgRate_sps=zeros(3,2);
		for BOOTrep3=1:paramsOUT.BOOTSTRAP_Navgs
			NumDrivenSpikes=NumDrivenSpikes+SACSCCmetrics{BOOTrep3}.NumDrivenSpikes/paramsOUT.BOOTSTRAP_Navgs;
			AvgRate_sps=AvgRate_sps+SACSCCmetrics{BOOTrep3}.AvgRate_sps/paramsOUT.BOOTSTRAP_Navgs;
		end
	end

	%%%%%%%%%%%%%%%%%%%%%
	%% PSD/CSD summations - Compute summed energy in ENVELOPE spectral
	%% densities over various frequency ranges
	%%%%%%%%%%%%%%%%%%%%%
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
	for i=1:size(paramsOUT.PSD_LHfreqs_Hz,1)
		[y,PSD_LHfreqs_inds(i,1)]=min(abs(SACSCCs.freqVEC-paramsOUT.PSD_LHfreqs_Hz(i,1)));
		[y,PSD_LHfreqs_inds(i,2)]=min(abs(SACSCCs.freqVEC-paramsOUT.PSD_LHfreqs_Hz(i,2)));
	end

	% Compute all sums for all spectral densities
	for i=1:size(paramsOUT.PSD_LHfreqs_Hz,1)
		SACSCCmetrics{BOOTrep}.sums.sumPSD_A(i) = sum(SACSCCs.PSDenv_A(PSD_LHfreqs_inds(i,1):PSD_LHfreqs_inds(i,2)));
		if ~isnan(NumDrivenSpikes(2,1))
			SACSCCmetrics{BOOTrep}.sums.sumPSD_B(i) = sum(SACSCCs.PSDenv_B(PSD_LHfreqs_inds(i,1):PSD_LHfreqs_inds(i,2)));
			SACSCCmetrics{BOOTrep}.sums.sumCSD_AB(i) = sum(SACSCCs.PSDenv_C(PSD_LHfreqs_inds(i,1):PSD_LHfreqs_inds(i,2)));
			SACSCCmetrics{BOOTrep}.sums.sumPSDrand_A(i) = sum(SACSCCs_rand.PSDenv_A(PSD_LHfreqs_inds(i,1):PSD_LHfreqs_inds(i,2)));
			SACSCCmetrics{BOOTrep}.sums.sumPSDrand_B(i) = sum(SACSCCs_rand.PSDenv_B(PSD_LHfreqs_inds(i,1):PSD_LHfreqs_inds(i,2)));
			SACSCCmetrics{BOOTrep}.sums.sumCSDrand_AB(i) = sum(SACSCCs_rand.PSDenv_C(PSD_LHfreqs_inds(i,1):PSD_LHfreqs_inds(i,2)));
		end
	end
	SACSCCmetrics{BOOTrep}.sums.PSD_LHfreqs_Hz=paramsOUT.PSD_LHfreqs_Hz;
		
	%%%%%%%%%%%%%%%%%%%%
	%% Compute metrics
	%%%%%%%%%%%%%%%%%%%%

	%% Characteristic Delays
	%% NEED TO IMPROVE - find largest, unless there is a "close-to-largest" that
	%% is closer to zero!
	% probably use some criterion for 2nd largest peak (if with 5%)???
	SACSCCmetrics{BOOTrep}.CDscc_usec=findCD_SCC(SACSCCs.SCC_AB_avg,SACSCCs.delays_usec);
	SACSCCmetrics{BOOTrep}.CDenv_usec=findCD_SCC(SACSCCs.SUMCORadj_AB,SACSCCs.delays_usec);  % use IFFTadjusted SUMCOR!
	SACSCCmetrics{BOOTrep}.CDtfs_usec=findCD_SCC(SACSCCs.DIFCOR_AB,SACSCCs.delays_usec);

	
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
		SACSCCmetrics{BOOTrep}.SCCpeak_AB(3)=SACSCCs.SCC_AB_avg(INDEX_CD);
	end
	%% User-passed list of delays to compute
	if isfield(paramsOUT,'UserDelays_usec')
		for i=1:length(paramsOUT.UserDelays_usec)
			% 4+) SACpeak (user delays)
			SACSCCmetrics{BOOTrep}.SACpeaks_legend{end+1}=sprintf('SACpeak_%d',round(paramsOUT.UserDelays_usec(i)));
			[yTEMP,INDEX_user]=min(abs(SACSCCs.delays_usec-paramsOUT.UserDelays_usec(i))); % find closest delay
			SACSCCmetrics{BOOTrep}.SACpeak_A(end+1)=SACSCCs.SAC_A_avg(INDEX_user);
			if ~isnan(NumDrivenSpikes(2,1))
				SACSCCmetrics{BOOTrep}.SACpeak_B(end+1)=SACSCCs.SAC_B_avg(INDEX_user);
				SACSCCmetrics{BOOTrep}.SCCpeak_AB(end+1)=SACSCCs.SCC_AB_avg(INDEX_user);
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
        SACSCCmetrics{BOOTrep}.DCpeak_AB(2)=SACSCCs.DIFCOR_AB(INDEX_0);
	end
	% 3) DCpeak (CD delay)
	SACSCCmetrics{BOOTrep}.DCpeaks_legend{3}='DCpeak_CD';
	INDEX_CD=find(SACSCCs.delays_usec==SACSCCmetrics{BOOTrep}.CDtfs_usec);
	SACSCCmetrics{BOOTrep}.DCpeak_A(3)=SACSCCs.DIFCOR_A(INDEX_CD);
	if ~isnan(NumDrivenSpikes(2,1))
		SACSCCmetrics{BOOTrep}.DCpeak_B(3)=SACSCCs.DIFCOR_B(INDEX_CD);
		SACSCCmetrics{BOOTrep}.DCpeak_AB(3)=SACSCCs.DIFCOR_AB(INDEX_CD);
	end
	%% User-passed list of delays to compute
	if isfield(paramsOUT,'UserDelays_usec')
		for i=1:length(paramsOUT.UserDelays_usec)
			% 4+) DCpeak (user delays)
			SACSCCmetrics{BOOTrep}.DCpeaks_legend{end+1}=sprintf('DCpeak_%d',round(paramsOUT.UserDelays_usec(i)));
			[yTEMP,INDEX_user]=min(abs(SACSCCs.delays_usec-paramsOUT.UserDelays_usec(i))); % find closest delay
			SACSCCmetrics{BOOTrep}.DCpeak_A(end+1)=SACSCCs.DIFCOR_A(INDEX_user);
			if ~isnan(NumDrivenSpikes(2,1))
				SACSCCmetrics{BOOTrep}.DCpeak_B(end+1)=SACSCCs.DIFCOR_B(INDEX_user);
				SACSCCmetrics{BOOTrep}.DCpeak_AB(end+1)=SACSCCs.DIFCOR_AB(INDEX_user);
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
		SACSCCmetrics{BOOTrep}.SCpeaks_AB(2)=1+(max(SACSCCs.SUMCOR_AB)-1)* ...
			SACSCCmetrics{BOOTrep}.sums.sumCSD_AB(CFsum_index)/SACSCCmetrics{BOOTrep}.sums.sumCSD_AB(NYQsum_index);
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
		SACSCCmetrics{BOOTrep}.SCpeaks_AB(3)=1+(max(SACSCCs.SUMCOR_AB)-1)* ...
			SACSCCmetrics{BOOTrep}.sums.sumCSD_AB(ENVsum_index)/SACSCCmetrics{BOOTrep}.sums.sumCSD_AB(NYQsum_index);
	end
	% 4) raw peaks of IFFTraw
	SACSCCmetrics{BOOTrep}.SCpeaks_legend{4}='IFFTraw';
	SACSCCmetrics{BOOTrep}.SCpeaks_A(4)=max(SACSCCs.SUMCORadj_A);
	if ~isnan(NumDrivenSpikes(2,1))
		SACSCCmetrics{BOOTrep}.SCpeaks_B(4)=max(SACSCCs.SUMCORadj_B);
		SACSCCmetrics{BOOTrep}.SCpeaks_C(4)=max(SACSCCs.SUMCORadj_C);
        SACSCCmetrics{BOOTrep}.SCpeaks_AB(4)=max(SACSCCs.SUMCORadj_AB);
	end
	% 5) IFFTraw (0 delay)
	SACSCCmetrics{BOOTrep}.SCpeaks_legend{5}='IFFTraw_0';
	INDEX_0=find(SACSCCs.delays_usec==0);
	SACSCCmetrics{BOOTrep}.SCpeaks_A(5)=SACSCCs.SUMCORadj_A(INDEX_0);
	if ~isnan(NumDrivenSpikes(2,1))
		SACSCCmetrics{BOOTrep}.SCpeaks_B(5)=SACSCCs.SUMCORadj_B(INDEX_0);
		SACSCCmetrics{BOOTrep}.SCpeaks_C(5)=SACSCCs.SUMCORadj_C(INDEX_0);
        SACSCCmetrics{BOOTrep}.SCpeaks_AB(5)=SACSCCs.SUMCORadj_AB(INDEX_0);
    end
    %%%Add SUMCOR peaks from 0-64, 0-5, 5-64, 64-300 and 0-300 here%%%%
    
    SACSCCmetrics{BOOTrep}.SCpeaks_legend{end+1}='IFFTraw_5@0';
    INDEX_0=find(SACSCCs.delays_usec==0);
	SACSCCmetrics{BOOTrep}.SCpeaks_A(end+1)=SACSCCs.SUMCORadj_A_5(INDEX_0);
	if ~isnan(NumDrivenSpikes(2,1))
		SACSCCmetrics{BOOTrep}.SCpeaks_B(end+1)=SACSCCs.SUMCORadj_B_5(INDEX_0);
        SACSCCmetrics{BOOTrep}.SCpeaks_C(end+1)=SACSCCs.SUMCORadj_C_5(INDEX_0);
		SACSCCmetrics{BOOTrep}.SCpeaks_AB(end+1)=SACSCCs.SUMCORadj_AB_5(INDEX_0);
    end

    SACSCCmetrics{BOOTrep}.SCpeaks_legend{end+1}='IFFTraw_64@0';
    INDEX_0=find(SACSCCs.delays_usec==0);
	SACSCCmetrics{BOOTrep}.SCpeaks_A(end+1)=SACSCCs.SUMCORadj_A_64(INDEX_0);
	if ~isnan(NumDrivenSpikes(2,1))
		SACSCCmetrics{BOOTrep}.SCpeaks_B(end+1)=SACSCCs.SUMCORadj_B_64(INDEX_0);
        SACSCCmetrics{BOOTrep}.SCpeaks_C(end+1)=SACSCCs.SUMCORadj_C_64(INDEX_0);
		SACSCCmetrics{BOOTrep}.SCpeaks_AB(end+1)=SACSCCs.SUMCORadj_AB_64(INDEX_0);
    end
    
    SACSCCmetrics{BOOTrep}.SCpeaks_legend{end+1}='IFFTraw_300@0';
    INDEX_0=find(SACSCCs.delays_usec==0);
	SACSCCmetrics{BOOTrep}.SCpeaks_A(end+1)=SACSCCs.SUMCORadj_A_300(INDEX_0);
	if ~isnan(NumDrivenSpikes(2,1))
		SACSCCmetrics{BOOTrep}.SCpeaks_B(end+1)=SACSCCs.SUMCORadj_B_300(INDEX_0);
        SACSCCmetrics{BOOTrep}.SCpeaks_C(end+1)=SACSCCs.SUMCORadj_C_300(INDEX_0);
		SACSCCmetrics{BOOTrep}.SCpeaks_AB(end+1)=SACSCCs.SUMCORadj_AB_300(INDEX_0);
    end
    
    SACSCCmetrics{BOOTrep}.SCpeaks_legend{end+1}='IFFTraw_5_64@0';
    INDEX_0=find(SACSCCs.delays_usec==0);
	SACSCCmetrics{BOOTrep}.SCpeaks_A(end+1)=SACSCCs.SUMCORadj_A_5_64(INDEX_0);
	if ~isnan(NumDrivenSpikes(2,1))
		SACSCCmetrics{BOOTrep}.SCpeaks_B(end+1)=SACSCCs.SUMCORadj_B_5_64(INDEX_0);
        SACSCCmetrics{BOOTrep}.SCpeaks_C(end+1)=SACSCCs.SUMCORadj_C_5_64(INDEX_0);
		SACSCCmetrics{BOOTrep}.SCpeaks_AB(end+1)=SACSCCs.SUMCORadj_AB_5_64(INDEX_0);
    end

    SACSCCmetrics{BOOTrep}.SCpeaks_legend{end+1}='IFFTraw_64_300@0';
    INDEX_0=find(SACSCCs.delays_usec==0);
	SACSCCmetrics{BOOTrep}.SCpeaks_A(end+1)=SACSCCs.SUMCORadj_A_64_300(INDEX_0);
	if ~isnan(NumDrivenSpikes(2,1))
		SACSCCmetrics{BOOTrep}.SCpeaks_B(end+1)=SACSCCs.SUMCORadj_B_64_300(INDEX_0);
        SACSCCmetrics{BOOTrep}.SCpeaks_C(end+1)=SACSCCs.SUMCORadj_C_64_300(INDEX_0);
		SACSCCmetrics{BOOTrep}.SCpeaks_AB(end+1)=SACSCCs.SUMCORadj_AB_64_300(INDEX_0);
    end
	% 6) IFFTraw (CD delay)
	SACSCCmetrics{BOOTrep}.SCpeaks_legend{end+1}='IFFTraw_CD';
	INDEX_CD=find(SACSCCs.delays_usec==SACSCCmetrics{BOOTrep}.CDenv_usec);
	SACSCCmetrics{BOOTrep}.SCpeaks_A(end+1)=SACSCCs.SUMCORadj_A(INDEX_CD);
	if ~isnan(NumDrivenSpikes(2,1))
		SACSCCmetrics{BOOTrep}.SCpeaks_B(end+1)=SACSCCs.SUMCORadj_B(INDEX_CD);
        SACSCCmetrics{BOOTrep}.SCpeaks_C(end+1)=SACSCCs.SUMCORadj_C(INDEX_CD);
		SACSCCmetrics{BOOTrep}.SCpeaks_AB(end+1)=SACSCCs.SUMCORadj_AB(INDEX_CD);
	end
	%% User-passed list of delays to compute
	if isfield(paramsOUT,'UserDelays_usec')
		for i=1:length(paramsOUT.UserDelays_usec)
			% 7+) IFFTraw (user delays)
			SACSCCmetrics{BOOTrep}.SCpeaks_legend{end+1}=sprintf('IFFTraw_%d',round(paramsOUT.UserDelays_usec(i)));
			[yTEMP,INDEX_user]=min(abs(SACSCCs.delays_usec-paramsOUT.UserDelays_usec(i))); % find closest delay
			SACSCCmetrics{BOOTrep}.SCpeaks_A(end+1)=SACSCCs.SUMCORadj_A(INDEX_user);
			if ~isnan(NumDrivenSpikes(2,1))
				SACSCCmetrics{BOOTrep}.SCpeaks_B(end+1)=SACSCCs.SUMCORadj_B(INDEX_user);
				SACSCCmetrics{BOOTrep}.SCpeaks_AB(end+1)=SACSCCs.SUMCORadj_AB(INDEX_user);
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
				& (SACSCCmetrics{BOOTrep}.DCpeak_B(DC_max_index)>=paramsOUT.minDCpeak_CCCtfs)
			% 1) using DCpeak_max
			SACSCCmetrics{BOOTrep}.CCCtfs_legend{1}='DCpeak_max';
			SACSCCmetrics{BOOTrep}.CCCtfs(1) = SACSCCmetrics{BOOTrep}.DCpeak_AB(DC_max_index) ...
				/ (sqrt(SACSCCmetrics{BOOTrep}.DCpeak_A(DC_max_index)*SACSCCmetrics{BOOTrep}.DCpeak_B(DC_max_index)));
			% 2) using DCpeak_0
			DC_0_index=find(strcmp(SACSCCmetrics{BOOTrep}.DCpeaks_legend,'DCpeak_0'));
			SACSCCmetrics{BOOTrep}.CCCtfs_legend{2}='DCpeak_0';
			SACSCCmetrics{BOOTrep}.CCCtfs(2) = SACSCCmetrics{BOOTrep}.DCpeak_AB(DC_0_index) ...
				/ (sqrt(SACSCCmetrics{BOOTrep}.DCpeak_A(DC_0_index)*SACSCCmetrics{BOOTrep}.DCpeak_B(DC_0_index)));
			% 3) using DCpeak_CD
			DC_CD_index=find(strcmp(SACSCCmetrics{BOOTrep}.DCpeaks_legend,'DCpeak_CD'));
			SACSCCmetrics{BOOTrep}.CCCtfs_legend{3}='DCpeak_CD';
			% NOTE: peaks for ACFs are taken at 0 delay (by definition),
			% rather than at specified CCF delay
			SACSCCmetrics{BOOTrep}.CCCtfs(3) = SACSCCmetrics{BOOTrep}.DCpeak_AB(DC_CD_index) ...
				/ (sqrt(SACSCCmetrics{BOOTrep}.DCpeak_A(DC_0_index)*SACSCCmetrics{BOOTrep}.DCpeak_B(DC_0_index)));
			% 4+) User-passed list of delays to compute
			if isfield(paramsOUT,'UserDelays_usec')
				for i=1:length(paramsOUT.UserDelays_usec)
					DC_user_index=find(strcmp(SACSCCmetrics{BOOTrep}.DCpeaks_legend,sprintf('DCpeak_%d',round(paramsOUT.UserDelays_usec(i)))));
					SACSCCmetrics{BOOTrep}.CCCtfs_legend{end+1}=sprintf('DCpeak_%d',round(paramsOUT.UserDelays_usec(i)));
					% NOTE: peaks for ACFs are taken at 0 delay (by definition),
					% rather than at specified CCF delay
					SACSCCmetrics{BOOTrep}.CCCtfs(end+1) = SACSCCmetrics{BOOTrep}.DCpeak_AB(DC_user_index) ...
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
			SACSCCmetrics{BOOTrep}.CCCenvs((i-1)*2+1) = max([0 (SACSCCmetrics{BOOTrep}.sums.sumCSD_AB(i)-SACSCCmetrics{BOOTrep}.sums.sumCSDrand_AB(i)) ])/ ...
				sqrt(max([0 (SACSCCmetrics{BOOTrep}.sums.sumPSD_A(i)-SACSCCmetrics{BOOTrep}.sums.sumPSDrand_A(i)) ])* ...
				max([0 (SACSCCmetrics{BOOTrep}.sums.sumPSD_B(i)-SACSCCmetrics{BOOTrep}.sums.sumPSDrand_B(i)) ]));
			% withBIAS
			SACSCCmetrics{BOOTrep}.CCCenvs((i-1)*2+2) = SACSCCmetrics{BOOTrep}.sums.sumCSD_AB(i)/ ...
				sqrt(SACSCCmetrics{BOOTrep}.sums.sumPSD_A(i)*SACSCCmetrics{BOOTrep}.sums.sumPSD_B(i));
		end
		% using raw SCpeaks
		RAWSCindex=find(strcmp(SACSCCmetrics{BOOTrep}.SCpeaks_legend,'raw'));
		SACSCCmetrics{BOOTrep}.CCCenvs_legend{end+1}='rawSC';
		SACSCCmetrics{BOOTrep}.CCCenvs(end+1) = (SACSCCmetrics{BOOTrep}.SCpeaks_AB(RAWSCindex)-1)/ ...
			(sqrt((SACSCCmetrics{BOOTrep}.SCpeaks_A(RAWSCindex)-1)*(SACSCCmetrics{BOOTrep}.SCpeaks_B(RAWSCindex)-1)));
		% using adj SCpeaks
		ADJSCindex=find(strcmp(SACSCCmetrics{BOOTrep}.SCpeaks_legend,'adj: 0-CF'));
		SACSCCmetrics{BOOTrep}.CCCenvs_legend{end+1}='adjSC';
		SACSCCmetrics{BOOTrep}.CCCenvs(end+1) = (SACSCCmetrics{BOOTrep}.SCpeaks_AB(ADJSCindex)-1)/ ...
			(sqrt((SACSCCmetrics{BOOTrep}.SCpeaks_A(ADJSCindex)-1)*(SACSCCmetrics{BOOTrep}.SCpeaks_B(ADJSCindex)-1)));
		% using IFFTraw SCpeaks
		IFFTSCindex=find(strcmp(SACSCCmetrics{BOOTrep}.SCpeaks_legend,'IFFTraw'));
		SACSCCmetrics{BOOTrep}.CCCenvs_legend{end+1}='IFFTrawSC';
		SACSCCmetrics{BOOTrep}.CCCenvs(end+1) = (SACSCCmetrics{BOOTrep}.SCpeaks_AB(IFFTSCindex)-1)/ ...
			(sqrt((SACSCCmetrics{BOOTrep}.SCpeaks_A(IFFTSCindex)-1)*(SACSCCmetrics{BOOTrep}.SCpeaks_B(IFFTSCindex)-1)));
		% using IFFT(0 delay) SCpeaks
		IFFTSC_0_index=find(strcmp(SACSCCmetrics{BOOTrep}.SCpeaks_legend,'IFFTraw_0'));
		SACSCCmetrics{BOOTrep}.CCCenvs_legend{end+1}='IFFTrawSC_0';
		SACSCCmetrics{BOOTrep}.CCCenvs(end+1) = (SACSCCmetrics{BOOTrep}.SCpeaks_AB(IFFTSC_0_index)-1)/ ...
			(sqrt((SACSCCmetrics{BOOTrep}.SCpeaks_A(IFFTSC_0_index)-1)*(SACSCCmetrics{BOOTrep}.SCpeaks_B(IFFTSC_0_index)-1)));
		% using IFFT(CD delay) SCpeaks
		IFFTSC_CD_index=find(strcmp(SACSCCmetrics{BOOTrep}.SCpeaks_legend,'IFFTraw_CD'));
		SACSCCmetrics{BOOTrep}.CCCenvs_legend{end+1}='IFFTrawSC_CD';
		% NOTE: peaks for ACFs are taken at 0 delay (by definition), rather
		% than at specified CCF delay 
		SACSCCmetrics{BOOTrep}.CCCenvs(end+1) = (SACSCCmetrics{BOOTrep}.SCpeaks_AB(IFFTSC_CD_index)-1)/ ...
			(sqrt((SACSCCmetrics{BOOTrep}.SCpeaks_A(IFFTSC_0_index)-1)*(SACSCCmetrics{BOOTrep}.SCpeaks_B(IFFTSC_0_index)-1)));
		% User-passed list of delays to compute
		if isfield(paramsOUT,'UserDelays_usec')
			for i=1:length(paramsOUT.UserDelays_usec)
				IFFTSC_user_index=find(strcmp(SACSCCmetrics{BOOTrep}.SCpeaks_legend,sprintf('IFFTraw_%d',round(paramsOUT.UserDelays_usec(i)))));
				SACSCCmetrics{BOOTrep}.CCCenvs_legend{end+1}=sprintf('IFFTrawSC_%d',round(paramsOUT.UserDelays_usec(i)));
				% NOTE: peaks for ACFs are taken at 0 delay (by definition),
				% rather than at specified CCF delay
				SACSCCmetrics{BOOTrep}.CCCenvs(end+1) = (SACSCCmetrics{BOOTrep}.SCpeaks_AB(IFFTSC_user_index)-1)/ ...
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
        SACSCCfunctions{BOOTrep}.rand.PSDenv_AB=SACSCCs_rand.PSDenv_AB;
	end

	SACSCCmetrics{BOOTrep}.NumDrivenSpikes=NumDrivenSpikes;
	SACSCCmetrics{BOOTrep}.AvgRate_sps=AvgRate_sps;

% 	if PLOT_15panel
% 		plot_CCCanal_0(SACSCCfunctions{BOOTrep},SACSCCmetrics{BOOTrep},paramsOUT)
% 	end

end


% for BOOTrep10=1:length(SACSCCfunctions)
% 	SACSCCfunctions{BOOTrep10}=rmfield(SACSCCfunctions{BOOTrep10},'PSDenv_A');
% 	SACSCCfunctions{BOOTrep10}=rmfield(SACSCCfunctions{BOOTrep10},'PSDenv_B');
% 	SACSCCfunctions{BOOTrep10}=rmfield(SACSCCfunctions{BOOTrep10},'PSDenv_C');
% 	SACSCCfunctions{BOOTrep10}=rmfield(SACSCCfunctions{BOOTrep10},'rand');
% 	SACSCCfunctions{BOOTrep10}=rmfield(SACSCCfunctions{BOOTrep10},'freqVEC');
% end



% %% TALLY AVG
% CCCtfsVECs=cell(size(SACSCCmetrics{1}.CCCtfs));
% CCCtfsAVGs=zeros(size(SACSCCmetrics{1}.CCCtfs));
% CCCtfsSTDs=zeros(size(SACSCCmetrics{1}.CCCtfs));
% for i=1:length(CCCtfsVECs)
% 	CCCtfsVECs{i}=zeros(1,paramsOUT.BOOTSTRAP_Navgs);
% 	for BOOTrep4=1:paramsOUT.BOOTSTRAP_Navgs
% 		CCCtfsVECs{i}(BOOTrep4)=SACSCCmetrics{BOOTrep4}.CCCtfs(i);
% 	end
% 	CCCtfsAVGs(i)=mean(CCCtfsVECs{i});
% 	CCCtfsSTDs(i)=std(CCCtfsVECs{i});
% end

% CCCenvVECs=cell(size(SACSCCmetrics{1}.CCCenvs));
% CCCenvAVGs=zeros(size(SACSCCmetrics{1}.CCCenvs));
% CCCenvSTDs=zeros(size(SACSCCmetrics{1}.CCCenvs));
% for i=1:length(CCCenvVECs)
% 	CCCenvVECs{i}=zeros(1,paramsOUT.BOOTSTRAP_Navgs);
% 	for BOOTrep4=1:paramsOUT.BOOTSTRAP_Navgs
% 		CCCenvVECs{i}(BOOTrep4)=SACSCCmetrics{BOOTrep4}.CCCenvs(i);
% 	end
% 	CCCenvAVGs(i)=mean(CCCenvVECs{i});
% 	CCCenvSTDs(i)=std(CCCenvVECs{i});
% end

% SACSCCmetrics{paramsOUT.BOOTSTRAP_Navgs+1}.CCCenvVECs=CCCenvVECs;
% SACSCCmetrics{paramsOUT.BOOTSTRAP_Navgs+1}.CCCenvAVGs=CCCenvAVGs;
% SACSCCmetrics{paramsOUT.BOOTSTRAP_Navgs+1}.CCCenvSTDs=CCCenvSTDs;
% SACSCCmetrics{paramsOUT.BOOTSTRAP_Navgs+1}.CCCtfsVECs=CCCtfsVECs;
% SACSCCmetrics{paramsOUT.BOOTSTRAP_Navgs+1}.CCCtfsAVGs=CCCtfsAVGs;
% SACSCCmetrics{paramsOUT.BOOTSTRAP_Navgs+1}.CCCtfsSTDs=CCCtfsSTDs;

if PLOT_15panel
	plot_SACSCCanal_SNRenv(SACSCCfunctions,SACSCCmetrics,paramsOUT)
    Library.saveFigureAs([resultDir 'doc/sac' resultPostfix{1} '_50ms.eps'])
    Library.saveFigureAs([resultDir 'doc/png/sac' resultPostfix{1} '_50ms.png'])
    plot_CCCanal_NRSA(SACSCCfunctions,SACSCCmetrics,paramsOUT)
    Library.saveFigureAs([resultDir 'doc/sacScc' resultPostfix{1} '_50ms.eps'])
    Library.saveFigureAs([resultDir 'doc/png/sacScc' resultPostfix{1} '_50ms.png'])
end



return;

function [SACSCCs,AvgRate_sps]=SACSCCanal(ST_A_plus,ST_A_minus,ST_B_plus,ST_B_minus,ST_C_plus,ST_C_minus,paramsOUT);

MAXdelay_ind=round(paramsOUT.MAXdelay_sec/paramsOUT.DELAYbinwidth_sec);  % old XLIM=250

%% COLUMN 1: CONDITION=A
%% Compute SAC (A+ and A-) %%%%%%%%%%%%%%%%%%%%
[SAC_A_plus,delays_usec,AvgRate_sps(1,1),xxNspikes] = Library.ShufAutoCorr(ST_A_plus,paramsOUT.DELAYbinwidth_sec,paramsOUT.SCCdur_sec);
[SAC_A_minus,delays_usec,AvgRate_sps(1,2),xxNspikes] = Library.ShufAutoCorr(ST_A_minus,paramsOUT.DELAYbinwidth_sec,paramsOUT.SCCdur_sec);
% SAME for all conditions
ZEROind=find(delays_usec==0);
SACinds=(ZEROind-MAXdelay_ind:ZEROind+MAXdelay_ind);
SACdelays_usec=delays_usec(SACinds);
% SAC limited-data Window CORRECTION
TEMP=linspace(1,0,ZEROind);
WindowCORRECTION=[TEMP fliplr(TEMP(1:end-1))];
% % USED TO TAKE OUT WINDOW CORRECTION to test
% beep
% disp('NO WINDOW CORRECTION APPLIED')
% WindowCORRECTION=zeros(size(WindowCORRECTION));

SAC_A_plus=SAC_A_plus+WindowCORRECTION;
SAC_A_minus=SAC_A_minus+WindowCORRECTION;
% smooth and AVG correlograms
SAC_A_plus=trifilt(SAC_A_plus(SACinds),paramsOUT.TriFiltWidthSAC);
SAC_A_minus=trifilt(SAC_A_minus(SACinds),paramsOUT.TriFiltWidthSAC);
SAC_A_avg=(SAC_A_plus+SAC_A_minus)/2;

%% Compute XpAC (A+/A- and A-/A+) %%%%%%%%%%%%%%%%%%%%
[XpAC_A_plus,delays_usec,xxAVGrate,xxNspikes] = Library.ShufCrossCorr({ST_A_plus ST_A_minus},paramsOUT.DELAYbinwidth_sec,paramsOUT.SCCdur_sec);
% SAC limited-data Window CORRECTION
XpAC_A_plus=XpAC_A_plus+WindowCORRECTION;
% smooth and AVG correlograms
XpAC_A_plus=trifilt(XpAC_A_plus(SACinds),paramsOUT.TriFiltWidthSAC);
XpAC_A_minus=fliplr(XpAC_A_plus);
XpAC_A_avg=(XpAC_A_plus+XpAC_A_minus)/2;

%% Compute DIFCOR and SUMCOR for Condition A
DIFCOR_A=trifilt(SAC_A_avg-XpAC_A_avg,paramsOUT.TriFiltWidthDC);  % Diffcorr =  Avg(SAC) - Avg(XpAC)
SUMCOR_A=trifilt((SAC_A_avg+XpAC_A_avg)/2,paramsOUT.TriFiltWidthSC); %Sumcorr = 1/2(Avg(SAC) + Avg(XpAC))

%% Remove TFS artifact (centered at 2*CF) from SUMCOR
% zero out above CF
FFTtemp=fft((SUMCOR_A-1),paramsOUT.Nfft_psd);
Fs_PSD=1/paramsOUT.DELAYbinwidth_sec;  % Sampling Rate for PSDs
freqVEC=(0:paramsOUT.Nfft_psd-1)/paramsOUT.Nfft_psd*Fs_PSD;
[y,CF_index]=min(abs(freqVEC-paramsOUT.SACSCC_CF_Hz)); % use SACSCC_CF_Hz
FFTadj=zeros(size(FFTtemp));
FFTadj(1:CF_index)=FFTtemp(1:CF_index);
FFTadj((length(FFTtemp)-CF_index+1):end)=FFTtemp((length(FFTtemp)-CF_index+1):end); %keep negative freqs
adjSC=ifft(FFTadj)+1;
SUMCORadj_A=real(adjSC(1:length(SUMCOR_A)));

%%%%%SUMCOR IFFT 0-5 Hz%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[y,CF_index_5]=min(abs(freqVEC-5)); % use CF_A
FFTadj_5=zeros(size(FFTtemp));
FFTadj_5(1:CF_index_5)=FFTtemp(1:CF_index_5);
FFTadj_5((length(FFTtemp)-CF_index_5+1):end)=FFTtemp((length(FFTtemp)-CF_index_5+1):end); %keep negative freqs
adjSC_5=ifft(FFTadj_5)+1;
SUMCORadj_A_5=real(adjSC_5(1:length(SUMCOR_A)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%SUMCOR IFFT 0-64 Hz%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[y,CF_index_64]=min(abs(freqVEC-64)); % use CF_A
FFTadj_64=zeros(size(FFTtemp));
FFTadj_64(1:CF_index_64)=FFTtemp(1:CF_index_64);
FFTadj_64((length(FFTtemp)-CF_index_64+1):end)=FFTtemp((length(FFTtemp)-CF_index_64+1):end); %keep negative freqs
adjSC_64=ifft(FFTadj_64)+1;
SUMCORadj_A_64=real(adjSC_64(1:length(SUMCOR_A)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%SUMCOR IFFT 0-300 Hz%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[y,CF_index_300]=min(abs(freqVEC-300)); % use CF_A
FFTadj_300=zeros(size(FFTtemp));
FFTadj_300(1:CF_index_300)=FFTtemp(1:CF_index_300);
FFTadj_300((length(FFTtemp)-CF_index_300+1):end)=FFTtemp((length(FFTtemp)-CF_index_300+1):end); %keep negative freqs
adjSC_300=ifft(FFTadj_300)+1;
SUMCORadj_A_300=real(adjSC_300(1:length(SUMCOR_A)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%SUMCOR IFFT 5-64 Hz%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FFTadj_5_64=zeros(size(FFTtemp));
FFTadj_5_64(CF_index_5:CF_index_64)=FFTtemp(CF_index_5:CF_index_64);
FFTadj_5_64((length(FFTtemp)-CF_index_64+1):(length(FFTtemp)-CF_index_5+1))=FFTtemp((length(FFTtemp)-CF_index_64+1):(length(FFTtemp)-CF_index_5+1)); %keep negative freqs
adjSC_5_64=ifft(FFTadj_5_64)+1;
SUMCORadj_A_5_64=real(adjSC_5_64(1:length(SUMCOR_A)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%SUMCOR IFFT 64-300 Hz%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FFTadj_64_300=zeros(size(FFTtemp));
FFTadj_64_300(CF_index_64:CF_index_300)=FFTtemp(CF_index_64:CF_index_300);
FFTadj_64_300((length(FFTtemp)-CF_index_300+1):(length(FFTtemp)-CF_index_64+1))=FFTtemp((length(FFTtemp)-CF_index_300+1):(length(FFTtemp)-CF_index_64+1)); %keep negative freqs
adjSC_64_300=ifft(FFTadj_64_300)+1;
SUMCORadj_A_64_300=real(adjSC_64_300(1:length(SUMCOR_A)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Store ENVELOPE POWER SPECTRAL DENSITY for Condition A
PSDenv_A=abs(FFTadj);  % use freqVEC to store


%% COLUMN 2: CONDITION=B
%% Compute SAC (B+ and B-) %%%%%%%%%%%%%%%%%%%%
[SAC_B_plus,delays_usec,AvgRate_sps(2,1),xxNspikes] = Library.ShufAutoCorr(ST_B_plus,paramsOUT.DELAYbinwidth_sec,paramsOUT.SCCdur_sec);
[SAC_B_minus,delays_usec,AvgRate_sps(2,2),xxNspikes] = Library.ShufAutoCorr(ST_B_minus,paramsOUT.DELAYbinwidth_sec,paramsOUT.SCCdur_sec);
% SAC limited-data Window CORRECTION
SAC_B_plus=SAC_B_plus+WindowCORRECTION;
SAC_B_minus=SAC_B_minus+WindowCORRECTION;
% smooth and AVG correlograms
SAC_B_plus=trifilt(SAC_B_plus(SACinds),paramsOUT.TriFiltWidthSAC);
SAC_B_minus=trifilt(SAC_B_minus(SACinds),paramsOUT.TriFiltWidthSAC);
SAC_B_avg=(SAC_B_plus+SAC_B_minus)/2;

%% Compute XpAC (B+/B- and B-/B+) %%%%%%%%%%%%%%%%%%%%
[XpAC_B_plus,delays_usec,xxAVGrate,xxNspikes] = Library.ShufCrossCorr({ST_B_plus ST_B_minus},paramsOUT.DELAYbinwidth_sec,paramsOUT.SCCdur_sec);
% SAC limited-data Window CORRECTION
XpAC_B_plus=XpAC_B_plus+WindowCORRECTION;
% smooth and AVG correlograms
XpAC_B_plus=trifilt(XpAC_B_plus(SACinds),paramsOUT.TriFiltWidthSAC);
XpAC_B_minus=fliplr(XpAC_B_plus);
XpAC_B_avg=(XpAC_B_plus+XpAC_B_minus)/2;

%% Compute DIFCOR and SUMCOR for Condition B
DIFCOR_B=trifilt(SAC_B_avg-XpAC_B_avg,paramsOUT.TriFiltWidthDC);  % Diffcorr =  Avg(SAC) - Avg(XpAC)
SUMCOR_B=trifilt((SAC_B_avg+XpAC_B_avg)/2,paramsOUT.TriFiltWidthSC); %Sumcorr = 1/2(Avg(SAC) + Avg(XpAC))

%% Remove TFS artifact (centered at 2*CF) from SUMCOR
% zero out above CF
FFTtemp=fft((SUMCOR_B-1),paramsOUT.Nfft_psd);
[y,CF_index]=min(abs(freqVEC-paramsOUT.SACSCC_CF_Hz)); % use SACSCC_CF_Hz
FFTadj=zeros(size(FFTtemp));
FFTadj(1:CF_index)=FFTtemp(1:CF_index);
FFTadj((length(FFTtemp)-CF_index+1):end)=FFTtemp((length(FFTtemp)-CF_index+1):end); %keep negative freqs
adjSC=ifft(FFTadj)+1;
SUMCORadj_B=real(adjSC(1:length(SUMCOR_B)));

%%%%%SUMCOR IFFT 0-5 Hz%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[y,CF_index_5]=min(abs(freqVEC-5)); % use CF_B
FFTadj_5=zeros(size(FFTtemp));
FFTadj_5(1:CF_index_5)=FFTtemp(1:CF_index_5);
FFTadj_5((length(FFTtemp)-CF_index_5+1):end)=FFTtemp((length(FFTtemp)-CF_index_5+1):end); %keep negative freqs
adjSC_5=ifft(FFTadj_5)+1;
SUMCORadj_B_5=real(adjSC_5(1:length(SUMCOR_B)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[y,CF_index_64]=min(abs(freqVEC-64)); % use CF_A
FFTadj_64=zeros(size(FFTtemp));
FFTadj_64(1:CF_index_64)=FFTtemp(1:CF_index_64);
FFTadj_64((length(FFTtemp)-CF_index_64+1):end)=FFTtemp((length(FFTtemp)-CF_index_64+1):end); %keep negative freqs
adjSC_64=ifft(FFTadj_64)+1;
SUMCORadj_B_64=real(adjSC_64(1:length(SUMCOR_B)));

[y,CF_index_300]=min(abs(freqVEC-300)); % use CF_A
FFTadj_300=zeros(size(FFTtemp));
FFTadj_300(1:CF_index_300)=FFTtemp(1:CF_index_300);
FFTadj_300((length(FFTtemp)-CF_index_300+1):end)=FFTtemp((length(FFTtemp)-CF_index_300+1):end); %keep negative freqs
adjSC_300=ifft(FFTadj_300)+1;
SUMCORadj_B_300=real(adjSC_300(1:length(SUMCOR_B)));

%%%%%SUMCOR IFFT 5-64 Hz%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FFTadj_5_64=zeros(size(FFTtemp));
FFTadj_5_64(CF_index_5:CF_index_64)=FFTtemp(CF_index_5:CF_index_64);
FFTadj_5_64((length(FFTtemp)-CF_index_64+1):(length(FFTtemp)-CF_index_5+1))=FFTtemp((length(FFTtemp)-CF_index_64+1):(length(FFTtemp)-CF_index_5+1)); %keep negative freqs
adjSC_5_64=ifft(FFTadj_5_64)+1;
SUMCORadj_B_5_64=real(adjSC_5_64(1:length(SUMCOR_B)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FFTadj_64_300=zeros(size(FFTtemp));
FFTadj_64_300(CF_index_64:CF_index_300)=FFTtemp(CF_index_64:CF_index_300);
FFTadj_64_300((length(FFTtemp)-CF_index_300+1):(length(FFTtemp)-CF_index_64+1))=FFTtemp((length(FFTtemp)-CF_index_300+1):(length(FFTtemp)-CF_index_64+1)); %keep negative freqs
adjSC_64_300=ifft(FFTadj_64_300)+1;
SUMCORadj_B_64_300=real(adjSC_64_300(1:length(SUMCOR_B)));

%% Compute ENVELOPE POWER SPECTRAL DENSITY for Condition B
PSDenv_B=abs(FFTadj);  % use freqVEC to store
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% COLUMN 3: CONDITION=C
%% Compute SAC (C+ and C-) %%%%%%%%%%%%%%%%%%%%
[SAC_C_plus,delays_usec,AvgRate_sps(3,1),xxNspikes] = Library.ShufAutoCorr(ST_C_plus,paramsOUT.DELAYbinwidth_sec,paramsOUT.SCCdur_sec);
[SAC_C_minus,delays_usec,AvgRate_sps(3,2),xxNspikes] = Library.ShufAutoCorr(ST_C_minus,paramsOUT.DELAYbinwidth_sec,paramsOUT.SCCdur_sec);
% SAC limited-data Window CORRECTION
SAC_C_plus=SAC_C_plus+WindowCORRECTION;
SAC_C_minus=SAC_C_minus+WindowCORRECTION;
% smooth and AVG correlograms
SAC_C_plus=trifilt(SAC_C_plus(SACinds),paramsOUT.TriFiltWidthSAC);
SAC_C_minus=trifilt(SAC_C_minus(SACinds),paramsOUT.TriFiltWidthSAC);
SAC_C_avg=(SAC_C_plus+SAC_C_minus)/2;

%% Compute XpAC (C+/C- and C-/C+) %%%%%%%%%%%%%%%%%%%%
[XpAC_C_plus,delays_usec,xxAVGrate,xxNspikes] = Library.ShufCrossCorr({ST_C_plus ST_C_minus},paramsOUT.DELAYbinwidth_sec,paramsOUT.SCCdur_sec);
% SAC limited-data Window CORRECTION
XpAC_C_plus=XpAC_C_plus+WindowCORRECTION;
% smooth and AVG correlograms
XpAC_C_plus=trifilt(XpAC_C_plus(SACinds),paramsOUT.TriFiltWidthSAC);
XpAC_C_minus=fliplr(XpAC_C_plus);
XpAC_C_avg=(XpAC_C_plus+XpAC_C_minus)/2;

%% Compute DIFCOR and SUMCOR for Condition C
DIFCOR_C=trifilt(SAC_C_avg-XpAC_C_avg,paramsOUT.TriFiltWidthDC);  % Diffcorr =  Avg(SAC) - Avg(XpAC)
SUMCOR_C=trifilt((SAC_C_avg+XpAC_C_avg)/2,paramsOUT.TriFiltWidthSC); %Sumcorr = 1/2(Avg(SAC) + Avg(XpAC))

%% Remove TFS artifact (centered at 2*CF) from SUMCOR
% zero out above CF
FFTtemp=fft((SUMCOR_C-1),paramsOUT.Nfft_psd);
[y,CF_index]=min(abs(freqVEC-paramsOUT.SACSCC_CF_Hz)); % use SACSCC_CF_Hz
FFTadj=zeros(size(FFTtemp));
FFTadj(1:CF_index)=FFTtemp(1:CF_index);
FFTadj((length(FFTtemp)-CF_index+1):end)=FFTtemp((length(FFTtemp)-CF_index+1):end); %keep negative freqs
adjSC=ifft(FFTadj)+1;
SUMCORadj_C=real(adjSC(1:length(SUMCOR_C)));

%%%%%SUMCOR IFFT 0-5 Hz%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[y,CF_index_5]=min(abs(freqVEC-5)); % use CF_C
FFTadj_5=zeros(size(FFTtemp));
FFTadj_5(1:CF_index_5)=FFTtemp(1:CF_index_5);
FFTadj_5((length(FFTtemp)-CF_index_5+1):end)=FFTtemp((length(FFTtemp)-CF_index_5+1):end); %keep negative freqs
adjSC_5=ifft(FFTadj_5)+1;
SUMCORadj_C_5=real(adjSC_5(1:length(SUMCOR_C)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[y,CF_index_64]=min(abs(freqVEC-64)); % use CF_A
FFTadj_64=zeros(size(FFTtemp));
FFTadj_64(1:CF_index_64)=FFTtemp(1:CF_index_64);
FFTadj_64((length(FFTtemp)-CF_index_64+1):end)=FFTtemp((length(FFTtemp)-CF_index_64+1):end); %keep negative freqs
adjSC_64=ifft(FFTadj_64)+1;
SUMCORadj_C_64=real(adjSC_64(1:length(SUMCOR_C)));

[y,CF_index_300]=min(abs(freqVEC-300)); % use CF_A
FFTadj_300=zeros(size(FFTtemp));
FFTadj_300(1:CF_index_300)=FFTtemp(1:CF_index_300);
FFTadj_300((length(FFTtemp)-CF_index_300+1):end)=FFTtemp((length(FFTtemp)-CF_index_300+1):end); %keep negative freqs
adjSC_300=ifft(FFTadj_300)+1;
SUMCORadj_C_300=real(adjSC_300(1:length(SUMCOR_C)));

%%%%%SUMCOR IFFT 5-64 Hz%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FFTadj_5_64=zeros(size(FFTtemp));
FFTadj_5_64(CF_index_5:CF_index_64)=FFTtemp(CF_index_5:CF_index_64);
FFTadj_5_64((length(FFTtemp)-CF_index_64+1):(length(FFTtemp)-CF_index_5+1))=FFTtemp((length(FFTtemp)-CF_index_64+1):(length(FFTtemp)-CF_index_5+1)); %keep negative freqs
adjSC_5_64=ifft(FFTadj_5_64)+1;
SUMCORadj_C_5_64=real(adjSC_5_64(1:length(SUMCOR_C)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FFTadj_64_300=zeros(size(FFTtemp));
FFTadj_64_300(CF_index_64:CF_index_300)=FFTtemp(CF_index_64:CF_index_300);
FFTadj_64_300((length(FFTtemp)-CF_index_300+1):(length(FFTtemp)-CF_index_64+1))=FFTtemp((length(FFTtemp)-CF_index_300+1):(length(FFTtemp)-CF_index_64+1)); %keep negative freqs
adjSC_64_300=ifft(FFTadj_64_300)+1;
SUMCORadj_C_64_300=real(adjSC_64_300(1:length(SUMCOR_C)));

%% Compute ENVELOPE POWER SPECTRAL DENSITY for Condition C
PSDenv_C=abs(FFTadj);  % use freqVEC to store
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% COLUMN 4: CONDITION=AB
%% Compute SCC (A+/B+ and A-/B-) %%%%%%%%%%%%%%%%%%%%
[SCC_AB_plus,delays_usec,xxAVGrate,xxNspikes] = Library.ShufCrossCorr({ST_A_plus ST_B_plus},paramsOUT.DELAYbinwidth_sec,paramsOUT.SCCdur_sec);
[SCC_AB_minus,delays_usec,xxAVGrate,xxNspikes] = Library.ShufCrossCorr({ST_A_minus ST_B_minus},paramsOUT.DELAYbinwidth_sec,paramsOUT.SCCdur_sec);
% SAC limited-data Window CORRECTION
SCC_AB_plus=SCC_AB_plus+WindowCORRECTION;
SCC_AB_minus=SCC_AB_minus+WindowCORRECTION;
% smooth and AVG correlograms
SCC_AB_plus=trifilt(SCC_AB_plus(SACinds),paramsOUT.TriFiltWidthSAC);
SCC_AB_minus=trifilt(SCC_AB_minus(SACinds),paramsOUT.TriFiltWidthSAC);
SCC_AB_avg=(SCC_AB_plus+SCC_AB_minus)/2; % That will be plotted in first row of plot

%% Compute XpCC (A+/B- and A-/B+) %%%%%%%%%%%%%%%%%%%%
[XpCC_AB_plus,delays_usec,xxAVGrate,xxNspikes] = Library.ShufCrossCorr({ST_A_plus ST_B_minus},paramsOUT.DELAYbinwidth_sec,paramsOUT.SCCdur_sec);
[XpCC_AB_minus,delays_usec,xxAVGrate,xxNspikes] = Library.ShufCrossCorr({ST_A_minus ST_B_plus},paramsOUT.DELAYbinwidth_sec,paramsOUT.SCCdur_sec);
% SAC limited-data Window CORRECTION
XpCC_AB_plus=XpCC_AB_plus+WindowCORRECTION;
XpCC_AB_minus=XpCC_AB_minus+WindowCORRECTION;
% smooth and AVG correlograms
XpCC_AB_plus=trifilt(XpCC_AB_plus(SACinds),paramsOUT.TriFiltWidthSAC);
XpCC_AB_minus=trifilt(XpCC_AB_minus(SACinds),paramsOUT.TriFiltWidthSAC);
XpCC_AB_avg=(XpCC_AB_plus+XpCC_AB_minus)/2; % that will be plotted in second row of plot

%% Compute DIFCOR and SUMCOR for Condition AB
DIFCOR_AB=trifilt(SCC_AB_avg-XpCC_AB_avg,paramsOUT.TriFiltWidthDC);  % Diffcorr =  Avg(SCC) - Avg(XpCC) that will be plotted in third row of plot
SUMCOR_AB=trifilt((SCC_AB_avg+XpCC_AB_avg)/2,paramsOUT.TriFiltWidthSC); %Sumcorr = 1/2(Avg(SCC) + Avg(XpCC))

%% Remove TFS artifact (centered at 2*CF) from SUMCOR
% zero out above CF
FFTtemp=fft((SUMCOR_AB-1),paramsOUT.Nfft_psd);
[y,CF_index]=min(abs(freqVEC-paramsOUT.SACSCC_CF_Hz)); % use min([CF_A CF_B])
FFTadj=zeros(size(FFTtemp));
FFTadj(1:CF_index)=FFTtemp(1:CF_index);
FFTadj((length(FFTtemp)-CF_index+1):end)=FFTtemp((length(FFTtemp)-CF_index+1):end); %keep negative freqs
adjSC=ifft(FFTadj)+1;
SUMCORadj_AB=real(adjSC(1:length(SUMCOR_AB))); % will be plotted in

%%%%%SUMCOR IFFT 0-5 Hz%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[y,CF_index_5]=min(abs(freqVEC-5)); % use CF_B
FFTadj_5=zeros(size(FFTtemp));
FFTadj_5(1:CF_index_5)=FFTtemp(1:CF_index_5);
FFTadj_5((length(FFTtemp)-CF_index_5+1):end)=FFTtemp((length(FFTtemp)-CF_index_5+1):end); %keep negative freqs
adjSC_5=ifft(FFTadj_5)+1;
SUMCORadj_AB_5=real(adjSC_5(1:length(SUMCOR_AB)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[y,CF_index_64]=min(abs(freqVEC-64)); % use CF_A
FFTadj_64=zeros(size(FFTtemp));
FFTadj_64(1:CF_index_64)=FFTtemp(1:CF_index_64);
FFTadj_64((length(FFTtemp)-CF_index_64+1):end)=FFTtemp((length(FFTtemp)-CF_index_64+1):end); %keep negative freqs
adjSC_64=ifft(FFTadj_64)+1;
SUMCORadj_AB_64=real(adjSC_64(1:length(SUMCOR_AB)));

[y,CF_index_300]=min(abs(freqVEC-300)); % use CF_A
FFTadj_300=zeros(size(FFTtemp));
FFTadj_300(1:CF_index_300)=FFTtemp(1:CF_index_300);
FFTadj_300((length(FFTtemp)-CF_index_300+1):end)=FFTtemp((length(FFTtemp)-CF_index_300+1):end); %keep negative freqs
adjSC_300=ifft(FFTadj_300)+1;
SUMCORadj_AB_300=real(adjSC_300(1:length(SUMCOR_AB)));

%%%%%SUMCOR IFFT 5-64 Hz%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FFTadj_5_64=zeros(size(FFTtemp));
FFTadj_5_64(CF_index_5:CF_index_64)=FFTtemp(CF_index_5:CF_index_64);
FFTadj_5_64((length(FFTtemp)-CF_index_64+1):(length(FFTtemp)-CF_index_5+1))=FFTtemp((length(FFTtemp)-CF_index_64+1):(length(FFTtemp)-CF_index_5+1)); %keep negative freqs
adjSC_5_64=ifft(FFTadj_5_64)+1;
SUMCORadj_AB_5_64=real(adjSC_5_64(1:length(SUMCOR_AB)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FFTadj_64_300=zeros(size(FFTtemp));
FFTadj_64_300(CF_index_64:CF_index_300)=FFTtemp(CF_index_64:CF_index_300);
FFTadj_64_300((length(FFTtemp)-CF_index_300+1):(length(FFTtemp)-CF_index_64+1))=FFTtemp((length(FFTtemp)-CF_index_300+1):(length(FFTtemp)-CF_index_64+1)); %keep negative freqs
adjSC_64_300=ifft(FFTadj_64_300)+1;
SUMCORadj_AB_64_300=real(adjSC_64_300(1:length(SUMCOR_AB)));
    

%% Compute ENVELOPE CROSS SPECTRAL DENSITY for Condition AB
PSDenv_AB=abs(FFTtemp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
%% Store all params, output for return;
% % % SACSCCfunctions: structure with ALL (including plus and minus) functions returned
% % SACSCCs=struct('delays_usec',SACdelays_usec,'SAC_A_plus',SAC_A_plus,'SAC_A_minus',SAC_A_minus,'SAC_A_avg',SAC_A_avg,'SAC_B_plus',SAC_B_plus, ...
% % 	'SAC_B_minus',SAC_B_minus,'SAC_B_avg',SAC_B_avg,'SCC_AB_plus',SCC_AB_plus,'SCC_AB_minus',SCC_AB_minus,'SCC_AB_avg',SCC_AB_avg, ...
% % 	'XpAC_A_plus',XpAC_A_plus,'XpAC_A_minus',XpAC_A_minus,'XpAC_A_avg',XpAC_A_avg,'XpAC_B_plus',XpAC_B_plus,'XpAC_B_minus',XpAC_B_minus, ...
% % 	'XpAC_B_avg',XpAC_B_avg,'XpCC_AB_plus',XpCC_AB_plus,'XpCC_AB_minus',XpCC_AB_minus,'XpCC_AB_avg',XpCC_AB_avg,'SUMCOR_A',SUMCOR_A, ...
% % 	'SUMCOR_B',SUMCOR_B,'SUMCOR_AB',SUMCOR_AB,'DIFCOR_A',DIFCOR_A,'DIFCOR_B',DIFCOR_B,'DIFCOR_AB',DIFCOR_AB,'PSDenv_A',PSDenv_A,'PSDenv_B',PSDenv_B, ...
% % 	'PSDenv_C',PSDenv_C,'freqVEC',freqVEC);
% SACSCCfunctions: structure with ALL useful functions returned
SACSCCs=struct('delays_usec',SACdelays_usec,'SAC_A_avg',SAC_A_avg,'SAC_B_avg',SAC_B_avg,'SAC_C_avg',SAC_C_avg,'SCC_AB_avg',SCC_AB_avg,'XpAC_A_avg',XpAC_A_avg, ...
	'XpAC_B_avg',XpAC_B_avg,'XpAC_C_avg',XpAC_C_avg,'XpCC_AB_avg',XpCC_AB_avg,'SUMCOR_A',SUMCOR_A, ...
	'SUMCOR_B',SUMCOR_B,'SUMCOR_C',SUMCOR_C,'SUMCOR_AB',SUMCOR_AB,'SUMCORadj_A',SUMCORadj_A,'SUMCORadj_B',SUMCORadj_B,'SUMCORadj_C',SUMCORadj_C,'SUMCORadj_AB',SUMCORadj_AB, ...
	'DIFCOR_A',DIFCOR_A,'DIFCOR_B',DIFCOR_B,'DIFCOR_C',DIFCOR_C,'DIFCOR_AB',DIFCOR_AB, ...
    'SUMCORadj_A_64',SUMCORadj_A_64,'SUMCORadj_B_64',SUMCORadj_B_64,'SUMCORadj_C_64',SUMCORadj_C_64,'SUMCORadj_AB_64',SUMCORadj_AB_64, ...
    'SUMCORadj_A_5',SUMCORadj_A_5,'SUMCORadj_B_5',SUMCORadj_B_5,'SUMCORadj_C_5',SUMCORadj_C_5,'SUMCORadj_AB_5',SUMCORadj_AB_5, ...
    'SUMCORadj_A_300',SUMCORadj_A_300,'SUMCORadj_B_300',SUMCORadj_B_300,'SUMCORadj_C_300',SUMCORadj_C_300,'SUMCORadj_AB_300',SUMCORadj_AB_300, ...
    'SUMCORadj_A_64_300',SUMCORadj_A_64_300,'SUMCORadj_B_64_300',SUMCORadj_B_64_300,'SUMCORadj_C_64_300',SUMCORadj_C_64_300,'SUMCORadj_AB_64_300',SUMCORadj_AB_64_300, ...
    'SUMCORadj_A_5_64',SUMCORadj_A_5_64,'SUMCORadj_B_5_64',SUMCORadj_B_5_64,'SUMCORadj_C_5_64',SUMCORadj_C_5_64,'SUMCORadj_AB_5_64',SUMCORadj_AB_5_64,'PSDenv_A',PSDenv_A,'PSDenv_B',PSDenv_B,'PSDenv_C',PSDenv_C,'PSDenv_AB',PSDenv_AB,'PSD_freqVEC',freqVEC,'freqVEC',freqVEC);

return;


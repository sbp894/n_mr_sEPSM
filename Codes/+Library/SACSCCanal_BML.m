function [SACSCCfunctions,SACSCCmetrics,paramsOUT] = SACSCCanal_BML(SpikeTrains,paramsIN,PLOT_15panel)
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

if ~exist('PLOT_15panel','var'),    PLOT_15panel=1;   end

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
stimdur_sec=paramsIN.durA_msec/1000;
paramsOUT.stimdur_sec=stimdur_sec;
% Limit the number of spikes in the SAC/SCC analysis
if ~isfield(paramsOUT,'MAXspikes')
	paramsOUT.MAXspikes=3700;  % If used, windowSTs will cut out extra REPs beyond MAXspikes - used to avoid memory limits
end
% CF: use minimum of all CFs to be most conservative for adjSCpeak
CF_Hz=paramsOUT.CF_A_Hz;  
paramsOUT.SACSCC_CF_Hz=CF_Hz;
%% BOOTSTRAPPING params
paramsOUT.BOOTSTRAP_percentdata=0.80;  % Percent of AN reps to include in each RUN

%% BOOTSTRAPPING - SETUP spike REPS to use for each RUN
clear Nreps
for i=1:1
	for j=1:2
		Nreps(i,j)=length(SpikeTrains{i,j});
		Ninclude(i,j)=ceil(paramsOUT.BOOTSTRAP_percentdata*Nreps(i,j)); % some rounding issues with ceil/floor - have to do it this way
		Nexclude(i,j)=max([1 Nreps(i,j)-Ninclude(i,j)]);  % at least 1 rep excluded to make sure runs are different
		Navgs(i,j)=floor(Nreps(i,j)/Nexclude(i,j));  % Number of RUNS for Averaging CCCs, or PSDs
	end
end
paramsOUT.BOOTSTRAP_Navgs=min(min(Navgs));  % Number of RUNS for Averaging CCCs, or PSDs
for i=1:1
	for j=1:2
		BOOTinds{i,j}=cell(1,paramsOUT.BOOTSTRAP_Navgs);
	end
end

for BOOTrep1=1:paramsOUT.BOOTSTRAP_Navgs
	for i=1:1
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
		paramsOUT.SCCdur_sec=stimdur_sec-paramsOUT.SCC_onsetIGNORE_sec;

		NumReps(1,1)=length(ST_A_plus); NumReps(1,2)=length(ST_A_minus); 
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		disp(sprintf('NUMBER OF SPIKES (A_plus, A_minus) = %d, %d',NumDrivenSpikes(1,1),NumDrivenSpikes(1,2)))
		disp(sprintf('NUMBER OF REPS (A_plus, A_minus) = %d, %d',NumReps(1,1),NumReps(1,2)))

		disp(sprintf('   ... Computing SACs/SCCs for AN spikes: REP %d of %d ...',BOOTrep,paramsOUT.BOOTSTRAP_Navgs))
		[SACSCCs,AvgRate_sps]=SACSCCanal_inBML(ST_A_plus,ST_A_minus,paramsOUT);

        %% RANDOMIZE spike times within same window used above (ignore 1st 50 ms, and cut to shorter of 2 stim)
        STrand_A_plus=randomizeSTs(ST_A_plus,paramsOUT.SCC_onsetIGNORE_sec,stimdur_sec);
        STrand_A_minus=randomizeSTs(ST_A_minus,paramsOUT.SCC_onsetIGNORE_sec,stimdur_sec);
        
        disp(sprintf('   ... Computing SACs/SCCs for RANDOM spikes (for Noise Floor) ...'))
        [SACSCCs_rand,AvgRate_rand_sps]=SACSCCanal_inBML(STrand_A_plus,STrand_A_minus,paramsOUT);

    else	% For LAST loop index, take AVG of all functions for computations
		% AVG everything! Then run rest of code ASIS!!
		% currently it has the last REP already
		for BOOTrep2=1:paramsOUT.BOOTSTRAP_Navgs-1
			SACSCCs.SAC_A_avg=SACSCCs.SAC_A_avg+SACSCCfunctions{BOOTrep2}.SAC_A_avg;
			SACSCCs.XpAC_A_avg=SACSCCs.XpAC_A_avg+SACSCCfunctions{BOOTrep2}.XpAC_A_avg;
			SACSCCs.SUMCOR_A=SACSCCs.SUMCOR_A+SACSCCfunctions{BOOTrep2}.SUMCOR_A;
			SACSCCs.SUMCORadj_A=SACSCCs.SUMCORadj_A+SACSCCfunctions{BOOTrep2}.SUMCORadj_A;
			SACSCCs.DIFCOR_A=SACSCCs.DIFCOR_A+SACSCCfunctions{BOOTrep2}.DIFCOR_A;
			SACSCCs.PSDenv_A=SACSCCs.PSDenv_A+SACSCCfunctions{BOOTrep2}.PSDenv_A;
			% RAND
			SACSCCs_rand.PSDenv_A=SACSCCs_rand.PSDenv_A+SACSCCfunctions{BOOTrep2}.rand.PSDenv_A;
		end
		SACSCCs.SAC_A_avg=SACSCCs.SAC_A_avg/paramsOUT.BOOTSTRAP_Navgs;
		SACSCCs.XpAC_A_avg=SACSCCs.XpAC_A_avg/paramsOUT.BOOTSTRAP_Navgs;
		SACSCCs.SUMCOR_A=SACSCCs.SUMCOR_A/paramsOUT.BOOTSTRAP_Navgs;
		SACSCCs.SUMCORadj_A=SACSCCs.SUMCORadj_A/paramsOUT.BOOTSTRAP_Navgs;
		SACSCCs.DIFCOR_A=SACSCCs.DIFCOR_A/paramsOUT.BOOTSTRAP_Navgs;
		%%%%%%%%%%%%%%%%%%%%%%%
		% PSDs - Best to AVG in PSD domain, not in time domain to reduce variance (e.g., Welch technique)
		SACSCCs.PSDenv_A=SACSCCs.PSDenv_A/paramsOUT.BOOTSTRAP_Navgs;
		SACSCCs_rand.PSDenv_A=SACSCCs_rand.PSDenv_A/paramsOUT.BOOTSTRAP_Navgs;
		
		NumDrivenSpikes=zeros(1,2);
		AvgRate_sps=zeros(1,2);
		for BOOTrep3=1:paramsOUT.BOOTSTRAP_Navgs
			NumDrivenSpikes=NumDrivenSpikes+SACSCCmetrics{BOOTrep3}.NumDrivenSpikes/paramsOUT.BOOTSTRAP_Navgs;
			AvgRate_sps=AvgRate_sps+SACSCCmetrics{BOOTrep3}.AvgRate_sps/paramsOUT.BOOTSTRAP_Navgs;
		end
	end

		
	%%%%%%%%%%%%%%%%%%%%
	%% Compute metrics
	%%%%%%%%%%%%%%%%%%%%

	%% SAC/DC/SC Peak Heights
	%%%%%%%%%
	%% SAC peaks - this is CI from Joris et al 2006 (HR) [don't subtract 1]
	% 1) SACpeak (pure MAX)
	SACSCCmetrics{BOOTrep}.SACpeaks_legend{1}='SACpeak_max';
	SACSCCmetrics{BOOTrep}.SACpeak_A(1)=max(SACSCCs.SAC_A_avg);
	% 2) SACpeak (0 delay)
	SACSCCmetrics{BOOTrep}.SACpeaks_legend{2}='SACpeak_0';
	INDEX_0=find(SACSCCs.delays_usec==0);
	SACSCCmetrics{BOOTrep}.SACpeak_A(2)=SACSCCs.SAC_A_avg(INDEX_0);

	%%%%%%%%%
	%% DC peaks - 
	% 1) DCpeak (pure MAX)
	SACSCCmetrics{BOOTrep}.DCpeaks_legend{1}='DCpeak_max';
	SACSCCmetrics{BOOTrep}.DCpeak_A(1)=max(SACSCCs.DIFCOR_A);
	% 2) DCpeak (0 delay)
	SACSCCmetrics{BOOTrep}.DCpeaks_legend{2}='DCpeak_0';
	SACSCCmetrics{BOOTrep}.DCpeak_A(2)=SACSCCs.DIFCOR_A(INDEX_0);

	%%%%%%%%%
	%% SUMCOR peaks (don't subtract 1; Louage et al 2004)
	% 1) raw peaks
	SACSCCmetrics{BOOTrep}.SCpeaks_legend{1}='raw';
	SACSCCmetrics{BOOTrep}.SCpeaks_A(1)=max(SACSCCs.SUMCOR_A);
	% 4) raw peaks of IFFTraw
	SACSCCmetrics{BOOTrep}.SCpeaks_legend{4}='IFFTraw';
	SACSCCmetrics{BOOTrep}.SCpeaks_A(4)=max(SACSCCs.SUMCORadj_A);
	% 5) IFFTraw (0 delay)
	SACSCCmetrics{BOOTrep}.SCpeaks_legend{5}='IFFTraw_0';
	INDEX_0=find(SACSCCs.delays_usec==0);
	SACSCCmetrics{BOOTrep}.SCpeaks_A(5)=SACSCCs.SUMCORadj_A(INDEX_0);


	%% Book-keeping
	SACSCCfunctions{BOOTrep}=SACSCCs;
    SACSCCfunctions{BOOTrep}.rand.PSDenv_A=SACSCCs_rand.PSDenv_A;

	SACSCCmetrics{BOOTrep}.NumDrivenSpikes=NumDrivenSpikes;
	SACSCCmetrics{BOOTrep}.AvgRate_sps=AvgRate_sps;

end

if PLOT_15panel
	plot_SACSCCanal_BML(SACSCCfunctions,SACSCCmetrics,paramsOUT)
end

return;

function [SACSCCs,AvgRate_sps]=SACSCCanal_inBML(ST_A_plus,ST_A_minus,paramsOUT);

MAXdelay_ind=round(paramsOUT.MAXdelay_sec/paramsOUT.DELAYbinwidth_sec);  % old XLIM=250

%% COLUMN 1: CONDITION=A
%% Compute SAC (A+ and A-) %%%%%%%%%%%%%%%%%%%%
[SAC_A_plus,delays_usec,AvgRate_sps(1,1),xxNspikes] = ShufAutoCorr(ST_A_plus,paramsOUT.DELAYbinwidth_sec,paramsOUT.SCCdur_sec);
[SAC_A_minus,delays_usec,AvgRate_sps(1,2),xxNspikes] = ShufAutoCorr(ST_A_minus,paramsOUT.DELAYbinwidth_sec,paramsOUT.SCCdur_sec);
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
[XpAC_A_plus,delays_usec,xxAVGrate,xxNspikes] = ShufCrossCorr({ST_A_plus ST_A_minus},paramsOUT.DELAYbinwidth_sec,paramsOUT.SCCdur_sec);
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

%% Store ENVELOPE POWER SPECTRAL DENSITY for Condition A
PSDenv_A=abs(FFTadj);  % use freqVEC to store

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
%% Store all params, output for return;
% SACSCCfunctions: structure with ALL useful functions returned
SACSCCs=struct('delays_usec',SACdelays_usec,'SAC_A_avg',SAC_A_avg,'XpAC_A_avg',XpAC_A_avg, ...
	'SUMCOR_A',SUMCOR_A, ...
	'SUMCORadj_A',SUMCORadj_A, ...
	'DIFCOR_A',DIFCOR_A,'PSDenv_A',PSDenv_A,'PSD_freqVEC',freqVEC);

return;


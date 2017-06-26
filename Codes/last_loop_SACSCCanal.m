function [SACSCCfunctions,SACSCCmetrics,paramsOUT]=last_loop_SACSCCanal(SACSCCfunctions,SACSCCmetrics,paramsOUT)

SACSCCs=SACSCCfunctions{1};
SACSCCs_rand=SACSCCfunctions{1}.rand;
paramsOUT.BOOTSTRAP_Navgs=length(SACSCCfunctions);

for BOOTrep2=2:paramsOUT.BOOTSTRAP_Navgs
%     SACSCCs.SAC_A_avg=SACSCCs.SAC_A_avg+SACSCCfunctions{BOOTrep2}.SAC_A_avg;
%     SACSCCs.SAC_B_avg=SACSCCs.SAC_B_avg+SACSCCfunctions{BOOTrep2}.SAC_B_avg;
%     SACSCCs.SAC_C_avg=SACSCCs.SAC_C_avg+SACSCCfunctions{BOOTrep2}.SAC_C_avg;
%     SACSCCs.SCC_AC_avg=SACSCCs.SCC_AC_avg+SACSCCfunctions{BOOTrep2}.SCC_AC_avg;
%     SACSCCs.XpAC_A_avg=SACSCCs.XpAC_A_avg+SACSCCfunctions{BOOTrep2}.XpAC_A_avg;
%     SACSCCs.XpAC_B_avg=SACSCCs.XpAC_B_avg+SACSCCfunctions{BOOTrep2}.XpAC_B_avg;
%     SACSCCs.XpAC_C_avg=SACSCCs.XpAC_C_avg+SACSCCfunctions{BOOTrep2}.XpAC_C_avg;
%     SACSCCs.XpCC_AC_avg=SACSCCs.XpCC_AC_avg+SACSCCfunctions{BOOTrep2}.XpCC_AC_avg;
%     SACSCCs.SUMCOR_A=SACSCCs.SUMCOR_A+SACSCCfunctions{BOOTrep2}.SUMCOR_A;
%     SACSCCs.SUMCOR_B=SACSCCs.SUMCOR_B+SACSCCfunctions{BOOTrep2}.SUMCOR_B;
%     SACSCCs.SUMCOR_C=SACSCCs.SUMCOR_C+SACSCCfunctions{BOOTrep2}.SUMCOR_C;
%     SACSCCs.SUMCOR_AC=SACSCCs.SUMCOR_AC+SACSCCfunctions{BOOTrep2}.SUMCOR_AC;
%     SACSCCs.SUMCORadj_A=SACSCCs.SUMCORadj_A+SACSCCfunctions{BOOTrep2}.SUMCORadj_A;
%     SACSCCs.SUMCORadj_B=SACSCCs.SUMCORadj_B+SACSCCfunctions{BOOTrep2}.SUMCORadj_B;
%     SACSCCs.SUMCORadj_C=SACSCCs.SUMCORadj_C+SACSCCfunctions{BOOTrep2}.SUMCORadj_C;
%     SACSCCs.SUMCORadj_AC=SACSCCs.SUMCORadj_AC+SACSCCfunctions{BOOTrep2}.SUMCORadj_AC;
%     SACSCCs.DIFCOR_A=SACSCCs.DIFCOR_A+SACSCCfunctions{BOOTrep2}.DIFCOR_A;
%     SACSCCs.DIFCOR_B=SACSCCs.DIFCOR_B+SACSCCfunctions{BOOTrep2}.DIFCOR_B;
%     SACSCCs.DIFCOR_C=SACSCCs.DIFCOR_C+SACSCCfunctions{BOOTrep2}.DIFCOR_C;
%     SACSCCs.DIFCOR_AC=SACSCCs.DIFCOR_AC+SACSCCfunctions{BOOTrep2}.DIFCOR_AC;
    %% PSDenv
    SACSCCs.PSDenv_A=SACSCCs.PSDenv_A+SACSCCfunctions{BOOTrep2}.PSDenv_A;
    SACSCCs.PSDenv_B=SACSCCs.PSDenv_B+SACSCCfunctions{BOOTrep2}.PSDenv_B;
    SACSCCs.PSDenv_C=SACSCCs.PSDenv_C+SACSCCfunctions{BOOTrep2}.PSDenv_C;
%     SACSCCs.PSDenv_AC=SACSCCs.PSDenv_AC+SACSCCfunctions{BOOTrep2}.PSDenv_AC;
    % RAND
    SACSCCs_rand.PSDenv_A=SACSCCs_rand.PSDenv_A+SACSCCfunctions{BOOTrep2}.rand.PSDenv_A;
    SACSCCs_rand.PSDenv_B=SACSCCs_rand.PSDenv_B+SACSCCfunctions{BOOTrep2}.rand.PSDenv_B;
    SACSCCs_rand.PSDenv_C=SACSCCs_rand.PSDenv_C+SACSCCfunctions{BOOTrep2}.rand.PSDenv_C;
%     SACSCCs_rand.PSDenv_AC=SACSCCs_rand.PSDenv_AC+SACSCCfunctions{BOOTrep2}.rand.PSDenv_AC;
    %% PSDtfs
%     SACSCCs.PSDtfs_A=SACSCCs.PSDtfs_A+SACSCCfunctions{BOOTrep2}.PSDtfs_A;
%     SACSCCs.PSDtfs_B=SACSCCs.PSDtfs_B+SACSCCfunctions{BOOTrep2}.PSDtfs_B;
%     SACSCCs.PSDtfs_C=SACSCCs.PSDtfs_C+SACSCCfunctions{BOOTrep2}.PSDtfs_C;
%     SACSCCs.PSDtfs_AC=SACSCCs.PSDtfs_AC+SACSCCfunctions{BOOTrep2}.PSDtfs_AC;
    % RAND
%     SACSCCs_rand.PSDtfs_A=SACSCCs_rand.PSDtfs_A+SACSCCfunctions{BOOTrep2}.rand.PSDtfs_A;
%     SACSCCs_rand.PSDtfs_B=SACSCCs_rand.PSDtfs_B+SACSCCfunctions{BOOTrep2}.rand.PSDtfs_B;
%     SACSCCs_rand.PSDtfs_C=SACSCCs_rand.PSDtfs_C+SACSCCfunctions{BOOTrep2}.rand.PSDtfs_C;
%     SACSCCs_rand.PSDtfs_AC=SACSCCs_rand.PSDtfs_AC+SACSCCfunctions{BOOTrep2}.rand.PSDtfs_AC;
    
end
% SACSCCs.SAC_A_avg=SACSCCs.SAC_A_avg/paramsOUT.BOOTSTRAP_Navgs;
% SACSCCs.SAC_B_avg=SACSCCs.SAC_B_avg/paramsOUT.BOOTSTRAP_Navgs;
% SACSCCs.SAC_C_avg=SACSCCs.SAC_C_avg/paramsOUT.BOOTSTRAP_Navgs;
% SACSCCs.SCC_AC_avg=SACSCCs.SCC_AC_avg/paramsOUT.BOOTSTRAP_Navgs;
% SACSCCs.XpAC_A_avg=SACSCCs.XpAC_A_avg/paramsOUT.BOOTSTRAP_Navgs;
% SACSCCs.XpAC_B_avg=SACSCCs.XpAC_B_avg/paramsOUT.BOOTSTRAP_Navgs;
% SACSCCs.XpAC_C_avg=SACSCCs.XpAC_C_avg/paramsOUT.BOOTSTRAP_Navgs;
% SACSCCs.XpCC_AC_avg=SACSCCs.XpCC_AC_avg/paramsOUT.BOOTSTRAP_Navgs;
% SACSCCs.SUMCOR_A=SACSCCs.SUMCOR_A/paramsOUT.BOOTSTRAP_Navgs;
% SACSCCs.SUMCOR_B=SACSCCs.SUMCOR_B/paramsOUT.BOOTSTRAP_Navgs;
% SACSCCs.SUMCOR_C=SACSCCs.SUMCOR_C/paramsOUT.BOOTSTRAP_Navgs;
% SACSCCs.SUMCOR_AC=SACSCCs.SUMCOR_AC/paramsOUT.BOOTSTRAP_Navgs;
% SACSCCs.SUMCORadj_A=SACSCCs.SUMCORadj_A/paramsOUT.BOOTSTRAP_Navgs;
% SACSCCs.SUMCORadj_B=SACSCCs.SUMCORadj_B/paramsOUT.BOOTSTRAP_Navgs;
% SACSCCs.SUMCORadj_C=SACSCCs.SUMCORadj_C/paramsOUT.BOOTSTRAP_Navgs;
% SACSCCs.SUMCORadj_AC=SACSCCs.SUMCORadj_AC/paramsOUT.BOOTSTRAP_Navgs;
% SACSCCs.DIFCOR_A=SACSCCs.DIFCOR_A/paramsOUT.BOOTSTRAP_Navgs;
% SACSCCs.DIFCOR_B=SACSCCs.DIFCOR_B/paramsOUT.BOOTSTRAP_Navgs;
% SACSCCs.DIFCOR_C=SACSCCs.DIFCOR_C/paramsOUT.BOOTSTRAP_Navgs;
% SACSCCs.DIFCOR_AC=SACSCCs.DIFCOR_AC/paramsOUT.BOOTSTRAP_Navgs;
%%%%%%%%%%%%%%%%%%%%%%%
% PSDs - Best to AVG in PSD domain, not in time domain to reduce variance (e.g., Welch technique)
SACSCCs.PSDenv_A=SACSCCs.PSDenv_A/paramsOUT.BOOTSTRAP_Navgs;
SACSCCs.PSDenv_B=SACSCCs.PSDenv_B/paramsOUT.BOOTSTRAP_Navgs;
SACSCCs.PSDenv_C=SACSCCs.PSDenv_C/paramsOUT.BOOTSTRAP_Navgs;
% SACSCCs.PSDenv_AC=SACSCCs.PSDenv_AC/paramsOUT.BOOTSTRAP_Navgs;
SACSCCs_rand.PSDenv_A=SACSCCs_rand.PSDenv_A/paramsOUT.BOOTSTRAP_Navgs;
SACSCCs_rand.PSDenv_B=SACSCCs_rand.PSDenv_B/paramsOUT.BOOTSTRAP_Navgs;
SACSCCs_rand.PSDenv_C=SACSCCs_rand.PSDenv_C/paramsOUT.BOOTSTRAP_Navgs;
% SACSCCs_rand.PSDenv_AC=SACSCCs_rand.PSDenv_AC/paramsOUT.BOOTSTRAP_Navgs;
% SACSCCs.PSDtfs_A=SACSCCs.PSDtfs_A/paramsOUT.BOOTSTRAP_Navgs;
% SACSCCs.PSDtfs_B=SACSCCs.PSDtfs_B/paramsOUT.BOOTSTRAP_Navgs;
% SACSCCs.PSDtfs_C=SACSCCs.PSDtfs_C/paramsOUT.BOOTSTRAP_Navgs;
% SACSCCs.PSDtfs_AC=SACSCCs.PSDtfs_AC/paramsOUT.BOOTSTRAP_Navgs;
% SACSCCs_rand.PSDtfs_A=SACSCCs_rand.PSDtfs_A/paramsOUT.BOOTSTRAP_Navgs;
% SACSCCs_rand.PSDtfs_B=SACSCCs_rand.PSDtfs_B/paramsOUT.BOOTSTRAP_Navgs;
% SACSCCs_rand.PSDtfs_C=SACSCCs_rand.PSDtfs_C/paramsOUT.BOOTSTRAP_Navgs;
% SACSCCs_rand.PSDtfs_AC=SACSCCs_rand.PSDtfs_AC/paramsOUT.BOOTSTRAP_Navgs;
% 		%%%%%%%%%%%%%%%%%%%%%%%
% 		%%% Recompute PSDs/CSD based on AVG SCs
% 		%%%%%%%%%%%%%%%%%%%%%%%
% 		SACSCCs.PSDenv_A=abs(fft((SACSCCs.SUMCOR_A-1),paramsOUT.Nfft_psd));
% 		SACSCCs.PSDenv_B=abs(fft((SACSCCs.SUMCOR_B-1),paramsOUT.Nfft_psd));
% 		SACSCCs.PSDenv_C=abs(fft((SACSCCs.SUMCOR_AC-1),paramsOUT.Nfft_psd));

% % % % NumDrivenSpikes=zeros(3,2);
% % % % AvgRate_sps=zeros(3,2);
% % % % for BOOTrep3=1:paramsOUT.BOOTSTRAP_Navgs
% % % %     NumDrivenSpikes=NumDrivenSpikes+SACSCCmetrics{BOOTrep3}.NumDrivenSpikes/paramsOUT.BOOTSTRAP_Navgs;
% % % %     AvgRate_sps=AvgRate_sps+SACSCCmetrics{BOOTrep3}.AvgRate_sps/paramsOUT.BOOTSTRAP_Navgs;
% % % % end

%% PSD/CSD summations - Compute summed energy in ENVELOPE spectral (densities over various frequency ranges)
% create complete list of LH frequencies over which sums are needed for
% SCpeak_adjusted and/or CCCenv
% % % % CF_Hz=min([paramsOUT.CF_A_Hz paramsOUT.CF_B_Hz paramsOUT.CF_C_Hz]);  
% % % % NYQ_Hz=0.5*(1/paramsOUT.DELAYbinwidth_sec);
% % % % % default list
% % % % default_PSD_LHfreqs_Hz=[0 NYQ_Hz; 0 CF_Hz; 0 300; 10 300; 10 CF_Hz; 0 50; 0 100; 0 64; 0 150];
% % % % % add user-passed list, and create unique list
% % % % if isfield(paramsOUT,'PSD_LHfreqs_Hz')
% % % %     paramsOUT.PSD_LHfreqs_Hz=unique([default_PSD_LHfreqs_Hz; paramsOUT.PSD_LHfreqs_Hz],'rows');
% % % % else
% % % %     paramsOUT.PSD_LHfreqs_Hz=default_PSD_LHfreqs_Hz;
% % % % end

% % % % % find INDs in freqVEC_Hz for relevant cutoffs
% % % % PSD_LHfreqs_inds=zeros(size(paramsOUT.PSD_LHfreqs_Hz,1),2);
% % % % for i=1:size(paramsOUT.PSD_LHfreqs_Hz,1)
% % % %     [~,PSD_LHfreqs_inds(i,1)]=min(abs(SACSCCs.freqVEC-paramsOUT.PSD_LHfreqs_Hz(i,1)));
% % % %     [~,PSD_LHfreqs_inds(i,2)]=min(abs(SACSCCs.freqVEC-paramsOUT.PSD_LHfreqs_Hz(i,2)));
% % % % end
% % % % 
% % % % % Compute all sums for all spectral densities
% % % % for i=1:size(paramsOUT.PSD_LHfreqs_Hz,1)
% % % %     SACSCCmetrics{end}.sums.sumPSD_A(i) = sum(SACSCCs.PSDenv_A(PSD_LHfreqs_inds(i,1):PSD_LHfreqs_inds(i,2)));
% % % %     if ~isnan(NumDrivenSpikes(2,1))
% % % %         SACSCCmetrics{end}.sums.sumPSD_B(i) = sum(SACSCCs.PSDenv_B(PSD_LHfreqs_inds(i,1):PSD_LHfreqs_inds(i,2)));
% % % %         SACSCCmetrics{end}.sums.sumCSD_AC(i) = sum(SACSCCs.PSDenv_C(PSD_LHfreqs_inds(i,1):PSD_LHfreqs_inds(i,2)));
% % % %         SACSCCmetrics{end}.sums.sumPSDrand_A(i) = sum(SACSCCs_rand.PSDenv_A(PSD_LHfreqs_inds(i,1):PSD_LHfreqs_inds(i,2)));
% % % %         SACSCCmetrics{end}.sums.sumPSDrand_B(i) = sum(SACSCCs_rand.PSDenv_B(PSD_LHfreqs_inds(i,1):PSD_LHfreqs_inds(i,2)));
% % % %         SACSCCmetrics{end}.sums.sumCSDrand_AC(i) = sum(SACSCCs_rand.PSDenv_C(PSD_LHfreqs_inds(i,1):PSD_LHfreqs_inds(i,2)));
% % % %     end
% % % % end
% % % % SACSCCmetrics{end}.sums.PSD_LHfreqs_Hz=paramsOUT.PSD_LHfreqs_Hz;

% TFS
% % % % NYQ_Hz=0.5*(1/paramsOUT.DELAYbinwidth_sec);
% % % % % default list
% % % % %     default_PSD_LHfreqs_Hz=[0 NYQ_Hz];
% % % % % find INDs in freqVEC_Hz for relevant cutoffs
% % % % for i=1:size(paramsOUT.PSD_LHfreqs_Hz,1)
% % % %     [~,PSD_LHfreqs_inds(i,1)]=min(abs(SACSCCs.freqVEC-paramsOUT.PSD_LHfreqs_Hz(i,1)));
% % % %     [~,PSD_LHfreqs_inds(i,2)]=min(abs(SACSCCs.freqVEC-paramsOUT.PSD_LHfreqs_Hz(i,2)));
% % % % end
% % % % % Compute all sums for all spectral densities
% % % % for i=1:size(paramsOUT.PSD_LHfreqs_Hz,1)
% % % % %     SACSCCmetrics{end}.sums.sumPSDtfs_A(i) = sum(SACSCCs.PSDtfs_A(PSD_LHfreqs_inds(i,1):PSD_LHfreqs_inds(i,2)));
% % % %     if ~isnan(NumDrivenSpikes(2,1))
% % % % %         SACSCCmetrics{end}.sums.sumPSDtfs_B(i) = sum(SACSCCs.PSDtfs_B(PSD_LHfreqs_inds(i,1):PSD_LHfreqs_inds(i,2)));
% % % % %         SACSCCmetrics{end}.sums.sumCSDtfs_C(i) = sum(SACSCCs.PSDtfs_C(PSD_LHfreqs_inds(i,1):PSD_LHfreqs_inds(i,2)));
% % % % %         SACSCCmetrics{end}.sums.sumPSDrandTfs_A(i) = sum(SACSCCs_rand.PSDtfs_A(PSD_LHfreqs_inds(i,1):PSD_LHfreqs_inds(i,2)));
% % % % %         SACSCCmetrics{end}.sums.sumPSDrandTfs_B(i) = sum(SACSCCs_rand.PSDtfs_B(PSD_LHfreqs_inds(i,1):PSD_LHfreqs_inds(i,2)));
% % % % %         SACSCCmetrics{end}.sums.sumCSDrandTfs_AC(i) = sum(SACSCCs_rand.PSDtfs_C(PSD_LHfreqs_inds(i,1):PSD_LHfreqs_inds(i,2)));
% % % %     end
% % % % end
% % % % SACSCCmetrics{end}.sums.PSD_LHfreqsTfs_Hz=paramsOUT.PSD_LHfreqs_Hz;


%%%%%%%%%%%%%%%%%%%%
%% Compute metrics
%%%%%%%%%%%%%%%%%%%%

%% Characteristic Delays
%% NEED TO IMPROVE - find largest, unless there is a "close-to-largest" that
%% is closer to zero!
% probably use some criterion for 2nd largest peak (if with 5%)???
% SACSCCmetrics{end}.CDscc_usec=Library.findCD_SCC(SACSCCs.SCC_AC_avg,SACSCCs.delays_usec);
% SACSCCmetrics{end}.CDenv_usec=Library.findCD_SCC(SACSCCs.SUMCORadj_AC,SACSCCs.delays_usec);  % use IFFTadjusted SUMCOR!
% SACSCCmetrics{end}.CDtfs_usec=Library.findCD_SCC(SACSCCs.DIFCOR_AC,SACSCCs.delays_usec);


%% SAC/DC/SC Peak Heights
%%%%%%%%%
%% SAC peaks - this is CI from Joris et al 2006 (HR) [don't subtract 1]
% 1) SACpeak (pure MAX)
% SACSCCmetrics{end}.SACpeaks_legend{1}='SACpeak_max';
% SACSCCmetrics{end}.SACpeak_A(1)=max(SACSCCs.SAC_A_avg);
% SACSCCmetrics{end}.SACpeak_B(1)=max(SACSCCs.SAC_B_avg);
% SACSCCmetrics{end}.SACpeak_C(1)=max(SACSCCs.SAC_C_avg);
% % 2) SACpeak (0 delay)
% SACSCCmetrics{end}.SACpeaks_legend{2}='SACpeak_0';
% INDEX_0=find(SACSCCs.delays_usec==0);
% SACSCCmetrics{end}.SACpeak_A(2)=SACSCCs.SAC_A_avg(INDEX_0);
% if ~isnan(NumDrivenSpikes(2,1))
%     SACSCCmetrics{end}.SACpeak_B(2)=SACSCCs.SAC_B_avg(INDEX_0);
%     SACSCCmetrics{end}.SACpeak_C(2)=SACSCCs.SAC_C_avg(INDEX_0);
% end
% % 3) SACpeak (CD delay)
% SACSCCmetrics{end}.SACpeaks_legend{3}='SACpeak_CD';
% INDEX_CD=find(SACSCCs.delays_usec==SACSCCmetrics{end}.CDscc_usec);
% SACSCCmetrics{end}.SACpeak_A(3)=SACSCCs.SAC_A_avg(INDEX_CD);
% if ~isnan(NumDrivenSpikes(2,1))
%     SACSCCmetrics{end}.SACpeak_B(3)=SACSCCs.SAC_B_avg(INDEX_CD);
%     SACSCCmetrics{end}.SCCpeak_AC(3)=SACSCCs.SCC_AC_avg(INDEX_CD);
% end
%% User-passed list of delays to compute
% if isfield(paramsOUT,'UserDelays_usec')
%     for i=1:length(paramsOUT.UserDelays_usec)
%         % 4+) SACpeak (user delays)
%         SACSCCmetrics{end}.SACpeaks_legend{end+1}=sprintf('SACpeak_%d',round(paramsOUT.UserDelays_usec(i)));
%         [~,INDEX_user]=min(abs(SACSCCs.delays_usec-paramsOUT.UserDelays_usec(i))); % find closest delay
%         SACSCCmetrics{end}.SACpeak_A(end+1)=SACSCCs.SAC_A_avg(INDEX_user);
%         if ~isnan(NumDrivenSpikes(2,1))
%             SACSCCmetrics{end}.SACpeak_B(end+1)=SACSCCs.SAC_B_avg(INDEX_user);
%             SACSCCmetrics{end}.SCCpeak_AC(end+1)=SACSCCs.SCC_AC_avg(INDEX_user);
%         end
%     end
% end

%%%%%%%%%
%% DC peaks -
% 1) DCpeak (pure MAX)
% SACSCCmetrics{end}.DCpeaks_legend{1}='DCpeak_max';
% SACSCCmetrics{end}.DCpeak_A(1)=max(SACSCCs.DIFCOR_A);
% if ~isnan(NumDrivenSpikes(2,1))
%     SACSCCmetrics{end}.DCpeak_B(1)=max(SACSCCs.DIFCOR_B);
%     SACSCCmetrics{end}.DCpeak_C(1)=max(SACSCCs.DIFCOR_C);
% end
% 2) DCpeak (0 delay)
% SACSCCmetrics{end}.DCpeaks_legend{2}='DCpeak_0';
% SACSCCmetrics{end}.DCpeak_A(2)=SACSCCs.DIFCOR_A(INDEX_0);
% if ~isnan(NumDrivenSpikes(2,1))
%     SACSCCmetrics{end}.DCpeak_B(2)=SACSCCs.DIFCOR_B(INDEX_0);
%     SACSCCmetrics{end}.DCpeak_C(2)=SACSCCs.DIFCOR_C(INDEX_0);
%     SACSCCmetrics{end}.DCpeak_AC(2)=SACSCCs.DIFCOR_AC(INDEX_0);
% end
% % 3) DCpeak (CD delay)
% SACSCCmetrics{end}.DCpeaks_legend{3}='DCpeak_CD';
% INDEX_CD=find(SACSCCs.delays_usec==SACSCCmetrics{end}.CDtfs_usec);
% SACSCCmetrics{end}.DCpeak_A(3)=SACSCCs.DIFCOR_A(INDEX_CD);
% if ~isnan(NumDrivenSpikes(2,1))
%     SACSCCmetrics{end}.DCpeak_B(3)=SACSCCs.DIFCOR_B(INDEX_CD);
%     SACSCCmetrics{end}.DCpeak_AC(3)=SACSCCs.DIFCOR_AC(INDEX_CD);
% end
%% User-passed list of delays to compute
% if isfield(paramsOUT,'UserDelays_usec')
%     for i=1:length(paramsOUT.UserDelays_usec)
%         % 4+) DCpeak (user delays)
%         SACSCCmetrics{end}.DCpeaks_legend{end+1}=sprintf('DCpeak_%d',round(paramsOUT.UserDelays_usec(i)));
%         [~,INDEX_user]=min(abs(SACSCCs.delays_usec-paramsOUT.UserDelays_usec(i))); % find closest delay
%         SACSCCmetrics{end}.DCpeak_A(end+1)=SACSCCs.DIFCOR_A(INDEX_user);
%         if ~isnan(NumDrivenSpikes(2,1))
%             SACSCCmetrics{end}.DCpeak_B(end+1)=SACSCCs.DIFCOR_B(INDEX_user);
%             SACSCCmetrics{end}.DCpeak_AC(end+1)=SACSCCs.DIFCOR_AC(INDEX_user);
%         end
%     end
% end

%%%%%%%%%
%% SUMCOR peaks (don't subtract 1; Louage et al 2004)
% 1) raw peaks
% SACSCCmetrics{end}.SCpeaks_legend{1}='raw';
% SACSCCmetrics{end}.SCpeaks_A(1)=max(SACSCCs.SUMCOR_A);
% if ~isnan(NumDrivenSpikes(2,1))
%     SACSCCmetrics{end}.SCpeaks_B(1)=max(SACSCCs.SUMCOR_B);
%     SACSCCmetrics{end}.SCpeaks_C(1)=max(SACSCCs.SUMCOR_C);
% end
% % 2) Adjusted SCpeak (Eq. 2): Adj SCpeak = rawSCpeak*sum[0,CF]/sum[0,Fs/2]
% SACSCCmetrics{end}.SCpeaks_legend{2}='adj: 0-CF';
% CFsum_index=find((paramsOUT.PSD_LHfreqs_Hz(:,1)==0)&(paramsOUT.PSD_LHfreqs_Hz(:,2)==CF_Hz));
% NYQsum_index=find((paramsOUT.PSD_LHfreqs_Hz(:,1)==0)&(paramsOUT.PSD_LHfreqs_Hz(:,2)==NYQ_Hz));
% SACSCCmetrics{end}.SCpeaks_A(2)=1+(max(SACSCCs.SUMCOR_A)-1)* ...
%     SACSCCmetrics{end}.sums.sumPSD_A(CFsum_index)/SACSCCmetrics{end}.sums.sumPSD_A(NYQsum_index);
% if ~isnan(NumDrivenSpikes(2,1))
%     SACSCCmetrics{end}.SCpeaks_B(2)=1+(max(SACSCCs.SUMCOR_B)-1)* ...
%         SACSCCmetrics{end}.sums.sumPSD_B(CFsum_index)/SACSCCmetrics{end}.sums.sumPSD_B(NYQsum_index);
%     SACSCCmetrics{end}.SCpeaks_AC(2)=1+(max(SACSCCs.SUMCOR_AC)-1)* ...
%         SACSCCmetrics{end}.sums.sumCSD_AC(CFsum_index)/SACSCCmetrics{end}.sums.sumCSD_AC(NYQsum_index);
% end
% % 3) alternative adjusted SC peak height [0,300]/[0,Fs/2]
% % just for FIG 2 in NM paper
% SACSCCmetrics{end}.SCpeaks_legend{3}='adj: 0-300';
% ENVsum_index=find((paramsOUT.PSD_LHfreqs_Hz(:,1)==0)&(paramsOUT.PSD_LHfreqs_Hz(:,2)==300));
% SACSCCmetrics{end}.SCpeaks_A(3)=1+(max(SACSCCs.SUMCOR_A)-1)* ...
%     SACSCCmetrics{end}.sums.sumPSD_A(ENVsum_index)/SACSCCmetrics{end}.sums.sumPSD_A(NYQsum_index);
% if ~isnan(NumDrivenSpikes(2,1))
%     SACSCCmetrics{end}.SCpeaks_B(3)=1+(max(SACSCCs.SUMCOR_B)-1)* ...
%         SACSCCmetrics{end}.sums.sumPSD_B(ENVsum_index)/SACSCCmetrics{end}.sums.sumPSD_B(NYQsum_index);
%     SACSCCmetrics{end}.SCpeaks_AC(3)=1+(max(SACSCCs.SUMCOR_AC)-1)* ...
%         SACSCCmetrics{end}.sums.sumCSD_AC(ENVsum_index)/SACSCCmetrics{end}.sums.sumCSD_AC(NYQsum_index);
% end
% % 4) raw peaks of IFFTraw
% SACSCCmetrics{end}.SCpeaks_legend{4}='IFFTraw';
% SACSCCmetrics{end}.SCpeaks_A(4)=max(SACSCCs.SUMCORadj_A);
% if ~isnan(NumDrivenSpikes(2,1))
%     SACSCCmetrics{end}.SCpeaks_B(4)=max(SACSCCs.SUMCORadj_B);
%     SACSCCmetrics{end}.SCpeaks_C(4)=max(SACSCCs.SUMCORadj_C);
%     SACSCCmetrics{end}.SCpeaks_AC(4)=max(SACSCCs.SUMCORadj_AC);
% end
% % 5) IFFTraw (0 delay)
% SACSCCmetrics{end}.SCpeaks_legend{5}='IFFTraw_0';
% INDEX_0=find(SACSCCs.delays_usec==0);
% SACSCCmetrics{end}.SCpeaks_A(5)=SACSCCs.SUMCORadj_A(INDEX_0);
% if ~isnan(NumDrivenSpikes(2,1))
%     SACSCCmetrics{end}.SCpeaks_B(5)=SACSCCs.SUMCORadj_B(INDEX_0);
%     SACSCCmetrics{end}.SCpeaks_C(5)=SACSCCs.SUMCORadj_C(INDEX_0);
%     SACSCCmetrics{end}.SCpeaks_AC(5)=SACSCCs.SUMCORadj_AC(INDEX_0);
% end
% %%%Add SUMCOR peaks from 0-64, 0-5, 5-64, 64-300 and 0-300 here%%%%
% 
% SACSCCmetrics{end}.SCpeaks_legend{end+1}='IFFTraw_5@0';
% INDEX_0=find(SACSCCs.delays_usec==0);
% SACSCCmetrics{end}.SCpeaks_A(end+1)=SACSCCs.SUMCORadj_A_5(INDEX_0);
% if ~isnan(NumDrivenSpikes(2,1))
%     SACSCCmetrics{end}.SCpeaks_B(end+1)=SACSCCs.SUMCORadj_B_5(INDEX_0);
%     SACSCCmetrics{end}.SCpeaks_C(end+1)=SACSCCs.SUMCORadj_C_5(INDEX_0);
%     SACSCCmetrics{end}.SCpeaks_AC(end+1)=SACSCCs.SUMCORadj_AC_5(INDEX_0);
% end
% 
% SACSCCmetrics{end}.SCpeaks_legend{end+1}='IFFTraw_64@0';
% INDEX_0=find(SACSCCs.delays_usec==0);
% SACSCCmetrics{end}.SCpeaks_A(end+1)=SACSCCs.SUMCORadj_A_64(INDEX_0);
% if ~isnan(NumDrivenSpikes(2,1))
%     SACSCCmetrics{end}.SCpeaks_B(end+1)=SACSCCs.SUMCORadj_B_64(INDEX_0);
%     SACSCCmetrics{end}.SCpeaks_C(end+1)=SACSCCs.SUMCORadj_C_64(INDEX_0);
%     SACSCCmetrics{end}.SCpeaks_AC(end+1)=SACSCCs.SUMCORadj_AC_64(INDEX_0);
% end
% 
% SACSCCmetrics{end}.SCpeaks_legend{end+1}='IFFTraw_300@0';
% INDEX_0=find(SACSCCs.delays_usec==0);
% SACSCCmetrics{end}.SCpeaks_A(end+1)=SACSCCs.SUMCORadj_A_300(INDEX_0);
% if ~isnan(NumDrivenSpikes(2,1))
%     SACSCCmetrics{end}.SCpeaks_B(end+1)=SACSCCs.SUMCORadj_B_300(INDEX_0);
%     SACSCCmetrics{end}.SCpeaks_C(end+1)=SACSCCs.SUMCORadj_C_300(INDEX_0);
%     SACSCCmetrics{end}.SCpeaks_AC(end+1)=SACSCCs.SUMCORadj_AC_300(INDEX_0);
% end
% 
% SACSCCmetrics{end}.SCpeaks_legend{end+1}='IFFTraw_5_64@0';
% INDEX_0=find(SACSCCs.delays_usec==0);
% SACSCCmetrics{end}.SCpeaks_A(end+1)=SACSCCs.SUMCORadj_A_5_64(INDEX_0);
% if ~isnan(NumDrivenSpikes(2,1))
%     SACSCCmetrics{end}.SCpeaks_B(end+1)=SACSCCs.SUMCORadj_B_5_64(INDEX_0);
%     SACSCCmetrics{end}.SCpeaks_C(end+1)=SACSCCs.SUMCORadj_C_5_64(INDEX_0);
%     SACSCCmetrics{end}.SCpeaks_AC(end+1)=SACSCCs.SUMCORadj_AC_5_64(INDEX_0);
% end
% 
% SACSCCmetrics{end}.SCpeaks_legend{end+1}='IFFTraw_64_300@0';
% INDEX_0=find(SACSCCs.delays_usec==0);
% SACSCCmetrics{end}.SCpeaks_A(end+1)=SACSCCs.SUMCORadj_A_64_300(INDEX_0);
% if ~isnan(NumDrivenSpikes(2,1))
%     SACSCCmetrics{end}.SCpeaks_B(end+1)=SACSCCs.SUMCORadj_B_64_300(INDEX_0);
%     SACSCCmetrics{end}.SCpeaks_C(end+1)=SACSCCs.SUMCORadj_C_64_300(INDEX_0);
%     SACSCCmetrics{end}.SCpeaks_AC(end+1)=SACSCCs.SUMCORadj_AC_64_300(INDEX_0);
% end
% % 6) IFFTraw (CD delay)
% SACSCCmetrics{end}.SCpeaks_legend{end+1}='IFFTraw_CD';
% INDEX_CD=find(SACSCCs.delays_usec==SACSCCmetrics{end}.CDenv_usec);
% SACSCCmetrics{end}.SCpeaks_A(end+1)=SACSCCs.SUMCORadj_A(INDEX_CD);
% if ~isnan(NumDrivenSpikes(2,1))
%     SACSCCmetrics{end}.SCpeaks_B(end+1)=SACSCCs.SUMCORadj_B(INDEX_CD);
%     SACSCCmetrics{end}.SCpeaks_C(end+1)=SACSCCs.SUMCORadj_C(INDEX_CD);
%     SACSCCmetrics{end}.SCpeaks_AC(end+1)=SACSCCs.SUMCORadj_AC(INDEX_CD);
% end
% %% User-passed list of delays to compute
% if isfield(paramsOUT,'UserDelays_usec')
%     for i=1:length(paramsOUT.UserDelays_usec)
%         % 7+) IFFTraw (user delays)
%         SACSCCmetrics{end}.SCpeaks_legend{end+1}=sprintf('IFFTraw_%d',round(paramsOUT.UserDelays_usec(i)));
%         [~,INDEX_user]=min(abs(SACSCCs.delays_usec-paramsOUT.UserDelays_usec(i))); % find closest delay
%         SACSCCmetrics{end}.SCpeaks_A(end+1)=SACSCCs.SUMCORadj_A(INDEX_user);
%         if ~isnan(NumDrivenSpikes(2,1))
%             SACSCCmetrics{end}.SCpeaks_B(end+1)=SACSCCs.SUMCORadj_B(INDEX_user);
%             SACSCCmetrics{end}.SCpeaks_AC(end+1)=SACSCCs.SUMCORadj_AC(INDEX_user);
%         end
%     end
% end
% 
% %%%%%%%%%%%%%%%%%
% %% Neural Cross Correlation Coefficients
% %%%%%%%%%%%%%%%%%
% if ~isnan(NumDrivenSpikes(2,1))
%     %% CCCtfs
%     DC_max_index=find(strcmp(SACSCCmetrics{end}.DCpeaks_legend,'DCpeak_max'));
%     if (SACSCCmetrics{end}.DCpeak_A(DC_max_index)>=paramsOUT.minDCpeak_CCCtfs) ...
%             && (SACSCCmetrics{end}.DCpeak_B(DC_max_index)>=paramsOUT.minDCpeak_CCCtfs)
%         % 1) using DCpeak_max
%         SACSCCmetrics{end}.CCCtfs_legend{1}='DCpeak_max';
%         SACSCCmetrics{end}.CCCtfs(1) = SACSCCmetrics{end}.DCpeak_AC(DC_max_index) ...
%             / (sqrt(SACSCCmetrics{end}.DCpeak_A(DC_max_index)*SACSCCmetrics{end}.DCpeak_B(DC_max_index)));
%         % 2) using DCpeak_0
%         DC_0_index=find(strcmp(SACSCCmetrics{end}.DCpeaks_legend,'DCpeak_0'));
%         SACSCCmetrics{end}.CCCtfs_legend{2}='DCpeak_0';
%         SACSCCmetrics{end}.CCCtfs(2) = SACSCCmetrics{end}.DCpeak_AC(DC_0_index) ...
%             / (sqrt(SACSCCmetrics{end}.DCpeak_A(DC_0_index)*SACSCCmetrics{end}.DCpeak_B(DC_0_index)));
%         % 3) using DCpeak_CD
%         DC_CD_index= strcmp(SACSCCmetrics{end}.DCpeaks_legend,'DCpeak_CD');
%         SACSCCmetrics{end}.CCCtfs_legend{3}='DCpeak_CD';
%         % NOTE: peaks for ACFs are taken at 0 delay (by definition),
%         % rather than at specified CCF delay
%         SACSCCmetrics{end}.CCCtfs(3) = SACSCCmetrics{end}.DCpeak_AC(DC_CD_index) ...
%             / (sqrt(SACSCCmetrics{end}.DCpeak_A(DC_0_index)*SACSCCmetrics{end}.DCpeak_B(DC_0_index)));
%         % 4+) User-passed list of delays to compute
%         if isfield(paramsOUT,'UserDelays_usec')
%             for i=1:length(paramsOUT.UserDelays_usec)
%                 DC_user_index= strcmp(SACSCCmetrics{end}.DCpeaks_legend,sprintf('DCpeak_%d',round(paramsOUT.UserDelays_usec(i))));
%                 SACSCCmetrics{end}.CCCtfs_legend{end+1}=sprintf('DCpeak_%d',round(paramsOUT.UserDelays_usec(i)));
%                 % NOTE: peaks for ACFs are taken at 0 delay (by definition),
%                 % rather than at specified CCF delay
%                 SACSCCmetrics{end}.CCCtfs(end+1) = SACSCCmetrics{end}.DCpeak_AC(DC_user_index) ...
%                     / (sqrt(SACSCCmetrics{end}.DCpeak_A(DC_0_index)*SACSCCmetrics{end}.DCpeak_B(DC_0_index)));
%             end
%         end
%         
%     else
%         for i=length(SACSCCmetrics{end}.DCpeaks_legend)
%             SACSCCmetrics{end}.CCCtfs(i) = NaN;
%         end
%     end
%     
%     %% CCCenv
%     % first, go through whole LHfreq list, for subBIAS (subtract bias) and
%     % withBIAS (Don't subtract bias)
%     for i=1:size(paramsOUT.PSD_LHfreqs_Hz,1)
%         if paramsOUT.PSD_LHfreqs_Hz(i,2)==CF_Hz
%             SACSCCmetrics{end}.CCCenvs_legend{(i-1)*2+1}=sprintf('%.f-CF, subBIAS',paramsOUT.PSD_LHfreqs_Hz(i,1));
%             SACSCCmetrics{end}.CCCenvs_legend{(i-1)*2+2}=sprintf('%.f-CF, withBIAS',paramsOUT.PSD_LHfreqs_Hz(i,1));
%         elseif paramsOUT.PSD_LHfreqs_Hz(i,2)==NYQ_Hz
%             SACSCCmetrics{end}.CCCenvs_legend{(i-1)*2+1}=sprintf('%.f-NYQ, subBIAS',paramsOUT.PSD_LHfreqs_Hz(i,1));
%             SACSCCmetrics{end}.CCCenvs_legend{(i-1)*2+2}=sprintf('%.f-NYQ, withBIAS',paramsOUT.PSD_LHfreqs_Hz(i,1));
%         else
%             SACSCCmetrics{end}.CCCenvs_legend{(i-1)*2+1}=sprintf('%.f-%.f, subBIAS',paramsOUT.PSD_LHfreqs_Hz(i,1),paramsOUT.PSD_LHfreqs_Hz(i,2));
%             SACSCCmetrics{end}.CCCenvs_legend{(i-1)*2+2}=sprintf('%.f-%.f, withBIAS',paramsOUT.PSD_LHfreqs_Hz(i,1),paramsOUT.PSD_LHfreqs_Hz(i,2));
%         end
%         % subBIAS
%         SACSCCmetrics{end}.CCCenvs((i-1)*2+1) = max([0 (SACSCCmetrics{end}.sums.sumCSD_AC(i)-SACSCCmetrics{end}.sums.sumCSDrand_AC(i)) ])/ ...
%             sqrt(max([0 (SACSCCmetrics{end}.sums.sumPSD_A(i)-SACSCCmetrics{end}.sums.sumPSDrand_A(i)) ])* ...
%             max([0 (SACSCCmetrics{end}.sums.sumPSD_B(i)-SACSCCmetrics{end}.sums.sumPSDrand_B(i)) ]));
%         % withBIAS
%         SACSCCmetrics{end}.CCCenvs((i-1)*2+2) = SACSCCmetrics{end}.sums.sumCSD_AC(i)/ ...
%             sqrt(SACSCCmetrics{end}.sums.sumPSD_A(i)*SACSCCmetrics{end}.sums.sumPSD_B(i));
%     end
%     % using raw SCpeaks
%     RAWSCindex=find(strcmp(SACSCCmetrics{end}.SCpeaks_legend,'raw'));
%     SACSCCmetrics{end}.CCCenvs_legend{end+1}='rawSC';
%     SACSCCmetrics{end}.CCCenvs(end+1) = (SACSCCmetrics{end}.SCpeaks_AC(RAWSCindex)-1)/ ...
%         (sqrt((SACSCCmetrics{end}.SCpeaks_A(RAWSCindex)-1)*(SACSCCmetrics{end}.SCpeaks_B(RAWSCindex)-1)));
%     % using adj SCpeaks
%     ADJSCindex=find(strcmp(SACSCCmetrics{end}.SCpeaks_legend,'adj: 0-CF'));
%     SACSCCmetrics{end}.CCCenvs_legend{end+1}='adjSC';
%     SACSCCmetrics{end}.CCCenvs(end+1) = (SACSCCmetrics{end}.SCpeaks_AC(ADJSCindex)-1)/ ...
%         (sqrt((SACSCCmetrics{end}.SCpeaks_A(ADJSCindex)-1)*(SACSCCmetrics{end}.SCpeaks_B(ADJSCindex)-1)));
%     % using IFFTraw SCpeaks
%     IFFTSCindex=find(strcmp(SACSCCmetrics{end}.SCpeaks_legend,'IFFTraw'));
%     SACSCCmetrics{end}.CCCenvs_legend{end+1}='IFFTrawSC';
%     SACSCCmetrics{end}.CCCenvs(end+1) = (SACSCCmetrics{end}.SCpeaks_AC(IFFTSCindex)-1)/ ...
%         (sqrt((SACSCCmetrics{end}.SCpeaks_A(IFFTSCindex)-1)*(SACSCCmetrics{end}.SCpeaks_B(IFFTSCindex)-1)));
%     % using IFFT(0 delay) SCpeaks
%     IFFTSC_0_index=find(strcmp(SACSCCmetrics{end}.SCpeaks_legend,'IFFTraw_0'));
%     SACSCCmetrics{end}.CCCenvs_legend{end+1}='IFFTrawSC_0';
%     SACSCCmetrics{end}.CCCenvs(end+1) = (SACSCCmetrics{end}.SCpeaks_AC(IFFTSC_0_index)-1)/ ...
%         (sqrt((SACSCCmetrics{end}.SCpeaks_A(IFFTSC_0_index)-1)*(SACSCCmetrics{end}.SCpeaks_B(IFFTSC_0_index)-1)));
%     % using IFFT(CD delay) SCpeaks
%     IFFTSC_CD_index= strcmp(SACSCCmetrics{end}.SCpeaks_legend,'IFFTraw_CD');
%     SACSCCmetrics{end}.CCCenvs_legend{end+1}='IFFTrawSC_CD';
%     % NOTE: peaks for ACFs are taken at 0 delay (by definition), rather
%     % than at specified CCF delay
%     SACSCCmetrics{end}.CCCenvs(end+1) = (SACSCCmetrics{end}.SCpeaks_AC(IFFTSC_CD_index)-1)/ ...
%         (sqrt((SACSCCmetrics{end}.SCpeaks_A(IFFTSC_0_index)-1)*(SACSCCmetrics{end}.SCpeaks_B(IFFTSC_0_index)-1)));
%     % User-passed list of delays to compute
%     if isfield(paramsOUT,'UserDelays_usec')
%         for i=1:length(paramsOUT.UserDelays_usec)
%             IFFTSC_user_index= strcmp(SACSCCmetrics{end}.SCpeaks_legend,sprintf('IFFTraw_%d',round(paramsOUT.UserDelays_usec(i))));
%             SACSCCmetrics{end}.CCCenvs_legend{end+1}=sprintf('IFFTrawSC_%d',round(paramsOUT.UserDelays_usec(i)));
%             % NOTE: peaks for ACFs are taken at 0 delay (by definition),
%             % rather than at specified CCF delay
%             SACSCCmetrics{end}.CCCenvs(end+1) = (SACSCCmetrics{end}.SCpeaks_AC(IFFTSC_user_index)-1)/ ...
%                 (sqrt((SACSCCmetrics{end}.SCpeaks_A(IFFTSC_0_index)-1)*(SACSCCmetrics{end}.SCpeaks_B(IFFTSC_0_index)-1)));
%         end
%     end
%     
% end

%% Book-keeping
SACSCCfunctions{paramsOUT.BOOTSTRAP_Navgs+1}=SACSCCs;
% if ~isnan(NumDrivenSpikes(2,1))
SACSCCfunctions{paramsOUT.BOOTSTRAP_Navgs+1}.rand.PSDenv_A=SACSCCs_rand.PSDenv_A;
SACSCCfunctions{paramsOUT.BOOTSTRAP_Navgs+1}.rand.PSDenv_B=SACSCCs_rand.PSDenv_B;
SACSCCfunctions{paramsOUT.BOOTSTRAP_Navgs+1}.rand.PSDenv_C=SACSCCs_rand.PSDenv_C;
%     SACSCCfunctions{end}.rand.PSDenv_AC=SACSCCs_rand.PSDenv_AC;
%     SACSCCfunctions{end}.rand.PSDtfs_A=SACSCCs_rand.PSDtfs_A;
%     SACSCCfunctions{end}.rand.PSDtfs_B=SACSCCs_rand.PSDtfs_B;
%     SACSCCfunctions{end}.rand.PSDtfs_C=SACSCCs_rand.PSDtfs_C;
%     SACSCCfunctions{end}.rand.PSDtfs_AC=SACSCCs_rand.PSDtfs_AC;
% end

% SACSCCmetrics{end}.NumDrivenSpikes=NumDrivenSpikes;
% SACSCCmetrics{end}.AvgRate_sps=AvgRate_sps;

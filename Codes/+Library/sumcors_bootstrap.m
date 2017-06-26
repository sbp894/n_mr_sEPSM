function [PSDenv_STRUCT,PSDtfs_STRUCT,PowerMod_STRUCT,PowerTfs_STRUCT] = sumcors_bootstrap...
    (SpikeTrains,paramsIN, resultsDir,resultPostfix, BootstrapLoopMax, BootstrapLoopReport, winCorr0Add1Mul)
PSDenv_STRUCT=[];
PSDtfs_STRUCT=[];
PowerTfs_STRUCT=[];
% This function returns the FFT of the SUMCORs (envelope) for each stimulus derived from the
% spike output of one AN fiber. It also returns the envelope power within a
% certain number of modulation frequency bands.

%%% Input parmeters
% stim_SN: Sentence plus noise
% stim_N: Noise
% CF_kHz: Characteristic frequency of the nerve fiber
% fiberType: 1-LSR, 2-MSR, 3-HSR (Spontaneous rate)
% Cohc: Outer hair cell parameter
% Cihc: Inner hair cell parameter
% OALevel_dBSPL: Stimulus intensity (for both conditions)


%% *MH: HOW TO OPTIMIZE??
% BootstrapLoopMax=ExpControlParams.BootstrapLoopMax;%300;
% BootstrapLoopReport=ExpControlParams.BootstrapLoopReport;
% nPSDs2Avg=ExpControlParams.nPSDs2Avg;
PLOT_15panel=0;
% [SACSCCfunctions,SACSCCmetrics,paramsOUT] = SACSCCanal_SNRenv_parallel(SpikeTrains,paramsIN,PLOT_15panel,resultsDir,resultPostfix,BootstrapLoopMax);
[SACSCCfunctions,SACSCCmetrics,~] = SACSCCanal_SNRenv_parallel(SpikeTrains,paramsIN,PLOT_15panel,resultsDir,resultPostfix,BootstrapLoopMax, winCorr0Add1Mul);

if ~isempty(SACSCCfunctions)
    BootstrapLoopMax=length(SACSCCfunctions)-1;
    nPSDs2Avg=max(1, ceil(BootstrapLoopMax/2));
    % % % % Library.parsave([resultsDir 'sac' filesep 'sac' resultPostfix '.mat'],  SACSCCfunctions, SACSCCmetrics );
    
    %%% PSDenv computed in SACSCCanal_SNRenv and passed within SACSCCfunctions
    % SACSCCfunctions{end} gives the all the average values
% % % %     PSDenv_S=SACSCCfunctions{end}.PSDenv_A;
% % % %     PSDenv_N=SACSCCfunctions{end}.PSDenv_B;
% % % %     PSDenv_SN=SACSCCfunctions{end}.PSDenv_C;
% % % %     PSDenv_S_noisefloor=SACSCCfunctions{end}.rand.PSDenv_A;
% % % %     PSDenv_N_noisefloor=SACSCCfunctions{end}.rand.PSDenv_B;
% % % %     PSDenv_SN_noisefloor=SACSCCfunctions{end}.rand.PSDenv_C;
    PSDfreqVEC_Hz=SACSCCfunctions{end}.PSD_freqVEC;
    
    % tfs
    % % % % tfsCutOff = 1.5;
%     CF_A_kHz=paramsIN.CF_A_Hz/1e3;
    
    % % % % CF_B_kHz=CF_A_kHz;
%     CF_C_kHz=CF_A_kHz;
    
    % if CF_A_kHz < tfsCutOff && CF_B_kHz < tfsCutOff && CF_C_kHz < tfsCutOff
    %     PSDtfs_S=SACSCCfunctions{end}.PSDtfs_A;
    %     PSDtfs_N=SACSCCfunctions{end}.PSDtfs_B;
    %     PSDtfs_SN=SACSCCfunctions{end}.PSDtfs_C;
    %     PSDtfs_S_noisefloor=SACSCCfunctions{end}.rand.PSDtfs_A;
    %     PSDtfs_N_noisefloor=SACSCCfunctions{end}.rand.PSDtfs_B;
    %     PSDtfs_SN_noisefloor=SACSCCfunctions{end}.rand.PSDtfs_C;
    % else
    %     PSDtfs_S=0;
    %     PSDtfs_N=0;
    %     PSDtfs_SN=0;
    %     PSDtfs_S_noisefloor=0;
    %     PSDtfs_N_noisefloor=0;
    %     PSDtfs_SN_noisefloor=0;
    % end
    
    %% *MH - look at spectra
% % %     if paramsIN.plt
% % % %         figure(figHandles.PSDplot);
% % %         clf;
% % %         set (gcf, 'Units', 'normalized', 'Position', [.1 .1 .8 .8]);
% % %         plot(PSDfreqVEC_Hz,PSDenv_S,'b-');
% % %         hold on; plot(PSDfreqVEC_Hz,PSDenv_N,'r-');
% % %         plot(PSDfreqVEC_Hz,PSDenv_SN,'k-');
% % %         xlim([0 64])
% % %         plot(PSDfreqVEC_Hz,PSDenv_S_noisefloor,'b:');
% % %         plot(PSDfreqVEC_Hz,PSDenv_N_noisefloor,'r:');
% % %         plot(PSDfreqVEC_Hz,PSDenv_SN_noisefloor,'k:');
% % %         legend('S','N','SN','NF[S]','NF[N]','NF[SN]')
% % %         xlabel(sprintf('MODULATION FREQUENCY (Hz)'),'FontSize',12,'HorizontalAlignment','center')
% % %         ylabel(sprintf('POWER SPECTRAL DENSITY'),'FontSize',12)
% % %         
% % %         title(sprintf('CF=%.3f kHz; SRtype=%d (1: LSR; 2: MSR; 3: HSR);\nCohc=%.2f; Cihc=%.2f; Overall (Signal) Level: %.1f dB SPL; SNR = %.1f dB', ...
% % %             CF_C_kHz,paramsIN.FiberType,99,99,paramsIN.level,paramsIN.SNR2use_dB),'FontSize',12);
% % %         
% % %         Library.saveFigureAs([resultsDir 'eps' filesep 'eps' resultPostfix '.eps'])
% % %         Library.saveFigureAs([resultsDir 'png' filesep 'psd' resultPostfix '.tiff'])
% % %     end
    % Divide the FFT output into modulation filters
    ModFreqs = [1,2,4,8,16,32,64];
    PowerModS = zeros(BootstrapLoopReport,length(ModFreqs));
    PowerModN = zeros(BootstrapLoopReport,length(ModFreqs));
    PowerModSN = zeros(BootstrapLoopReport,length(ModFreqs));
    PowerModS_noisefloor = zeros(BootstrapLoopReport,length(ModFreqs));
    PowerModN_noisefloor = zeros(BootstrapLoopReport,length(ModFreqs));
    PowerModSN_noisefloor = zeros(BootstrapLoopReport,length(ModFreqs));
    
    %% Look at RAND PSDs as well to get noise floor.
    plt =0;
    qFactor=ones(size(ModFreqs));
    PSDsize=size(SACSCCfunctions{1}.PSDenv_A);
    for bs_var=1:BootstrapLoopReport
        
        ind=randsample(1:BootstrapLoopMax, nPSDs2Avg);
        PSDenv_A=zeros(PSDsize);
        PSDenv_B=zeros(PSDsize);
        PSDenv_C=zeros(PSDsize);
        PSDenv_A_NF=zeros(PSDsize);
        PSDenv_B_NF=zeros(PSDsize);
        PSDenv_C_NF=zeros(PSDsize);
        
        for psd_smooth_var=1:length(ind)
            PSDenv_A=PSDenv_A+SACSCCfunctions{ind(psd_smooth_var)}.PSDenv_A/length(ind);
            PSDenv_B=PSDenv_B+SACSCCfunctions{ind(psd_smooth_var)}.PSDenv_B/length(ind);
            PSDenv_C=PSDenv_C+SACSCCfunctions{ind(psd_smooth_var)}.PSDenv_C/length(ind);
            PSDenv_A_NF=PSDenv_A_NF+SACSCCfunctions{ind(psd_smooth_var)}.rand.PSDenv_A/length(ind);
            PSDenv_B_NF=PSDenv_B_NF+SACSCCfunctions{ind(psd_smooth_var)}.rand.PSDenv_B/length(ind);
            PSDenv_C_NF=PSDenv_C_NF+SACSCCfunctions{ind(psd_smooth_var)}.rand.PSDenv_C/length(ind);
        end
        
        binWeights = modFbank(PSDfreqVEC_Hz,ModFreqs, qFactor, plt);
        
        PowerModS(bs_var,:)=nansum(binWeights.* repmat(PSDenv_A, length(ModFreqs), 1),2);
        PowerModN(bs_var,:) = nansum(binWeights.* repmat(PSDenv_B, length(ModFreqs), 1),2);
        PowerModSN(bs_var,:) = nansum(binWeights.* repmat(PSDenv_C, length(ModFreqs), 1),2);
        PowerModS_noisefloor(bs_var,:) = nansum(binWeights.* repmat(PSDenv_A_NF, length(ModFreqs), 1),2);
        PowerModN_noisefloor(bs_var,:) = nansum(binWeights.* repmat(PSDenv_B_NF, length(ModFreqs), 1),2);
        PowerModSN_noisefloor(bs_var,:) = nansum(binWeights.* repmat(PSDenv_C_NF, length(ModFreqs), 1),2);
        
        %     for modf_var = 1:length(ModFreqs)
        %         PowerModS(bs_var,modf_var) = Library.MODENERGY(ModFreqs(modf_var), PSDfreqVEC_Hz, PSDenv_A);
        %         PowerModN(bs_var,modf_var) = Library.MODENERGY(ModFreqs(modf_var), PSDfreqVEC_Hz, PSDenv_B);
        %         PowerModSN(bs_var,modf_var) = Library.MODENERGY(ModFreqs(modf_var), PSDfreqVEC_Hz, PSDenv_C);
        %         PowerModS_noisefloor(bs_var,modf_var) = Library.MODENERGY(ModFreqs(modf_var), PSDfreqVEC_Hz, PSDenv_A_NF);
        %         PowerModN_noisefloor(bs_var,modf_var) = Library.MODENERGY(ModFreqs(modf_var), PSDfreqVEC_Hz, PSDenv_B_NF);
        %         PowerModSN_noisefloor(bs_var,modf_var) = Library.MODENERGY(ModFreqs(modf_var), PSDfreqVEC_Hz, PSDenv_C_NF);
        %     end
    end
    
    %% Plot modPs (these figures and PowerMod_STRUCT are the only things that we should look at)
% % % %     if paramsIN.plt
% % % %         lw=2;
% % % %         col='bkr';
% % % % %         figure(figHandles.modPPlots); clf; hold on;
% % % %         
% % % %         errorbar(mean(PowerModS), std(PowerModS),[col(1) '-'],'linewidth',lw);
% % % %         errorbar(mean(PowerModSN), std(PowerModSN),[col(2) '-'],'linewidth',lw);
% % % %         errorbar(mean(PowerModN), std(PowerModN),[col(3) '-'],'linewidth',lw);
% % % %         errorbar(mean(PowerModS_noisefloor), std(PowerModSN_noisefloor),[col(1) ':'],'linewidth',lw);
% % % %         errorbar(mean(PowerModSN_noisefloor), std(PowerModSN_noisefloor),[col(2) ':'],'linewidth',lw);
% % % %         errorbar(mean(PowerModN_noisefloor), std(PowerModSN_noisefloor),[col(3) ':'],'linewidth',lw);
% % % %         
% % % %         title(sprintf('CF=%.3f kHz; SRtype=%d (1: LSR; 2: MSR; 3: HSR);\nCohc=%.2f; Cihc=%.2f; Overall (Signal) Level: %.1f dB SPL; SNR = %.1f dB', ...
% % % %             CF_C_kHz,paramsIN.FiberType,99,99,paramsIN.level,paramsIN.SNR2use_dB),'FontSize',12);axis tight; set(gca, 'xticklabels', ModFreqs);...
% % % %             xlabel('modFreq'); ylabel('modPower'); set (gcf, 'Units', 'normalized', 'Position', [.1,.1,.8,.8]);
% % % %         legend('S','SN','N','NF[S]','NF[SN]','NF[N]');
% % % %         Library.saveFigureAs([resultsDir 'eps' filesep 'modPnew' resultPostfix '.eps'])
% % % %         Library.saveFigureAs([resultsDir 'png' filesep 'modPnew' resultPostfix '.tiff'])
% % % %     end
    
    %% Pass back structures to make cleaner (to fit in noise floor comps)
% % %     PSDenv_STRUCT.PSDenv_S=PSDenv_S;
% % %     PSDenv_STRUCT.PSDenv_N=PSDenv_N;
% % %     PSDenv_STRUCT.PSDenv_SN=PSDenv_SN;
% % %     PSDenv_STRUCT.PSDfreqVEC_Hz=PSDfreqVEC_Hz;
% % %     PSDenv_STRUCT.PSDenv_S_noisefloor=PSDenv_S_noisefloor;
% % %     PSDenv_STRUCT.PSDenv_N_noisefloor=PSDenv_N_noisefloor;
% % %     PSDenv_STRUCT.PSDenv_SN_noisefloor=PSDenv_SN_noisefloor;
    
    % PSDtfs_STRUCT.PSDtfs_S=PSDtfs_S;
    % PSDtfs_STRUCT.PSDtfs_N=PSDtfs_N;
    % PSDtfs_STRUCT.PSDtfs_SN=PSDtfs_SN;
    % PSDtfs_STRUCT.PSDfreqVEC_Hz=PSDfreqVEC_Hz;
    % PSDtfs_STRUCT.PSDtfs_S_noisefloor=PSDtfs_S_noisefloor;
    % PSDtfs_STRUCT.PSDtfs_N_noisefloor=PSDtfs_N_noisefloor;
    % PSDtfs_STRUCT.PSDtfs_SN_noisefloor=PSDtfs_SN_noisefloor;
    
    % if CF_A_kHz < tfsCutOff && CF_B_kHz < tfsCutOff && CF_C_kHz < tfsCutOff
    %     NYQ=0.5*(1/paramsOUT.DELAYbinwidth_sec);
    %     PowerTfs_STRUCT.PowerTfsS = Library.MODENERGY(NYQ, PSDfreqVEC_Hz, PSDtfs_S);
    %     PowerTfs_STRUCT.PowerTfsN = Library.MODENERGY(NYQ, PSDfreqVEC_Hz, PSDtfs_N);
    %     PowerTfs_STRUCT.PowerTfsSN = Library.MODENERGY(NYQ, PSDfreqVEC_Hz, PSDtfs_SN);
    %     PowerTfs_STRUCT.PowerTfsS_noisefloor = Library.MODENERGY(NYQ, PSDfreqVEC_Hz, PSDtfs_S_noisefloor);
    %     PowerTfs_STRUCT.PowerTfsN_noisefloor = Library.MODENERGY(NYQ, PSDfreqVEC_Hz, PSDtfs_N_noisefloor);
    %     PowerTfs_STRUCT.PowerTfsSN_noisefloor = Library.MODENERGY(NYQ, PSDfreqVEC_Hz, PSDtfs_SN_noisefloor);
    % else
    %     PowerTfs_STRUCT.PowerTfsS = 0;
    %     PowerTfs_STRUCT.PowerTfsN = 0;
    %     PowerTfs_STRUCT.PowerTfsSN = 0;
    %     PowerTfs_STRUCT.PowerTfsS_noisefloor = 0;
    %     PowerTfs_STRUCT.PowerTfsN_noisefloor = 0;
    %     PowerTfs_STRUCT.PowerTfsSN_noisefloor = 0;
    % end
    
    PowerMod_STRUCT.PowerModS=PowerModS;
    PowerMod_STRUCT.PowerModN=PowerModN;
    PowerMod_STRUCT.PowerModSN=PowerModSN;
    PowerMod_STRUCT.PowerModS_noisefloor=PowerModS_noisefloor;
    PowerMod_STRUCT.PowerModN_noisefloor=PowerModN_noisefloor;
    PowerMod_STRUCT.PowerModSN_noisefloor=PowerModSN_noisefloor;
    
    params=repmat(struct([]),BootstrapLoopMax,1);
    for i=1:BootstrapLoopMax
        if ~isempty(SACSCCmetrics{i})
            params(i).nLines=SACSCCmetrics{i}.nLines;
            params(i).NumDrivenSpikes=SACSCCmetrics{i}.NumDrivenSpikes;
            params(i).AvgRate_sps=SACSCCmetrics{i}.AvgRate_sps;
        end
    end
    PowerMod_STRUCT.params=params;
else
    PSDenv_STRUCT={};
    PSDtfs_STRUCT={};
    PowerMod_STRUCT={};
    PowerTfs_STRUCT={};
end
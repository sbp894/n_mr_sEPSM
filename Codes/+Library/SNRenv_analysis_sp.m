%%
%for now, it is defined only for 1 SNR

function ModEP=SNRenv_analysis_sp(Power_Mat,paramsIN)

PowerMod_S=zeros(length(Power_Mat),size(Power_Mat(1).PowerMod_STRUCT.PowerMod_S,2));
PowerMod_N=zeros(length(Power_Mat),size(Power_Mat(1).PowerMod_STRUCT.PowerMod_N,2));
PowerMod_SN=zeros(length(Power_Mat),size(Power_Mat(1).PowerMod_STRUCT.PowerMod_SN,2));

PowerMod_S_noisefloor=zeros(length(Power_Mat),size(Power_Mat(1).PowerMod_STRUCT.PowerMod_S_noisefloor,2));
PowerMod_N_noisefloor=zeros(length(Power_Mat),size(Power_Mat(1).PowerMod_STRUCT.PowerMod_N_noisefloor,2));
PowerMod_SN_noisefloor=zeros(length(Power_Mat),size(Power_Mat(1).PowerMod_STRUCT.PowerMod_SN_noisefloor,2));

for i=1:length(Power_Mat)
    PowerMod_STRUCT=Power_Mat(i).PowerMod_STRUCT;
    PowerMod_S(i,:)=PowerMod_STRUCT.PowerMod_S(end,:);
    PowerMod_N(i,:)=PowerMod_STRUCT.PowerMod_N(end,:);
    PowerMod_SN(i,:)=PowerMod_STRUCT.PowerMod_SN(end,:);
    PowerMod_S_noisefloor(i,:)=PowerMod_STRUCT.PowerMod_S_noisefloor(end,:);
    PowerMod_N_noisefloor(i,:)=PowerMod_STRUCT.PowerMod_N_noisefloor(end,:);
    PowerMod_SN_noisefloor(i,:)=PowerMod_STRUCT.PowerMod_SN_noisefloor(end,:);
end


NoiseFloor_PercTolerance = 10;  % Apply to S and SN conditions.

Nan_INDs_S=((PowerMod_S - PowerMod_S_noisefloor)./PowerMod_S_noisefloor*100 < NoiseFloor_PercTolerance);
PowerMod_S(Nan_INDs_S)=NaN;
if sum(sum(Nan_INDs_S))
%     fprintf('*** %d values in PowerMod_S were set to NaN (ignored) because they were not above the Noise Floor!',sum(sum(Nan_INDs_S)));
end

% PowerMod_SN
% PowerMod_SN_noisefloor
% (PowerMod_SN - PowerMod_SN_noisefloor)./PowerMod_SN_noisefloor*100

Nan_INDs_SN=((PowerMod_SN - PowerMod_SN_noisefloor)./PowerMod_SN_noisefloor*100 < NoiseFloor_PercTolerance);
PowerMod_SN(Nan_INDs_SN)=NaN;
if sum(sum(Nan_INDs_SN))
%     fprintf('*** %d values in PowerMod_SN were set to NaN (ignored) because they were not above the Noise Floor!',sum(sum(Nan_INDs_SN)));
end

%% don't apply to noise (?does this make sense?)
% Calculate the SNR of the envelope of speech by subtracting the SNR of
% noise envelope from the SNR of speech and noise.

SNR_env_SN_N=zeros(size(PowerMod_S,1),size(PowerMod_S,2));
SNR_env_SN_N_noisefloor=zeros(size(PowerMod_S,1),size(PowerMod_S,2));

for j = 1:size(PowerMod_S,2)
    for i = 1:size(PowerMod_S,1)
        SNR_env_SN_N(i,j) = Library.EnvelopeSNR(PowerMod_SN(i,j), PowerMod_N(i,j));      
        SNR_env_SN_N_noisefloor(i,j) = Library.EnvelopeSNR(PowerMod_SN_noisefloor(i,j), PowerMod_N_noisefloor(i,j));      
    end
end


%% Plot ENVpower vs. MF
ModEP_S=nanmean(10*log10(PowerMod_S),1);
ModEP_N=nanmean(10*log10(PowerMod_N),1);
ModEP_SN=nanmean(10*log10(PowerMod_SN),1);

figure(102); 
set (gcf, 'Units', 'normalized', 'Position', [.1 .1 .8 .8]);

plot(ModEP_S,'ks-','MarkerSize',8,'MarkerFaceColor','g','LineWidth',1.5)
hold on
plot(ModEP_SN,'ko-','MarkerSize',8,'MarkerFaceColor','k','LineWidth',1.5)
plot(ModEP_N,'k^-','MarkerSize',8,'MarkerFaceColor','r','LineWidth',1.5)
hold off
% title(sprintf('\nCFs=[%s]kHz; SRtype=%d (1: LSR; 2: MSR; 3: HSR);\nOHCloss=%.2f dB; IHCloss=%.2f dB; Overall (Signal) Level: %.1f dB SPL; SNR=%.1f dB',CFstr(1:end-2),fiberType,OHCloss_dB,IHCloss_dB,OALevel_dBSPL,SNR2use_dB));

title(sprintf('\nCFs=[%1.2f]kHz; SRtype=%d (1: LSR; 2: MSR; 3: HSR);SNR=%.1f dB',paramsIN.CF_A_Hz/1e3,paramsIN.FiberType,paramsIN.SNR2use_dB));

ylabel('Envelope Power (dB)')
xlabel('Modulation Band CF (Hz)')
set(gca,'XTickLabel',{'1','2','4','8','16','32','64'})
legend('Clean Speech','Noisy Speech','Noise alone')

TotalSNRenv_SN_N_dB_noisefloor = 10*log10(sqrt(nansum(nansum(SNR_env_SN_N_noisefloor.^2)))); %% Total SNR of the envelope for speech across all channels and modulation frequencies

if sum(~isnan(SNR_env_SN_N))
    TotalSNRenv_SN_N_dB = 10*log10(sqrt(nansum(nansum(SNR_env_SN_N.^2)))); %% Total SNR of the envelope for speech across all channels and modulation frequencies
else 
    TotalSNRenv_SN_N_dB = TotalSNRenv_SN_N_dB_noisefloor;
end

if isinf(sum(sum(TotalSNRenv_SN_N_dB_noisefloor)))
   TotalSNRenv_SN_N_dB_noisefloor(isinf(TotalSNRenv_SN_N_dB_noisefloor))=NaN;
end

% CFstr=sprintf('%.2f, ',paramsIN.CF_A_Hz);
% fprintf('\nCFs=[%s]kHz; SRtype=%d (1: LSR; 2: MSR; 3: HSR);',CFstr(1:end-2),paramsIN.FiberType);
% fprintf('Acoustic SNR = %.2f dB\n',paramsIN.SNR2use_dB);
% fprintf('TotalSNRenv_SN_N_dB = %.2f dB  [SNR Noise Floor = %.2f dB]\n',TotalSNRenv_SN_N_dB,TotalSNRenv_SN_N_dB_noisefloor);

ModEP=struct('ModEP_S',ModEP_S,'ModEP_N',ModEP_N,'ModEP_SN',ModEP_SN,'TotalSNRenv_SN_N_dB',TotalSNRenv_SN_N_dB,'TotalSNRenv_SN_N_dB_noisefloor',TotalSNRenv_SN_N_dB_noisefloor);
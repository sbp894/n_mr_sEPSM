%%
%for now, it is defined only for 1 SNR

function ModEP=SNRenv_analysis_sp3(Power_Mat,paramsIN)

PowerMod_STRUCT=Power_Mat.PowerMod_STRUCT;
PowerMod_S=PowerMod_STRUCT.PowerMod_S;
PowerMod_N=PowerMod_STRUCT.PowerMod_N;
PowerMod_SN=PowerMod_STRUCT.PowerMod_SN;
PowerMod_S_noisefloor=PowerMod_STRUCT.PowerMod_S_noisefloor;
PowerMod_N_noisefloor=PowerMod_STRUCT.PowerMod_N_noisefloor;
PowerMod_SN_noisefloor=PowerMod_STRUCT.PowerMod_SN_noisefloor;

PowerMod_S=(PowerMod_S-PowerMod_S_noisefloor)./PowerMod_S_noisefloor;
PowerMod_N=(PowerMod_N-PowerMod_N_noisefloor)./PowerMod_N_noisefloor;
PowerMod_SN=(PowerMod_SN-PowerMod_SN_noisefloor)./PowerMod_SN_noisefloor;

NoiseFloor_PercTolerance = 10;  % Apply to S and SN conditions.
Noisefloor_value=.1;

Nan_INDs_S=PowerMod_S<NoiseFloor_PercTolerance/100;
PowerMod_S(Nan_INDs_S)=Noisefloor_value;
if sum(sum(Nan_INDs_S))
    %     fprintf('*** %d values in PowerMod_S were set to NaN (ignored) because they were not above the Noise Floor!',sum(sum(Nan_INDs_S)));
end

% PowerMod_SN
% PowerMod_SN_noisefloor
% (PowerMod_SN - PowerMod_SN_noisefloor)./PowerMod_SN_noisefloor*100

Nan_INDs_SN=((PowerMod_SN - PowerMod_SN_noisefloor)./PowerMod_SN_noisefloor*100 < NoiseFloor_PercTolerance);
PowerMod_SN(Nan_INDs_SN)=Noisefloor_value;
if sum(sum(Nan_INDs_SN))
    %     fprintf('*** %d values in PowerMod_SN were set to NaN (ignored) because they were not above the Noise Floor!',sum(sum(Nan_INDs_SN)));
end


Nan_INDs_N=((PowerMod_N - PowerMod_N_noisefloor)./PowerMod_N_noisefloor*100 < NoiseFloor_PercTolerance);
PowerMod_N(Nan_INDs_N)=Noisefloor_value;
if sum(sum(Nan_INDs_N))
    %     fprintf('*** %d values in PowerMod_N were set to NaN (ignored) because they were not above the Noise Floor!',sum(sum(Nan_INDs_SN)));
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
% ModEP_S=10*log10(PowerMod_S);
% ModEP_N=10*log10(PowerMod_N);
% ModEP_SN=10*log10(PowerMod_SN);

ModEP_S=10*log10(PowerMod_S);
ModEP_N=10*log10(PowerMod_N);
ModEP_SN=10*log10(PowerMod_SN);

nanmean(ModEP_N)

figure(102);
% set (gcf, 'Units', 'normalized', 'Position', [.1 .1 .8 .8]);

errorbar(nanmean(ModEP_S), nanstd(ModEP_S),'ks-','MarkerSize',8,'MarkerFaceColor','g','LineWidth',1.5)
hold on
errorbar(nanmean(ModEP_SN), nanstd(ModEP_SN),'ko-','MarkerSize',8,'MarkerFaceColor','k','LineWidth',1.5)
errorbar(nanmean(ModEP_N), nanstd(ModEP_N),'k^-','MarkerSize',8,'MarkerFaceColor','r','LineWidth',1.5)
hold off
% title(sprintf('\nCFs=[%s]kHz; SRtype=%d (1: LSR; 2: MSR; 3: HSR);\nOHCloss=%.2f dB; IHCloss=%.2f dB; Overall (Signal) Level: %.1f dB SPL; SNR=%.1f dB',CFstr(1:end-2),fiberType,OHCloss_dB,IHCloss_dB,OALevel_dBSPL,SNR2use_dB));

title(sprintf('\nCFs=[%1.2f]kHz; SRtype=%d (1: LSR; 2: MSR; 3: HSR);SNR=%.1f dB',paramsIN.CF_A_Hz/1e3,paramsIN.FiberType,paramsIN.SNR2use_dB));

axis tight;
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

ModEP=struct('ModEP_S',ModEP_S,'ModEP_N',ModEP_N,'ModEP_SN',ModEP_SN, ...
    'TotalSNRenv_SN_N_dB',TotalSNRenv_SN_N_dB,'TotalSNRenv_SN_N_dB_noisefloor',TotalSNRenv_SN_N_dB_noisefloor);
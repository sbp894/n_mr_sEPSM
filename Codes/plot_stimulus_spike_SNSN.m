function plot_stimulus_spike_SNSN(SpikeTrains,plot_assist,resultsDir,resultPostfix)

global figHandles;

CurStims=plot_assist.CurStims;
ANmodel_Fs_Hz=plot_assist.Fs;
paramsIN=plot_assist.paramsIN;


Tbin=0.01;

stim_A_model=CurStims{1,1}{1};
stim_B_model=CurStims{2,1}{1};
stim_C_model=CurStims{3,1}{1};

spiketimes_Ap=get_type_polar_spikes(SpikeTrains{1,1});
spiketimes_An=get_type_polar_spikes(SpikeTrains{1,2});
spiketimes_Bp=get_type_polar_spikes(SpikeTrains{2,1});
spiketimes_Bn=get_type_polar_spikes(SpikeTrains{2,2});
spiketimes_Cp=get_type_polar_spikes(SpikeTrains{3,1});
spiketimes_Cn=get_type_polar_spikes(SpikeTrains{3,2});

figure(figHandles.meanRates);
clf;
set (gcf, 'Units', 'normalized', 'Position', [.1 .1 .8 .8]);
set(gcf,'Visible','Off');
ymaxTEMP=-Inf;
Tmax=max([length(stim_A_model)/ANmodel_Fs_Hz,length(stim_B_model)/ANmodel_Fs_Hz,length(stim_C_model)/ANmodel_Fs_Hz,...
    max(spiketimes_Ap(:,1)),max(spiketimes_An(:,1)),max(spiketimes_Bp(:,1)),max(spiketimes_Bn(:,1)),max(spiketimes_Cp(:,1)),...
    max(spiketimes_Cn(:,1))]);

timeVEC1_sec=0:1/ANmodel_Fs_Hz:1.1*Tmax;

subplot(411)
plot(timeVEC1_sec,[stim_C_model; zeros(length(timeVEC1_sec)-length(stim_C_model),1)],'b');
hold on; 
plot(timeVEC1_sec,[stim_B_model; zeros(length(timeVEC1_sec)-length(stim_B_model),1)],'r'); 
plot(timeVEC1_sec,[stim_A_model; zeros(length(timeVEC1_sec)-length(stim_A_model),1)],'g');
legend('SN','N','S');

title(sprintf('CF=%.3f kHz; SRtype=%d (1: LSR; 2: MSR; 3: HSR); SNR = %.1f dB', paramsIN.CF_C_Hz,paramsIN.FiberType,paramsIN.SNR2use_dB));
ylabel(sprintf('Stimulus Amplitude\n(Pascals)'))

subplot(412)
% timeVEC1_sec=((1:length(spiketimes_Ap))-1)/ANmodel_Fs_Hz;

% plot(spiketimes_Ap(:,1),spiketimes_Ap(:,2),'k');
histogram(spiketimes_Ap,ceil(Tmax/Tbin));
% hold on; 
% plot(spiketimes_An(:,1),spiketimes_An(:,2),'g');
histogram(spiketimes_An,ceil(Tmax/Tbin));
ylabel(sprintf('Synapse piketimesput\n(spikes/sec)'))
title('Signal');
% legend('POS','NEG');
ymaxTEMP=max([max(ylim) ymaxTEMP]);

subplot(413)
% plot(spiketimes_Bp(:,1),spiketimes_Bp(:,2),'k');
histogram(spiketimes_Bp,ceil(Tmax/Tbin));
% hold on; 
% plot(spiketimes_Bn(:,1),spiketimes_Bn(:,2),'r');
% histogram(spiketimes_Bn,ceil(Tmax/Tbin));
ylabel(sprintf('Synapse piketimesput\n(spikes/sec)'))
title('Noise')
% legend('POS','NEG')
ymaxTEMP=max([max(ylim) ymaxTEMP]);

subplot(414)
% plot(spiketimes_Cp(:,1),spiketimes_Cp(:,2),'k');
histogram(spiketimes_Cp,ceil(Tmax/Tbin));
% hold on; 
% plot(spiketimes_Cn(:,1),spiketimes_Cn(:,2),'b');
% histogram(spiketimes_Cn,ceil(Tmax/Tbin));
ylabel(sprintf('Synapse piketimesput\n(spikes/sec)'))
title('Signal + Noise')
xlabel('Time (sec)')
% legend('POS','NEG')
ymaxTEMP=max([max(ylim) ymaxTEMP]);

subplot(412); ylim([0 ymaxTEMP]); xlim([0 Tmax]);
subplot(413); ylim([0 ymaxTEMP]); xlim([0 Tmax]);
subplot(414); ylim([0 ymaxTEMP]); xlim([0 Tmax]);

% Library.saveFigureAs([resultsDir 'SpeechSpikePlots_png' filesep 'SSPlots' resultPostfix '.png']);
Library.saveFigureAs([resultsDir 'SpeechSpikePlots_eps' filesep 'SSPlots' resultPostfix '.eps']);
saveas(gcf,[resultsDir 'SpeechSpikePlots_bmp' filesep 'SSPlots' resultPostfix '.tiff']);
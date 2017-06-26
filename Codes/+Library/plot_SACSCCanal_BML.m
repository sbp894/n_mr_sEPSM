function plot_SACSCCanal_BML(SACSCCfunctions,SACSCCmetrics,paramsIN)
% File: plot_CCCanal(SACSCCfunctions,SACSCCmetrics,paramsIN)
%
% plot_CCCanal_0: used with CCCanal_0 - plots ANY individual rep
% plot_CCCanal_1: used with CCCanal_1 - plots 1 rep, but using RS.avgPSD/CSD
% plot_CCCanal_2: used with CCCanal_2 - plots rep 1, but also shows AVG/STD vals 
% plot_CCCanal_3: for CCCanal_3 - plots all BOOTSTRAP reps, along with all AVGs 
%
% M. Heinz May 23, 2008
% Plot 15-panel plot for SAC/SCC_CCC analysis of 4 spike trains.  Assumes
% all SAC/SCC/CCC analyiss has been done already by CCCanal.m, and the
% functions and metrics are just passed.


% Find indices and Freq Ranges for SCpeak and CCCenvs TOUSE
DCpeak_TOUSEindex=find(strcmp(paramsIN.DCpeak_TOUSE,SACSCCmetrics{end}.DCpeaks_legend));
% CCCtfs_TOUSEindex=find(strcmp(paramsIN.CCCtfs_TOUSE,SACSCCmetrics{end}.CCCtfs_legend));
SCpeak_TOUSEindex=find(strcmp(paramsIN.SCpeak_TOUSE,SACSCCmetrics{end}.SCpeaks_legend));
% CCCenv_TOUSEindex=find(strcmp(paramsIN.CCCenv_TOUSE,SACSCCmetrics{end}.CCCenvs_legend));
DC0peak_TOUSEindex=find(strcmp('DCpeak_0',SACSCCmetrics{end}.DCpeaks_legend));
SC0peak_TOUSEindex=find(strcmp('IFFTraw_0',SACSCCmetrics{end}.SCpeaks_legend));
% if ceil(CCCenv_TOUSEindex/2)>length(SACSCCmetrics{end}.sums.PSD_LHfreqs_Hz)
% 	doSUMS=0;  % These are "add-ons" at the end of CCCenvs (e.g., "raw", "adjSC"), where sums from LOWmodHz to HIGHmod_Hz are undefined
% else
% 	doSUMS=1;
% end
% if doSUMS
% 	CCCenv_LOWmod_Hz=SACSCCmetrics{end}.sums.PSD_LHfreqs_Hz(ceil(CCCenv_TOUSEindex/2),1);
% 	CCCenv_HIGHmod_Hz=SACSCCmetrics{end}.sums.PSD_LHfreqs_Hz(ceil(CCCenv_TOUSEindex/2),2);
% else
% 	CCCenv_LOWmod_Hz=NaN;  
% 	CCCenv_HIGHmod_Hz=NaN;
% end

% Params for all reps and AVG
REPcolors={'b','r','g','c','y','m','b','r','g','c','y','m',};
AVGcolor='k';
Nreps=length(SACSCCfunctions)-1;
AVGind=length(SACSCCfunctions);

% User-specified plot LIMITS - can be specified from outside by including
% in paramsIN
if isfield(paramsIN,'XLIMIT_delay'), XLIMIT_delay=paramsIN.XLIMIT_delay;, else XLIMIT_delay=paramsIN.MAXdelay_sec*1e3;   end
if isfield(paramsIN,'XLIMIT_PSDhigh'), XLIMIT_PSDhigh=paramsIN.XLIMIT_PSDhigh;, else XLIMIT_PSDhigh=64;   end
if isfield(paramsIN,'YLIMIT_SClow'), YLIMIT_SClow=paramsIN.YLIMIT_SClow;, else YLIMIT_SClow=0;   end
% Set X/Y LIMITS to be consistent across panels
YLIMIT_SAC=-999;
for repIND=1:Nreps+1
	YLIMIT_SAC=max([YLIMIT_SAC max(SACSCCfunctions{repIND}.SAC_A_avg) max(SACSCCfunctions{repIND}.XpAC_A_avg)]);
end
YLIMIT_DC=-999;
for repIND=1:Nreps+1
	YLIMIT_DC=max([YLIMIT_DC max(abs(SACSCCfunctions{repIND}.DIFCOR_A))]);
end
YLIMIT_SC=-999;
for repIND=1:Nreps+1
	YLIMIT_SC=max([YLIMIT_SC max(SACSCCfunctions{repIND}.SUMCORadj_A)]);
end
YLIMIT_PSD=-999;
for repIND=1:Nreps+1
    YLIMIT_PSD=max([YLIMIT_PSD max(SACSCCfunctions{repIND}.PSDenv_A)]);
end


figure; clf
%% Shuffled Correlograms
subplot(4,3,1)
for repIND=1:Nreps
	plot(SACSCCfunctions{repIND}.delays_usec/1000,SACSCCfunctions{repIND}.SAC_A_avg,'Color',REPcolors{repIND},'LineWidth', 2.5); hold on
	plot(SACSCCfunctions{repIND}.delays_usec/1000,SACSCCfunctions{repIND}.XpAC_A_avg,'Color',REPcolors{repIND},'LineWidth', 1);
end
plot(SACSCCfunctions{AVGind}.delays_usec/1000,SACSCCfunctions{AVGind}.SAC_A_avg,'Color',AVGcolor,'LineWidth', 2.5);
plot(SACSCCfunctions{AVGind}.delays_usec/1000,SACSCCfunctions{AVGind}.XpAC_A_avg,'Color',AVGcolor,'LineWidth', 1);
plot(XLIMIT_delay*[-1 1],ones(1,2),'-k');
plot(zeros(1,2),YLIMIT_SAC*[0 1],'--k');
hold off
set(gca, 'Box', 'off', 'TickDir', 'out');
xlim(XLIMIT_delay*[-1 1]); ylim(YLIMIT_SAC*[0 1])
title(sprintf('CONDITION A [S]\n'),'FontSize',14);
text(0.02,0.95,sprintf('%d/%d',round(SACSCCmetrics{AVGind}.NumDrivenSpikes(1,1)),round(SACSCCmetrics{AVGind}.NumDrivenSpikes(1,2))), ...
	'units','norm','FontSize',8)
text(.5,1,sprintf('CF = %.3f kHz',paramsIN.CF_A_Hz/1000),'units','norm','VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',8)

%%DIFCORs
subplot(4,3,4);
DCpeak_A_TEXT='[';
for repIND=1:Nreps
	plot(SACSCCfunctions{repIND}.delays_usec/1000,SACSCCfunctions{repIND}.DIFCOR_A,'Color',REPcolors{repIND}); hold on
	DCpeak_A_TEXT=sprintf('%s%.2f, ',DCpeak_A_TEXT,SACSCCmetrics{repIND}.DCpeak_A(DC0peak_TOUSEindex));   % ACF always at ZERO
end
plot(SACSCCfunctions{AVGind}.delays_usec/1000,SACSCCfunctions{AVGind}.DIFCOR_A,'Color',AVGcolor,'Linewidth',2);
plot(zeros(1,2),YLIMIT_DC*[-1 1],'--k');
plot(XLIMIT_delay*[-1 1],zeros(1,2),'-k');hold off
xlim(XLIMIT_delay*[-1 1]); ylim(YLIMIT_DC*[-1 1]);
title(sprintf('DIFCOR_A (peak=%.2f %s] ["%s"] )',SACSCCmetrics{AVGind}.DCpeak_A(DC0peak_TOUSEindex),DCpeak_A_TEXT(1:end-2),paramsIN.DCpeak_TOUSE), ...
	'Interpreter','none','FontSize',8);
set(gca, 'Box', 'off', 'TickDir', 'out');
ylabel(sprintf('NORMALIZED\n# COINCIDENCES'),'FontSize',12)

%%SUMCORs
subplot(4,3,7);
SCpeak_A_TEXT='[';
for repIND=1:Nreps
	plot(SACSCCfunctions{repIND}.delays_usec/1000,SACSCCfunctions{repIND}.SUMCORadj_A,'Color',REPcolors{repIND}); hold on
	SCpeak_A_TEXT=sprintf('%s%.2f, ',SCpeak_A_TEXT,SACSCCmetrics{repIND}.SCpeaks_A(SC0peak_TOUSEindex));  % ACF always at ZERO
end
plot(SACSCCfunctions{AVGind}.delays_usec/1000,SACSCCfunctions{AVGind}.SUMCORadj_A,'Color',AVGcolor,'Linewidth',2);
plot(XLIMIT_delay*[-1 1],ones(1,2),'-k');
plot(zeros(1,2),YLIMIT_SC*[0 1],'--k'); hold off
xlim(XLIMIT_delay*[-1 1]); ylim([YLIMIT_SClow YLIMIT_SC]);
title(sprintf('SUMCOR_A (peak=%.2f %s] ["%s"] )',SACSCCmetrics{AVGind}.SCpeaks_A(SC0peak_TOUSEindex),SCpeak_A_TEXT(1:end-2),paramsIN.SCpeak_TOUSE), ...
	'Interpreter','none','FontSize',8);
set(gca, 'Box', 'off', 'TickDir', 'out');

%% PSDs
subplot(4,3,10);
for repIND=1:Nreps
	plot(SACSCCfunctions{repIND}.PSD_freqVEC,SACSCCfunctions{repIND}.PSDenv_A,'Color',REPcolors{repIND});  hold on
	plot(SACSCCfunctions{repIND}.PSD_freqVEC,SACSCCfunctions{repIND}.rand.PSDenv_A,'Color',REPcolors{repIND},'LineStyle',':');  
end
plot(SACSCCfunctions{AVGind}.PSD_freqVEC,SACSCCfunctions{AVGind}.PSDenv_A,'Color',AVGcolor,'Linewidth',2);  xlim([0 XLIMIT_PSDhigh]); ylim([0 YLIMIT_PSD]);
plot(SACSCCfunctions{AVGind}.PSD_freqVEC,SACSCCfunctions{AVGind}.rand.PSDenv_A,'Color',AVGcolor,'Linewidth',2,'LineStyle',':');  
% plot(ones(1,2)*CCCenv_LOWmod_Hz,[0 YLIMIT_PSD],'--k','linewidth',2);
% plot(ones(1,2)*CCCenv_HIGHmod_Hz,[0 YLIMIT_PSD],'--k','linewidth',2); hold off
xlim([0 XLIMIT_PSDhigh]); ylim([0 YLIMIT_PSD]); 
set(gca, 'Box', 'off', 'TickDir', 'out');
ylabel(sprintf('SPECTRAL DENSITY\nAMPLITUDE'),'FontSize',12,'HorizontalAlignment','center')

orient landscape

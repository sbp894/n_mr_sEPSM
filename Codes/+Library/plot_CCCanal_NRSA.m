function plot_CCCanal_NRSA(SACSCCfunctions,SACSCCmetrics,paramsIN)

% Find indices and Freq Ranges for SCpeak and CCCenvs TOUSE
DCpeak_TOUSEindex=find(strcmp(paramsIN.DCpeak_TOUSE,SACSCCmetrics{end}.DCpeaks_legend));
% CCCtfs_TOUSEindex=find(strcmp(paramsIN.CCCtfs_TOUSE,SACSCCmetrics{end}.CCCtfs_legend));
SCpeak_TOUSEindex=find(strcmp(paramsIN.SCpeak_TOUSE,SACSCCmetrics{end}.SCpeaks_legend));
CCCenv_TOUSEindex=find(strcmp(paramsIN.CCCenv_TOUSE,SACSCCmetrics{end}.CCCenvs_legend));
DC0peak_TOUSEindex=find(strcmp('DCpeak_0',SACSCCmetrics{end}.DCpeaks_legend));
SC0peak_TOUSEindex=find(strcmp('IFFTraw_0',SACSCCmetrics{end}.SCpeaks_legend));

REPcolors={'b','r','g','c','y','m',};
AVGcolor='k';
Nreps=length(SACSCCfunctions)-1;
AVGind=length(SACSCCfunctions);

% User-specified plot LIMITS - can be specified from outside by including
% in paramsIN
if isfield(paramsIN,'XLIMIT_delay') 
    XLIMIT_delay=paramsIN.XLIMIT_delay; 
else
    XLIMIT_delay = 20e3; % 50ms
end
if isfield(paramsIN,'YLIMIT_SClow')
    YLIMIT_SClow=paramsIN.YLIMIT_SClow;
else
    YLIMIT_SClow=0.8;   
end


ticDivider = 1e3; %show tics in ms
% Set X/Y LIMITS to be consistent across panels
YLIMIT_SAC=-999;
for repIND=1:Nreps+1
	YLIMIT_SAC=max([YLIMIT_SAC max(SACSCCfunctions{repIND}.SAC_A_avg) max(SACSCCfunctions{repIND}.XpAC_A_avg) max(SACSCCfunctions{repIND}.SAC_C_avg) ...
		max(SACSCCfunctions{repIND}.XpAC_C_avg) max(SACSCCfunctions{repIND}.SCC_AC_avg) max(SACSCCfunctions{repIND}.XpCC_AC_avg)]);
end
YLIMIT_DC=-999;
for repIND=1:Nreps+1
	YLIMIT_DC=max([YLIMIT_DC max(abs(SACSCCfunctions{repIND}.DIFCOR_A)) max(abs(SACSCCfunctions{repIND}.DIFCOR_C)) max(abs(SACSCCfunctions{repIND}.DIFCOR_AC))]);
end
YLIMIT_SC=-999;
for repIND=1:Nreps+1
	YLIMIT_SC=max([YLIMIT_SC max(SACSCCfunctions{repIND}.SUMCOR_A) max(SACSCCfunctions{repIND}.SUMCOR_C) max(SACSCCfunctions{repIND}.SUMCOR_AC)]);
end

TICKlength=0.02;
figure_prop_name = {'PaperPositionMode','units','Position'};
figure_prop_val =  { 'auto'            ,'inches', [3.75    0.25    6.75    5.2]};
figure(3); %clf
set(3,figure_prop_name,figure_prop_val,'Visible','Off');

% XLIMIT_delay=6;

%% Shuffled Correlograms
h1=subplot(4,3,1);
% plot(SACSCCfunctions{AVGind}.delays_usec/1000,SACSCCfunctions{AVGind}.SAC_A_avg,'Color',AVGcolor,'LineWidth', 2.5); hold on
plot(SACSCCfunctions{AVGind}.delays_usec,SACSCCfunctions{AVGind}.SAC_A_avg,'Color',AVGcolor,'LineWidth', 2.5); hold on
plot(SACSCCfunctions{AVGind}.delays_usec,SACSCCfunctions{AVGind}.XpAC_A_avg,'Color',AVGcolor,'LineWidth', 1);
plot(XLIMIT_delay*[-1 1],ones(1,2),'-k');
plot(zeros(1,2),YLIMIT_SAC*[0 1],'--k');
hold off
set(gca, 'Box', 'off', 'TickDir', 'out');
xlim(XLIMIT_delay*[-1 1]); ylim(YLIMIT_SAC*[0 1])
title(sprintf('CONDITION A (S)\n'),'FontSize',14);
% text(0.9,1.1,'SAC and XpAC','units','norm','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',12)
% xlim([-10 10]);
text(0.02,0.85,sprintf('CF=%d Hz',paramsIN.CF_A_Hz), ...
'units','norm','FontSize',12)
tickRowVec = [-XLIMIT_delay:2*XLIMIT_delay/4:XLIMIT_delay];
set(h1,'TickLength',[TICKlength 0.025],'XTick',tickRowVec,'XTickLabel',tickRowVec./ticDivider)

h2=subplot(4,3,2);
plot(SACSCCfunctions{AVGind}.delays_usec,SACSCCfunctions{AVGind}.SAC_C_avg,'Color',AVGcolor,'LineWidth', 2.5);hold on
plot(SACSCCfunctions{AVGind}.delays_usec,SACSCCfunctions{AVGind}.XpAC_C_avg,'Color',AVGcolor,'LineWidth', 1);
plot(XLIMIT_delay*[-1 1],ones(1,2),'-k');
plot(zeros(1,2),YLIMIT_SAC*[0 1],'--k');
hold off
xlim(XLIMIT_delay*[-1 1]); ylim(YLIMIT_SAC*[0 1])
set(gca, 'Box', 'off', 'TickDir', 'out');
title(sprintf(['CONDITION C (SN) \n SNR ' num2str(paramsIN.SNR2use_dB) 'dB']),'FontSize',14);
% xlim([-10 10]);
set(h2,'TickLength',[TICKlength 0.025],'XTick',tickRowVec,'XTickLabel',tickRowVec./ticDivider)

h3=subplot(4,3,3);
plot(SACSCCfunctions{AVGind}.delays_usec,SACSCCfunctions{AVGind}.SCC_AC_avg,'Color',AVGcolor,'LineWidth', 2.5);hold on
plot(SACSCCfunctions{AVGind}.delays_usec,SACSCCfunctions{AVGind}.XpCC_AC_avg,'Color',AVGcolor,'LineWidth', 1);
plot(XLIMIT_delay*[-1 1],ones(1,2),'-k');
plot(zeros(1,2),YLIMIT_SAC*[0 1],'--k');
xlim(XLIMIT_delay*[-1 1]); ylim(YLIMIT_SAC*[0 1]); hold off
set(gca, 'Box', 'off', 'TickDir', 'out');
title(sprintf('CROSS-CORR (S,SN)\n'),'FontSize',14);
% text(0.5,1.1,'SCC and XpCC','units','norm','VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',12)
% xlim([-10 10]);
set(h3,'TickLength',[TICKlength 0.025],'XTick',tickRowVec,'XTickLabel',tickRowVec./ticDivider)

%%DIFCORs
h4=subplot(4,3,4);
plot(SACSCCfunctions{AVGind}.delays_usec,SACSCCfunctions{AVGind}.DIFCOR_A,'Color',AVGcolor,'Linewidth',2);hold on
plot(zeros(1,2),YLIMIT_DC*[-1 1],'--k');
plot(XLIMIT_delay*[-1 1],zeros(1,2),'-k');hold off
xlim(XLIMIT_delay*[-1 1]); ylim(YLIMIT_DC*[-1 1]);
title(sprintf('DIFCOR A'), ...
	'Interpreter','none','FontSize',8);
set(gca, 'Box', 'off', 'TickDir', 'out');
ylabel(sprintf('NORMALIZED\n# COINCIDENCES'),'FontSize',10)
% xlim([-10 10]);
set(h4,'TickLength',[TICKlength 0.025],'XTick',tickRowVec,'XTickLabel',tickRowVec./ticDivider)

h5=subplot(4,3,5);
plot(SACSCCfunctions{AVGind}.delays_usec,SACSCCfunctions{AVGind}.DIFCOR_C,'Color',AVGcolor,'Linewidth',2);hold on
plot(zeros(1,2),YLIMIT_DC*[-1 1],'--k');
plot(XLIMIT_delay*[-1 1],zeros(1,2),'-k');hold off
xlim(XLIMIT_delay*[-1 1]); ylim(YLIMIT_DC*[-1 1]);
title(sprintf('DIFCOR C','Interpreter','none','FontSize',8));
set(gca, 'Box', 'off', 'TickDir', 'out');
%xlim([-10000 10000]);
set(h5,'TickLength',[TICKlength 0.025],'XTick',tickRowVec,'XTickLabel',tickRowVec./ticDivider)

h6=subplot(4,3,6);
plot(SACSCCfunctions{AVGind}.delays_usec,SACSCCfunctions{AVGind}.DIFCOR_AC,'Color',AVGcolor,'Linewidth',2);hold on
plot(zeros(1,2),YLIMIT_DC*[-1 1],'--k');
plot(XLIMIT_delay*[-1 1],zeros(1,2),'-k');hold off
xlim(XLIMIT_delay*[-1 1]); ylim(YLIMIT_DC*[-1 1]);
title(sprintf('DIFCOR AC'), ...
	'Interpreter','none','FontSize',8);
set(gca, 'Box', 'off', 'TickDir', 'out');
text(0.02,0.85,sprintf('rho_{TFS}=%.2f',SACSCCmetrics{end}.CCCtfs(strcmp('DCpeak_0',SACSCCmetrics{end}.DCpeaks_legend))), ...
'units','norm','FontSize',8)
%xlim([-10000 10000]);
set(h6,'TickLength',[TICKlength 0.025],'XTick',tickRowVec,'XTickLabel',tickRowVec./ticDivider)

%%SUMCORs
% XLIMIT_delay=100;
h7=subplot(4,3,7);
plot(SACSCCfunctions{AVGind}.delays_usec,SACSCCfunctions{AVGind}.SUMCORadj_A,'Color',AVGcolor,'Linewidth',2);hold on
plot(XLIMIT_delay*[-1 1],ones(1,2),'-k');
plot(zeros(1,2),YLIMIT_SC*[0 1],'--k'); hold off
xlim(XLIMIT_delay*[-1 1]); ylim([YLIMIT_SClow YLIMIT_SC]);
set(gca, 'Box', 'off', 'TickDir', 'out');
title(sprintf('SUMCOR A', ...
	'Interpreter','none','FontSize',8));
%xlim([-10000 10000]);
set(h7,'TickLength',[TICKlength 0.025],'XTick',tickRowVec,'XTickLabel',tickRowVec./ticDivider)

h8=subplot(4,3,8);
plot(SACSCCfunctions{AVGind}.delays_usec,SACSCCfunctions{AVGind}.SUMCORadj_C,'Color',AVGcolor,'Linewidth',2);hold on
plot(XLIMIT_delay*[-1 1],ones(1,2),'-k');
plot(zeros(1,2),YLIMIT_SC*[0 1],'--k'); hold off
xlim(XLIMIT_delay*[-1 1]); ylim([YLIMIT_SClow YLIMIT_SC]);
title(sprintf('SUMCOR C','Interpreter','none','FontSize',8));
set(gca, 'Box', 'off', 'TickDir', 'out');
%xlim([-10000 10000]);
xlabel(sprintf('DELAY (msec)'),'FontSize',10)
set(h8,'TickLength',[TICKlength 0.025],'XTick',tickRowVec,'XTickLabel',tickRowVec./ticDivider)

h9=subplot(4,3,9);
plot(SACSCCfunctions{AVGind}.delays_usec,SACSCCfunctions{AVGind}.SUMCORadj_AC,'Color',AVGcolor,'Linewidth',2);hold on
plot(XLIMIT_delay*[-1 1],ones(1,2),'-k');
plot(zeros(1,2),YLIMIT_SC*[0 1],'--k');hold off
xlim(XLIMIT_delay*[-1 1]);  ylim([YLIMIT_SClow YLIMIT_SC]);
set(gca, 'Box', 'off', 'TickDir', 'out');
title(sprintf('SUMCOR AC','Interpreter','none','FontSize',8));
%xlim([-10000 10000]);
text(0.02,0.85,sprintf('rho_{ENV}=%.2f',SACSCCmetrics{end}.CCCenvs(strcmp('IFFTrawSC_0',SACSCCmetrics{end}.CCCenvs_legend))), ...
'units','norm','FontSize',8)
set(h9,'TickLength',[TICKlength 0.025],'XTick',tickRowVec,'XTickLabel',tickRowVec./ticDivider)

h10=subplot(4,3,10);
plot(SACSCCfunctions{AVGind}.delays_usec,SACSCCfunctions{AVGind}.SUMCORadj_A_64,'Color',AVGcolor,'Linewidth',2);hold on
plot(XLIMIT_delay*[-1 1],ones(1,2),'-k');
plot(zeros(1,2),YLIMIT_SC*[0 1],'--k'); hold off
xlim(XLIMIT_delay*[-1 1]); ylim([YLIMIT_SClow YLIMIT_SC]);
set(gca, 'Box', 'off', 'TickDir', 'out');
modFRange = ' 0-64Hz';
title(sprintf(['SUMCOR AC adj' modFRange],'Interpreter','none','FontSize',8));
%xlim([-100 100]);
set(h10,'TickLength',[TICKlength 0.025],'XTick',tickRowVec,'XTickLabel',tickRowVec./ticDivider)

h11=subplot(4,3,11);
plot(SACSCCfunctions{AVGind}.delays_usec,SACSCCfunctions{AVGind}.SUMCORadj_C_64,'Color',AVGcolor,'Linewidth',2);hold on
plot(XLIMIT_delay*[-1 1],ones(1,2),'-k');
plot(zeros(1,2),YLIMIT_SC*[0 1],'--k'); hold off
xlim(XLIMIT_delay*[-1 1]); ylim([YLIMIT_SClow YLIMIT_SC]);
set(gca, 'Box', 'off', 'TickDir', 'out');
title(sprintf(['SUMCOR AC adj' modFRange],'Interpreter','none','FontSize',8));
%xlim([-100 100]);
set(h11,'TickLength',[TICKlength 0.025],'XTick',tickRowVec,'XTickLabel',tickRowVec./ticDivider)
xlabel(sprintf('DELAY (msec)'),'FontSize',10)

h12=subplot(4,3,12);
plot(SACSCCfunctions{AVGind}.delays_usec,SACSCCfunctions{AVGind}.SUMCORadj_AC_64,'Color',AVGcolor,'Linewidth',2);hold on
plot(XLIMIT_delay*[-1 1],ones(1,2),'-k');
plot(zeros(1,2),YLIMIT_SC*[0 1],'--k');hold off
xlim(XLIMIT_delay*[-1 1]);  ylim([YLIMIT_SClow YLIMIT_SC]);
set(gca, 'Box', 'off', 'TickDir', 'out');
title(sprintf(['SUMCOR AC adj' modFRange],'Interpreter','none','FontSize',8));
%xlim([-100 100]);
text(0.02,0.85,sprintf('rho_{ENV0-64}=%.2f',SACSCCmetrics{end}.CCCenvs(strcmp('0-64, subBIAS',SACSCCmetrics{end}.CCCenvs_legend))), ...
'units','norm','FontSize',8)
set(h12,'TickLength',[TICKlength 0.025],'XTick',tickRowVec,'XTickLabel',tickRowVec./ticDivider)

orient landscape

% %%
% clear;
% close all;
% clc;

%%
function plot_summary(ChinID)
% DataDir='Z:\Users\SP\Codes\SNRenv-OUTPUT\DataAnal\SP-2016_10_28-Q285_AN_Normal-1\';
% DataDir='Z:\Users\SP\Codes\SNRenv-OUTPUT\DataAnal\SP-2016_11_04-Q284_AN_normal-1\';

DataRepository='Z:\Users\SP\Codes\SNRenv-OUTPUT\DataAnal\';

temp=dir([DataRepository '*' num2str(ChinID) '*']);
DataDir=strcat(DataRepository,temp(end).name,filesep);

SummaryVar=load([DataDir 'SummaryVar.mat']);
SummaryVar=SummaryVar.SummaryVar;

CF_dB=unique(SummaryVar(:,[1 3]),'rows');
legstr=cell(size(CF_dB,1),1);
h1=nan(size(CF_dB,1),1);
for i=1:size(CF_dB,1)
    ind= prod(SummaryVar(:,[1 3])== repmat(CF_dB(i,:),size(SummaryVar,1),1),2)==1;
    SNR=SummaryVar(ind,2);
    SNRenv=SummaryVar(ind,5);
    SR=unique(SummaryVar(ind,4));
    
    figure(1);
    h1(i)=plot(SNR,SNRenv,'LINEWIDTH',1.2);
    hold on;
    plot(SNR,SNRenv,'*','color',get(h1(i),'color'));
    
    legstr{i}=sprintf('CF=%3.2d | dB= %d | SR= %2.1d',CF_dB(i,1), CF_dB(i,2), SR);
    
end
xlabel('SNR');
ylabel('SNRenv');
title(sprintf('Q-%d',ChinID));
legend(h1,legstr,'location','northwestoutside');
set (gcf, 'Units', 'normalized', 'Position', [.1 .1 .8 .8]);
saveas(gcf,[DataDir 'SummaryPlot.bmp']);

%%
figure(2);
% clf;
hold on;


SNRcolor='rgb';
MarkerType='ovsp';
plot(nan,nan,SNRcolor(3)); plot(nan,nan,SNRcolor(2)); plot(nan,nan,SNRcolor(1));

unique_CFs=unique(SummaryVar(:,1));
flag_dont_plot=0;

for cf_var=1:length(unique_CFs)
    indCF=SummaryVar(:,1)==unique_CFs(cf_var);
    unique_dBs=sort(unique(SummaryVar(indCF,3)));
    
    for dB_var=1:length(unique_dBs)
        db_cur=unique_dBs(dB_var);
        ind2check=find(and(indCF, SummaryVar(:,3)==db_cur)==1);
        markerType=MarkerType(dB_var);
        
        for snr_var=1:length(ind2check)
            indcur=ind2check(snr_var);
            snrcur=SummaryVar(indcur,2);
            
            
            switch snrcur
                case -6,
                    Color2Use=SNRcolor(1);
                case 0,
                    Color2Use=SNRcolor(2);
                case 6,
                    Color2Use=SNRcolor(3);
                otherwise, 
                    flag_dont_plot=1;
            end
            if ~flag_dont_plot
                plot(SummaryVar(indcur,4), SummaryVar(indcur,5), 'color', Color2Use, 'marker', markerType);
            else 
                flag_dont_plot=0;
            end
        end
                
    end
                
end
set (gcf, 'Units', 'normalized', 'Position', [.1 .1 .8 .8]);
xlabel('SR');
ylabel('SNRenv');
legend('+6','0','-6');
title(sprintf('Q-%d',ChinID));
saveas(gcf,[DataDir 'SR_SNRenv.bmp']);


%%
figure(3);
% clf;
hold on;


SNRcolor='rgb';
MarkerType='ovsp';
plot(nan,nan,SNRcolor(3)); plot(nan,nan,SNRcolor(2)); plot(nan,nan,SNRcolor(1));

unique_CFs=unique(SummaryVar(:,1));
flag_dont_plot=0;

for cf_var=1:length(unique_CFs)
    indCF=SummaryVar(:,1)==unique_CFs(cf_var);
    unique_dBs=sort(unique(SummaryVar(indCF,3)));
    
    for dB_var=1:length(unique_dBs)
        db_cur=unique_dBs(dB_var);
        ind2check=find(and(indCF, SummaryVar(:,3)==db_cur)==1);
        markerType=MarkerType(dB_var);
        
        for snr_var=1:length(ind2check)
            indcur=ind2check(snr_var);
            snrcur=SummaryVar(indcur,2);
            
            
            switch snrcur
                case -6,
                    Color2Use=SNRcolor(1);
                case 0,
                    Color2Use=SNRcolor(2);
                case 6,
                    Color2Use=SNRcolor(3);
                otherwise, 
                    flag_dont_plot=1;
            end
            if ~flag_dont_plot
                plot(SummaryVar(indcur,1), SummaryVar(indcur,5), 'color', Color2Use, 'marker', markerType);
            else 
                flag_dont_plot=0;
            end
        end
                
    end
                
end
set (gcf, 'Units', 'normalized', 'Position', [.1 .1 .8 .8]);
xlabel('CF (Hz)');
ylabel('SNRenv');
legend('+6','0','-6');
title(sprintf('Q-%d',ChinID));
saveas(gcf,[DataDir 'CF_SNRenv.bmp']);
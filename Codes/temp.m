clear;
close all;
clc;

resultsDir='D:\Study Stuff\Matlab\SNRenv-SimData-Updated4PSDavg\Output\Simulation\20170606-1\';
resultTxt='CF_%1.4fk_Sent_%i_%s_level_%1.2f_SNR%i_fType%1.0f';
lw=1.6;

%%
if resultsDir(end)~=filesep
    resultsDir=[resultsDir filesep];
end
SimDatafName=[resultsDir 'SimData.mat'];

%%
condition=load([resultsDir 'conditions.mat']);
csCell=condition.csCell;

cfs=unique(cell2mat(csCell(:,1)));
sents=3; %unique(cell2mat(csCell(:,2)));
levels=unique(cell2mat(csCell(:,3)));
noisetypes=unique(csCell(:,4));
snrs=unique(cell2mat(csCell(:,5)));
fibertypes=unique(cell2mat(csCell(:,6)));


% ModEP_SNR_CF=zeros(length(cfs),length(levels),length(snrs))-11;
% ModEP_NF_SNR_CF=zeros(length(cfs),length(levels),length(snrs))-11;% Noisefloor

for ntypeVar=1:length(noisetypes)
    for cfVar=1:3:length(cfs)
        figure;
        for levelVar=1:length(levels)
            for snrVar=1:2:length(snrs)
                subplot(3,2,(snrVar+1)/2);
                for sentVar=1:length(sents)
                    
                    for ftypeVar=1:length(fibertypes)
                        resultPostfix=sprintf(resultTxt,  cfs(cfVar)/1e3, sents(sentVar), noisetypes{ntypeVar}, levels(levelVar), snrs(snrVar), fibertypes(ftypeVar));
                        PSD_STRUCT=load([resultsDir 'psd' filesep 'psd' resultPostfix '.mat']);
                        plot(PSD_STRUCT.PSDfreqVEC_Hz, PSD_STRUCT.PSDenv_SN, 'linewidth',lw);
                        hold on;
                        plot(PSD_STRUCT.PSDfreqVEC_Hz, PSD_STRUCT.PSDenv_N, 'g');
                        
                    end
                    hold off;
                    title(sprintf('%s, snr=%i, cf=%1.2f', ... %, snr-brick=%1.2f, snr-from_finalpsd=%1.2f',...
                        noisetypes{ntypeVar}, snrs(snrVar), cfs(cfVar))); %, snrBrick, snrFinal));
                    xlim([0 20]);
                    ylabel('PSD'); xlabel('freq');
                    fprintf('%i/%i: %i/%i: %i/%i -- %s\n',ntypeVar, length(noisetypes), snrVar, length(snrs), cfVar, length(cfs), resultPostfix);
                    
                end
            end
        end
        legend('SN-LSR', 'N-LSR','SN-MSR', 'N-MSR','SN-HSR', 'N-HSR');
    end
end
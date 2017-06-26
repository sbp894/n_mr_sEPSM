%%
% Loads saved PSD files, reads the PowerMod_STRUCT structure that contains
% the boostrap estimats of modulation power.
 
 
 
%%
function parse_saved_data_for_SNRenv3(resultsDir,resultTxt)
 
% clear; close all; clc;
% resultsDir='D:\Study Stuff\Matlab\SNRenv-SimData-Updated4PSDavg\Output\Simulation\20170609-2\';
% resultTxt='CF_%1.4fk_Sent_%i_%s_level_%1.2f_SNR%i_fType%1.0f';
 
NumSTDsNFtol=3;
CompareWithNF=1;
 
%%
if resultsDir(end)~=filesep
    resultsDir=[resultsDir filesep];
end
SimDatafName=[resultsDir 'SimData.mat'];
 
%%
condition=load([resultsDir 'conditions.mat']);
ExpControlParams=load([resultsDir 'ExpControlParams.mat']);
ExpControlParams=ExpControlParams.ExpControlParams;
csCell=condition.csCell;
 
cfs=unique(cell2mat(csCell(:,1)));
sents=3; unique(cell2mat(csCell(:,2)));
levels=unique(cell2mat(csCell(:,3)));
noisetypes=unique(csCell(:,4));
snrs=unique(cell2mat(csCell(:,5)));
fibertypes=unique(cell2mat(csCell(:,6)));
 
% ModEP_SNR_CF=zeros(length(cfs),length(levels),length(snrs))-11;
% ModEP_NF_SNR_CF=zeros(length(cfs),length(levels),length(snrs))-11;% Noisefloor
 
 
%% Descrition of loops
% structure: first noise types (SSN or SAM)
% for each noise type, look at each level
% For each level, find SNRenv
% for each SNR, one SNRenv value based on averages across sentences
% For each sentence, do a weighted average across spont rates (fibertypes)
% For each fibertypes, sum across all CFs
 
%%
figOutDir=[fileparts(resultsDir(1:end-1)) filesep 'sEPSM_final_outputs' filesep];
if ~exist(figOutDir)
    mkdir(figOutDir);
end
figOutName=sprintf('CompareWithNF_%i__fixSPL_%i__CorrWindow_0add1mul_%i__sent_num_%i', CompareWithNF, ExpControlParams.fixSPL, ExpControlParams.winCorr0Add1Mul, sents);

%%
SRweights=[1 1 1]/3;
ModFreqs=[1 2 4 8 16 32 64];
clear veryFirstFlag;
% if ~exist(SimDatafName, 'file')
if 1
    SimData=repmat(struct('SNRenv',nan), length(noisetypes), length(levels), length(snrs), length(sents), length(cfs), length(fibertypes));
    
    for ntypeVar=1:length(noisetypes)
        for levelVar=1:length(levels)
            for snrVar=1:length(snrs)
                for sentVar=1:length(sents)
                    for cfVar=1:length(cfs)
                        for ftypeVar=1:length(fibertypes)
                            
                            resultPostfix=sprintf(resultTxt,  cfs(cfVar)/1e3, sents(sentVar), noisetypes{ntypeVar}, levels(levelVar), snrs(snrVar), fibertypes(ftypeVar));
                            PSD_STRUCT=load([resultsDir 'psd' filesep 'psd' resultPostfix '.mat']);
                            
                            PowerMod_SN=PSD_STRUCT.PowerMod_SN;
                            pModSN=mean(PowerMod_SN);
                            if CompareWithNF
                                PowerMod_SN_noisefloor=PSD_STRUCT.PowerMod_SN_noisefloor;
                                zSN=(mean(PowerMod_SN)-mean(PowerMod_SN_noisefloor))./(sqrt(std(PowerMod_SN).*std(PowerMod_SN_noisefloor)));
                                nanInds=zSN<NumSTDsNFtol;
                                pModSN(nanInds)=nan;
                            end
                            
                            PowerMod_N=PSD_STRUCT.PowerMod_N;
                            pModN=mean(PowerMod_N);
                            if CompareWithNF
                                PowerMod_N_noisefloor=PSD_STRUCT.PowerMod_N_noisefloor;
                                u_nNF=mean(PowerMod_N_noisefloor);
                                std_nNF=std(PowerMod_N_noisefloor);
                                zN=(mean(PowerMod_N)-u_nNF)./sqrt(std(PowerMod_N).*std_nNF);
                                nanInds=zN<NumSTDsNFtol;
                                %                                 pModN(nanInds)=max([pModN(nanInds); u_nNF(nanInds)])+NumSTDsNFtol*std_nNF(nanInds);
                                pModN(nanInds)=u_nNF(nanInds)+NumSTDsNFtol*std_nNF(nanInds);
                            end
                            
                            SNRenv=(pModSN-pModN)./pModN;
                            SNRenv(SNRenv<0)=nan;
                            
                            SimData(ntypeVar, levelVar, snrVar, sentVar, cfVar, ftypeVar).PowerModSN = pModSN;
                            SimData(ntypeVar, levelVar, snrVar, sentVar, cfVar, ftypeVar).PowerModN = pModN;
                            SimData(ntypeVar, levelVar, snrVar, sentVar, cfVar, ftypeVar).SNRenv= SNRenv;
                        end
                        
                        fprintf('%i/%i: %i/%i: %i/%i -- %s\n',ntypeVar, length(noisetypes), snrVar, length(snrs), cfVar, length(cfs), resultPostfix);
                    end
                end
            end
        end
    end
    
    save(SimDatafName,'SimData');
    
else
    load(SimDatafName);
end
 
 
%% Weighted across SRs
%  Assumed: one level and one sentence
OutputWeighted=repmat(struct('pModSN',[],'pModN',[],'SNRenv',[] ), length(noisetypes), length(snrs));
for ntypeVar=1:length(noisetypes)
    figure(ntypeVar);
    set(gcf,'units', 'normalized','position', [.1 .1 .8 .8]);
    subplotCount=1;
    for snrVar=1:length(snrs)
        PowerModSN=zeros(1, length(cfs)*length(ModFreqs));
        PowerModN=zeros(1, length(cfs)*length(ModFreqs));
        SNRenv=zeros(1, length(cfs)*length(ModFreqs));
        
        for ftypeVar=1:length(fibertypes)
            tempPowerModSN=[SimData(ntypeVar, 1, snrVar, 1, :, ftypeVar).PowerModSN]*SRweights(fibertypes(ftypeVar));
            tempPowerModSN(isnan(tempPowerModSN))=0;
            PowerModSN=PowerModSN+tempPowerModSN;
            
            tempPowerModSN=[SimData(ntypeVar, 1, snrVar, 1, :, ftypeVar).PowerModN]*SRweights(fibertypes(ftypeVar));
            tempPowerModSN(isnan(tempPowerModSN))=0;
            PowerModSN=PowerModSN+tempPowerModSN;
            
            tempSNRenv=[SimData(ntypeVar, 1, snrVar, 1, :, ftypeVar).SNRenv]*SRweights(fibertypes(ftypeVar));
            tempSNRenv(isnan(tempSNRenv))=0;
            SNRenv=SNRenv+tempSNRenv;
        end
        
        
        OutputWeighted(ntypeVar, snrVar).pModSN=reshape(PowerModSN, length([SimData(1, 1, 1, 1, 1, 1).PowerModSN]), length(cfs))';
        OutputWeighted(ntypeVar, snrVar).pModN=reshape(PowerModN, length([SimData(1, 1, 1, 1, 1).PowerModN]), length(cfs))';
        OutputWeighted(ntypeVar, snrVar).SNRenv=reshape(SNRenv, length([SimData(1, 1, 1, 1, 1).SNRenv]), length(cfs))';
        
        
        OutputWeighted(ntypeVar, snrVar).SNRenvSingle=nansum(nansum(OutputWeighted(ntypeVar, snrVar).SNRenv.^2));
        
        subplot(length(snrs),1,subplotCount);
        imagesc(OutputWeighted(ntypeVar, snrVar).SNRenv)
        set(gca, 'Ytick', 1:size(OutputWeighted(ntypeVar, snrVar).SNRenv,1), 'yticklabel',  strread(num2str(uint16(cfs')),'%s'),'Xtick', 1:length(ModFreqs), 'xticklabel',  strread(num2str(ModFreqs),'%s'))
        title(sprintf('%s, SNR=%i, SNR_e_n_v=%1.2f', noisetypes{ntypeVar}, snrs(snrVar), OutputWeighted(ntypeVar, snrVar).SNRenvSingle));
        colorbar;
        subplotCount=subplotCount+1;
        %         end
    end
    saveas(gcf, [figOutDir figOutName '__' noisetypes{ntypeVar} '.bmp'], 'bmp');
    saveas(gcf, [figOutDir figOutName '__' noisetypes{ntypeVar}]);
end
 
%%
markers='dos<>pv';
figure(length(noisetypes)+1);
set(gcf,'units', 'normalized','position', [.1 .1 .8 .8]);
lw=2;
legstr=cell(length(noisetypes),1);
hold on;
for ntypeVar=1:length(noisetypes)
    plot(snrs, 20*log10([OutputWeighted(ntypeVar, :).SNRenvSingle]), [ '-' markers(ntypeVar)],'linewidth',lw);
    %     plot(snrs, 20*log10([OutputWeighted(ntypeVar, :).SNRenvSingleNF]), [ 'g--' markers(ntypeVar)],'linewidth',lw);
    legstr{ntypeVar}=noisetypes{ntypeVar};
    %     legstr{2*ntypeVar}=[noisetypes{ntypeVar} '-NF'];
end
 
set(gca, 'XTick', snrs)
legend(legstr, 'location', 'northwest');
xlabel('Acoustic SNR (dB)');
ylabel('SNR_e_n_v (dB)');
title(strrep(figOutName, '_', '-'));
grid on;

saveas(gcf, [figOutDir figOutName '__SNRenv_vs_SNR.bmp'], 'bmp');
saveas(gcf, [figOutDir figOutName '__SNRenv_vs_SNR']);
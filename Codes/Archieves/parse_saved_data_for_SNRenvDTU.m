% close all;
% clc;

% function parse_saved_data_for_SNRenvDTU(resultsDir,resultTxt)
resultsDir='D:\Study Stuff\Matlab\SNRenv-SimData-Updated4PSDavg\Output\Simulation\20170605-1_imp\';
resultTxt='CF_%1.4fk_Sent_%i_%s_level_%1.2f_SNR%i_fType%1.0f';


%%
if resultsDir(end)~=filesep
    resultsDir=[resultsDir filesep];
end
SimDatafName=[resultsDir 'SimData.mat'];

%%
condition=load([resultsDir 'conditions.mat']);
csCell=condition.csCell;

cfs=unique(cell2mat(csCell(:,1)));
sents=1; %unique(cell2mat(csCell(:,2)));
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
SRweights=[.2 .2 .6];
ModFreqs=[1 2 4 8 16 32 64];
clear veryFirstFlag;
if ~exist(SimDatafName, 'file')
    SimData=repmat(struct('SNRenv',nan), length(noisetypes), length(levels), length(snrs), length(sents), length(cfs));
    
    for ntypeVar=1:length(noisetypes)
        for levelVar=1:length(levels)
            for snrVar=1:length(snrs)
                for sentVar=1:length(sents)
                    for cfVar=1:length(cfs)
                        
                        if exist('veryFirstFlag', 'var')
                            PSDenv_SN=0*PSDenv_SN;
                            PSDenv_SN_noisefloor=0*PSDenv_SN_noisefloor;
                            PSDenv_N=0*PSDenv_N;
                            PSDenv_N_noisefloor=0*PSDenv_N_noisefloor;
                        end
                        
                        for ftypeVar=1:length(fibertypes)
                            
                            resultPostfix=sprintf(resultTxt,  cfs(cfVar)/1e3, sents(sentVar), noisetypes{ntypeVar}, levels(levelVar), snrs(snrVar), fibertypes(ftypeVar));
                            PSD_STRUCT=load([resultsDir 'psd' filesep 'psd' resultPostfix '.mat'], 'PSDenv_SN','PSDenv_SN_noisefloor','PSDenv_N','PSDenv_N_noisefloor','PSDfreqVEC_Hz');
                            
                            if ~exist('veryFirstFlag', 'var')
                                veryFirstFlag=1;
                                
                                PSDenv_SN=SRweights(ftypeVar)*PSD_STRUCT.PSDenv_SN;
                                PSDenv_SN_noisefloor=SRweights(ftypeVar)*PSD_STRUCT.PSDenv_SN_noisefloor;
                                PSDenv_N=SRweights(ftypeVar)*PSD_STRUCT.PSDenv_N;
                                PSDenv_N_noisefloor=SRweights(ftypeVar)*PSD_STRUCT.PSDenv_N_noisefloor;
                            else
                                PSDenv_SN=PSDenv_SN+SRweights(ftypeVar)*PSD_STRUCT.PSDenv_SN;
                                PSDenv_SN_noisefloor=PSDenv_SN_noisefloor+SRweights(ftypeVar)*PSD_STRUCT.PSDenv_SN_noisefloor;
                                PSDenv_N=PSDenv_N+SRweights(ftypeVar)*PSD_STRUCT.PSDenv_N;
                                PSDenv_N_noisefloor=PSDenv_N_noisefloor+SRweights(ftypeVar)*PSD_STRUCT.PSDenv_N_noisefloor;
                            end
                        end
                        
                        plt =0;
                        qFactor=ones(size(ModFreqs));
                        
                        binWeights = modFbank(PSD_STRUCT.PSDfreqVEC_Hz,ModFreqs, qFactor, plt);
                        
                        pModSN= nansum(binWeights.* repmat(PSDenv_SN, length(ModFreqs), 1),2);
                        pModSN_NF=nansum(binWeights.* repmat(PSDenv_SN_noisefloor, length(ModFreqs), 1),2);
                        pModN=nansum(binWeights.* repmat(PSDenv_N, length(ModFreqs), 1),2);
                        pModN_NF=nansum(binWeights.* repmat(PSDenv_N_noisefloor, length(ModFreqs), 1),2);
                        
                        pModSNnorm=(pModSN-pModSN_NF)./pModSN_NF;
                        pModNnorm=(pModN-pModN_NF)./pModN_NF;
                        
                        SNRenv=(pModSN-pModN)./pModN;
                        SNRenv(SNRenv<0)=0;
                        SNRenvnorm=(pModSNnorm-pModNnorm)./pModNnorm;
                        SNRenvnorm(SNRenvnorm<0)=0;
                        SNRenvNF=(pModSN_NF-pModN_NF)./pModN_NF;
                        SNRenvNF(SNRenvNF<0)=0;
                        
                        SimData(ntypeVar, levelVar, snrVar, sentVar, cfVar).PowerModSN = pModSN;
                        SimData(ntypeVar, levelVar, snrVar, sentVar, cfVar).PowerModSN_NF = pModSN_NF;
                        SimData(ntypeVar, levelVar, snrVar, sentVar, cfVar).PowerModN = pModN;
                        SimData(ntypeVar, levelVar, snrVar, sentVar, cfVar).PowerModN_NF = pModN_NF;
                        
                        SimData(ntypeVar, levelVar, snrVar, sentVar, cfVar).PowerModSNnorm= pModSNnorm;
                        SimData(ntypeVar, levelVar, snrVar, sentVar, cfVar).PowerModNnorm= pModNnorm;
                        
                        SimData(ntypeVar, levelVar, snrVar, sentVar, cfVar).SNRenv= SNRenv;
                        SimData(ntypeVar, levelVar, snrVar, sentVar, cfVar).SNRenvnorm= SNRenvnorm;
                        SimData(ntypeVar, levelVar, snrVar, sentVar, cfVar).SNRenvNF= SNRenvNF;
                        
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
OutputWeighted=repmat(struct('pModSNnorm',[],'pModNnorm',[],'SNRenvnorm',[] ), length(noisetypes), length(snrs));
for ntypeVar=1:length(noisetypes)
    figure(ntypeVar);
    subplotCount=1;
    for snrVar=1:length(snrs)
%         OutputWeighted(ntypeVar, snrVar).pModSNnorm=reshape([SimData(ntypeVar, 1, snrVar, 1, :).PowerModSNnorm], length([SimData(1, 1, 1, 1, 1).PowerModSNnorm]), length(cfs))';
%         OutputWeighted(ntypeVar, snrVar).pModNnorm=reshape([SimData(ntypeVar, 1, snrVar, 1, :).PowerModNnorm], length([SimData(1, 1, 1, 1, 1).PowerModNnorm]), length(cfs))';
%         OutputWeighted(ntypeVar, snrVar).SNRenvnorm=reshape([SimData(ntypeVar, 1, snrVar, 1, :).SNRenvnorm], length([SimData(1, 1, 1, 1, 1).SNRenvnorm]), length(cfs))';
        OutputWeighted(ntypeVar, snrVar).SNRenv=reshape([SimData(ntypeVar, 1, snrVar, 1, :).SNRenv], length([SimData(1, 1, 1, 1, 1).SNRenv]), length(cfs))';
%         OutputWeighted(ntypeVar, snrVar).SNRenvNF=reshape([SimData(ntypeVar, 1, snrVar, 1, :).SNRenvNF], length([SimData(1, 1, 1, 1, 1).SNRenvNF]), length(cfs))';
        
%         OutputWeighted(ntypeVar, snrVar).SNRenvSinglenorm=nansum(nansum(OutputWeighted(ntypeVar, snrVar).SNRenvnorm.^2));
        OutputWeighted(ntypeVar, snrVar).SNRenvSingle=nansum(nansum(OutputWeighted(ntypeVar, snrVar).SNRenv.^2));
%         OutputWeighted(ntypeVar, snrVar).SNRenvSingleNF=nansum(nansum(OutputWeighted(ntypeVar, snrVar).SNRenvNF.^2));
        
        %         if rem(snrVar,2)
        subplot(3,4,subplotCount);
        imagesc(OutputWeighted(ntypeVar, snrVar).SNRenv)
        set(gca, 'Ytick', 1:size(OutputWeighted(ntypeVar, snrVar).SNRenv,1), 'yticklabel',  strread(num2str(uint16(cfs')),'%s'),'Xtick', 1:length(ModFreqs), 'xticklabel',  strread(num2str(ModFreqs),'%s'))
        title(sprintf('%s, SNR=%i, SNR_e_n_v=%1.2f', noisetypes{ntypeVar}, snrs(snrVar), OutputWeighted(ntypeVar, snrVar).SNRenvSingle));
        colorbar;
        subplotCount=subplotCount+1;
        %         end
    end
end

%%
markers='dos<>pv';
figure;
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
title('Averaged across CFs and SRs');
grid on;
% function parse_saved_data_for_SNRenv2(resultsDir,resultTxt)
resultsDir='D:\Study Stuff\Matlab\SNRenv-SimData-Updated4PSDavg\Output\Simulation\20170530-2\';
resultTxt='CF_%1.4fk_Sent_%i_%s_level_%1.2f_SNR%i';

%%
if resultsDir(end)~=filesep
    resultsDir=[resultsDir filesep];
end

condition=load([resultsDir 'conditions.mat']);
cs=condition.csCell;
list_of_cfs=cell2mat(cs(:,1));
unique_cfs=unique(list_of_cfs);
unique_snr=unique(cell2mat(cs(:,3)));
unique_levels=unique(cell2mat(cs(:,5)));


ModEP_SNR_CF=zeros(length(unique_cfs),length(unique_levels),length(unique_snr))-11;
ModEP_NF_SNR_CF=zeros(length(unique_cfs),length(unique_levels),length(unique_snr))-11;% Noisefloor

MinSNRs=[1, 3, 5]; % based on [-6 -3 0 3 6] to [-6 0 6] conversion

for cf_var=1:length(unique_cfs)
    cur_cf=unique_cfs(cf_var);
    ind1=find(list_of_cfs==cur_cf);
    db_for_cf=unique(cell2mat(cs(ind1,5)));
    
    for j=1:length(db_for_cf)
        ind=ind1(cell2mat(cs(ind1,5))==db_for_cf(j));
        db_var=find(unique_levels==db_for_cf(j));
        if rem(length(ind),length(unique_snr))==0 && length(ind)~=length(unique_snr) %% to separate out multiple dBs, take care of the case when nSNR is not equal to uniqueSNR
           ind=ind(end/2+1:end);
        elseif length(ind)==6 % Hardcoded as of now, 6 means 2 different SPLs 
            
        end
        
        for snr_var=1:length(ind)
            
            cs_var=ind(snr_var);
            
            
            switch sign(cs{cs_var,3})
                
                case -1
                    resultPostfixN=sprintf([resultTxt(1:end-5) 'xx' resultTxt(end-4:end)],       cs{cs_var,1}/1e3,  cs{cs_var,2},    cs{cs_var,4},  cs{cs_var,5}, cs{cs_var,3});
                case 0
                    resultPostfixN=sprintf([resultTxt(1:end-5) 'yy' resultTxt(end-4:end)],       cs{cs_var,1}/1e3,  cs{cs_var,2},    cs{cs_var,4},  cs{cs_var,5}, cs{cs_var,3});
                case 1
                    resultPostfixN=sprintf([resultTxt(1:end-5) 'zz' resultTxt(end-4:end)],       cs{cs_var,1}/1e3,  cs{cs_var,2},    cs{cs_var,4},  cs{cs_var,5}, cs{cs_var,3});
            end
            resultPostfix=sprintf(resultTxt,       cs{cs_var,1}/1e3,  cs{cs_var,2},  cs{cs_var,4}, cs{cs_var,5},    cs{cs_var,3});
            
            Power_Mat.PowerMod_STRUCT=load([resultsDir 'psd' filesep 'psd' resultPostfix '.mat'],...
                'PowerMod_S','PowerMod_N','PowerMod_SN','PowerMod_S_noisefloor','PowerMod_N_noisefloor','PowerMod_SN_noisefloor');
            Power_Mat.CF=cs{cs_var,1};
            
            disp(['----------------------------------\n' resultPostfix]);
            
            load([resultsDir 'paramsIN' filesep 'paramsIN' resultPostfix '.mat']);
            
            ModEP=Library.SNRenv_analysis_sp(Power_Mat,paramsIN);
            title(strrep(resultPostfixN,'_','|'));
            if length(unique_snr)~=length(ind)
                SNRind=MinSNRs(snr_var);
            else 
                SNRind=(snr_var);   
            end
            
            ModEP_SNR_CF(cf_var,db_var,SNRind)=ModEP.TotalSNRenv_SN_N_dB;
            ModEP_NF_SNR_CF(cf_var,db_var,SNRind)=ModEP.TotalSNRenv_SN_N_dB_noisefloor;
            
            Library.saveFigureAs([resultsDir 'envPowereps' filesep 'envPower' resultPostfixN '.eps']);
            Library.saveFigureAs([resultsDir 'envPowerpng' filesep 'envPower' resultPostfixN '.png']);
            Library.parsave([resultsDir 'ModEP' filesep resultPostfixN '.mat'],ModEP);
            
            clear paramsIN;
        end
    end
end


%%
figure;
clf;
hold on;
ModEP_SNR_CF(isinf(ModEP_SNR_CF))=NaN;
ModEP_NF_SNR_CF(isinf(ModEP_NF_SNR_CF))=NaN;

ModEP_NF_SNR_CF((ModEP_NF_SNR_CF-ModEP_SNR_CF)>0)=ModEP_SNR_CF((ModEP_NF_SNR_CF-ModEP_SNR_CF)>0);
count=0;
legend_str=cell(1);
TotalModEP=zeros(size(unique_snr));

for cf_var=1:size(ModEP_SNR_CF,1)
    for db_var=1:size(ModEP_SNR_CF,2)
        if ModEP_SNR_CF(cf_var,db_var,1)~=-11
            count=count+1;
            plot(unique_snr,squeeze(ModEP_SNR_CF(cf_var,db_var,:))-squeeze(ModEP_NF_SNR_CF(cf_var,db_var,:)),'-d','LINEWIDTH',2);
            TotalModEP=TotalModEP+squeeze(ModEP_SNR_CF(cf_var,db_var,:))-squeeze(ModEP_NF_SNR_CF(cf_var,db_var,:));
            legend_str{count}=sprintf('CF=%1.2fHz @%1.2f dB SPL',unique_cfs(cf_var),unique_levels(db_var));
        end
    end
end

plot(unique_snr,TotalModEP,'-o','LINEWIDTH',3);
% leg_in=strsplit(legend_str);
% leg_in=leg_in(~strcmp(leg_in,''));

legend([legend_str {' Sum'}],'location','eastoutside');
xlabel('SNR (dB)');
ylabel('ModEP (dB)')
title(' TOTAL SNRenv computed vs SNR');
Library.saveFigureAs([resultsDir 'envPowereps' filesep 'TotalSNRenvPLOT.eps']);
Library.saveFigureAs([resultsDir 'envPowerpng' filesep 'TotalSNRenvPLOT.png']);
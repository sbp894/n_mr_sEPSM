function parse_saved_data_for_SNRenv(resultsDir)

condition=load([resultsDir 'conditions.mat']);
cs=condition.csCell;
list_of_snr=cell2mat(cs(:,3));
unique_snr=unique(list_of_snr);
unique_cfs=unique(cell2mat(cs(:,1)));


ModEP_SNR_CF=zeros(length(unique_snr),length(unique_cfs));
ModEP_NF_SNR_CF=zeros(length(unique_snr),length(unique_cfs)); % Noisefloor

for snr_var=1:length(unique_snr)
    cur_snr=unique_snr(snr_var);
    ind=find((list_of_snr==cur_snr));
    
    %     Power_Mat(length(ind))=struct('PowerMod_STRUCT',[]); %#ok<*AGROW>
    %     resultPostfix=cell(length(ind),1);
    
    for cf_var=1:length(ind)
        cs_var=ind(cf_var);
        %         resultPostfix{cf_var}=sprintf('CF_%1.2fk_Sent_%i_SNR%i_Noise_%s',       cs{cs_var,1}/1e3,  cs{cs_var,2},    cs{cs_var,3},  cs{cs_var,4});
        %         Power_Mat(cf_var).PowerMod_STRUCT=load([resultsDir 'psd' filesep 'psd' resultPostfix '.mat'],'PowerMod_S','PowerMod_N','PowerMod_SN','PowerMod_S_noisefloor','PowerMod_N_noisefloor','PowerMod_SN_noisefloor');
        %         Power_Mat(cf_var).CF=cs{ind(cf_var),1};
        %
        resultPostfix=sprintf('CF_%1.2fk_Sent_%i_SNR%i_Noise_%s',       cs{cs_var,1}/1e3,  cs{cs_var,2},    cs{cs_var,3},  cs{cs_var,4});
        
        Power_Mat.PowerMod_STRUCT=load([resultsDir 'psd' filesep 'psd' resultPostfix '.mat'],'PowerMod_S','PowerMod_N','PowerMod_SN','PowerMod_S_noisefloor','PowerMod_N_noisefloor','PowerMod_SN_noisefloor');
        Power_Mat.CF=cs{cs_var,1};
        
        load([resultsDir 'paramsIN' filesep 'paramsIN' resultPostfix '.mat']);
        
        ModEP=Library.SNRenv_analysis_sp(Power_Mat,paramsIN);
        ModEP_SNR_CF(snr_var,cf_var)=ModEP.TotalSNRenv_SN_N_dB;
        ModEP_NF_SNR_CF(snr_var,cf_var)=ModEP.TotalSNRenv_SN_N_dB_noisefloor;
        
        Library.saveFigureAs([resultsDir 'envPowereps' filesep 'envPower' resultPostfix '.eps']);
        Library.saveFigureAs([resultsDir 'envPowerpng' filesep 'envPower' resultPostfix '.png']);
        Library.parsave([resultsDir 'ModEP' filesep resultPostfix '.mat'],ModEP);
        
        clear paramsIN;
    end
    
end

legend_str='';
figure;
hold on;

for cf_var=1:size(ModEP_SNR_CF,2)
    plot(unique_snr,ModEP_SNR_CF(:,cf_var)-ModEP_NF_SNR_CF(:,cf_var),'-d','LINEWIDTH',2);
    legend_str=sprintf('%s CF=%1.2fHz ',legend_str,unique_cfs(cf_var));
%     strcat(legend_str,['CF=' num2str(unique_cfs(cf_var)),' ']);
end
plot(unique_snr,nansum(ModEP_SNR_CF,2)-nansum(ModEP_NF_SNR_CF,2),'-o','LINEWIDTH',3);
leg_in=strsplit(legend_str);
leg_in=leg_in(~strcmp(leg_in,''));

legend([leg_in ' Sum']);
xlabel('SNR (dB)');
ylabel('ModEP (dB)')
title(' TOTAL SNRenv computed vs SNR');
Library.saveFigureAs([resultsDir 'envPowereps' filesep 'TotalSNRenvPLOT.eps']);
Library.saveFigureAs([resultsDir 'envPowerpng' filesep 'TotalSNRenvPLOT.png']);
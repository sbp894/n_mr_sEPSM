% %%
% clear;
% close all;
% clc;

%%
function create_summary(ChinID)

DataRepository=[fileparts(pwd) '\Output\DataAnal\'];

temp=dir([DataRepository '*' num2str(ChinID) '*']);
DataDir=strcat(DataRepository,temp(end).name,filesep);% cs=load([DataDir 'conditions.mat']);

% csCell=cs.csCell;
SACDir=[DataDir 'sac\'];


%%
allFiles=dir([SACDir '*.mat']);
SummaryVar=repmat(struct('BF_TC',nan),length(allFiles),1);
ModFreq = [1,2,4,8,16,32,64];

%% Read Unit Files
MATDataDir=[fileparts(pwd) '\MATData\' DataDir(length(fileparts(fileparts((DataDir))))+2:end-3) filesep];
unitFiles=dir([MATDataDir 'Unit_*.mat']);
tIgnore=.1;

calibFName=dir([MATDataDir '*calib*']);
CalibData=load([MATDataDir calibFName(end).name]);
CalibData=CalibData.data.CalibData(:,1:2);
TFiltWidthTC=5;
lw=1.5;

% CF_SR_SAT_Q10_var=nan(length(unitFiles),4);
CF_SR_SAT_Q10_1=nan(length(unitFiles),1);
CF_SR_SAT_Q10_2=nan(length(unitFiles),1);
CF_SR_SAT_Q10_3=nan(length(unitFiles),1);
CF_SR_SAT_Q10_4=nan(length(unitFiles),1);

parfor unit_var=1:length(unitFiles)
    data1=load([MATDataDir unitFiles(unit_var).name]);
    data1=data1.data;
    CF_SR_SAT_Q10_1(unit_var)=data1.BFmod;
    
    SRFname=[MATDataDir  sprintf('*_u%1d_%02d_SR*',data1.track,data1.No)];
    SNRenvFname=[MATDataDir  sprintf('*_u%1d_%02d_SNRenv*',data1.track,data1.No)];
    RLFFname=[MATDataDir  sprintf('*_u%1d_%02d_RLV*',data1.track,data1.No)];
    TCFname=[MATDataDir  sprintf('*_u%1d_%02d_tc*',data1.track,data1.No)];
    
    if ~isempty(dir(SRFname))
        temp=dir(SRFname);
        data=load([MATDataDir temp(1).name]);
        data=data.data;
        if isempty(data.spikes{1,1})
            CF_SR_SAT_Q10_2(unit_var)=0;
            fprintf('The spont rate should be 0 for %s\n',SRFname);
        else
            CF_SR_SAT_Q10_2(unit_var)=size(data.spikes{1,1},1)/max(data.spikes{1,1}(:,1));
        end
    elseif ~isempty(dir(SNRenvFname))
        temp=dir(SNRenvFname);
        data=load([MATDataDir temp(1).name]);
        data=data.data;
        CF_SR_SAT_Q10_2(unit_var)=size(data.spikes{1,1},1)/max(data.spikes{1,1}(:,1));
        spkData=data.spikes{1,1};
        tStimOffStart=(data.Hardware.Trigger.StmOn+tIgnore*1e3)/1e3;
        spkData=spkData(spkData(:,1)<=data.Stimuli.fully_presented_stimuli,2);
        CF_SR_SAT_Q10_2(unit_var)=sum(spkData>tStimOffStart)/data.Stimuli.fully_presented_stimuli/(data.Hardware.Trigger.StmOff/1e3-tIgnore);
    end
    
    if ~isempty(dir(RLFFname))
        temp=dir(RLFFname);
        data=load([MATDataDir temp(1).name]);
        data=data.data;
        spkData=data.spikes{1,1};
        dBs2run=nan(data.Stimuli.fully_presented_lines,1);
        rates=nan(data.Stimuli.fully_presented_lines,1);
        
        %%
        PlotRLV=0;
        plotFittedRLV=0;
        %%
        
        for line_var=1:length(rates)
            dBs2run(line_var)=line_var; % actual value doesn't matter. We are interested in sat rate
            spkCurLine=spkData(spkData(:,1)==line_var,2);
            if isempty(spkCurLine)
                rates(line_var)=0;
            else
                spkCurLine=spkCurLine(spkCurLine<data.Hardware.Trigger.StmOn/1e3);
                rates(line_var)=length(spkCurLine)/(data.Hardware.Trigger.StmOn/1e3);
            end
        end
        [RLVparams,~,~]=NELfuns.fitRLfun(dBs2run,rates,PlotRLV,plotFittedRLV);
        CF_SR_SAT_Q10_3(unit_var)=RLVparams.R_Sat;
        %         pause;
        %         clf;
        
    else
        CF_SR_SAT_Q10_3(unit_var)=nan;
    end
    
    if ~isempty(dir(TCFname))
        temp=dir(TCFname);
        TCdata=load([MATDataDir temp(1).name]);
        TCdata=TCdata.data.TcData;
        TCdata=TCdata(TCdata(:,1)~=0,:); % Get rid of all 0 freqs
        for i=1:size(TCdata,1)
            TCdata(i,3)=NELfuns.CalibInterp(TCdata(i,1),CalibData)-TCdata(i,2);
        end
        TCdata(:,4)=Library.trifilt(TCdata(:,3)',TFiltWidthTC)';
        [CF_SR_SAT_Q10_4(unit_var),~,~,~] = NELfuns.findQ10(TCdata(:,1),TCdata(:,4),CF_SR_SAT_Q10_1(unit_var));
    else
        CF_SR_SAT_Q10_4(unit_var)=nan;
    end
    
end

CF_SR_SAT_Q10_var=[CF_SR_SAT_Q10_1, CF_SR_SAT_Q10_2, CF_SR_SAT_Q10_3, CF_SR_SAT_Q10_4];
SpikeStimulusData=load([DataDir 'SpikeStimulusData.mat']);

for file_var=1:size(SummaryVar,1)
    curFName=allFiles(file_var).name;
    %     curSACMAT=load([SACDir curFName]);
    
    %% Read and Assign CF, SNR, dB
    curParamsIN=load([DataDir 'paramsIN\paramsIN' curFName(4:end)]);
    curParamsIN=curParamsIN.paramsIN;
    SummaryVar(file_var).fName=curFName;
    SummaryVar(file_var).DataDir=DataDir;
    SummaryVar(file_var).ChinID=ChinID;
    SummaryVar(file_var).BF_TC=curParamsIN.CF_A_Hz;
    SummaryVar(file_var).SNR=curParamsIN.SNR2use_dB;
    SummaryVar(file_var).SPL=curParamsIN.level;
    SummaryVar(file_var).fName=curFName(7:end);
    %% Assign Spont Rate
    SummaryVar(file_var).SR=nanmean(CF_SR_SAT_Q10_var(CF_SR_SAT_Q10_var(:,1)==curParamsIN.CF_A_Hz/1e3,2));
    SummaryVar(file_var).SatR=nanmean(CF_SR_SAT_Q10_var(CF_SR_SAT_Q10_var(:,1)==curParamsIN.CF_A_Hz/1e3,3));
    SummaryVar(file_var).Q10_TC=nanmean(CF_SR_SAT_Q10_var(CF_SR_SAT_Q10_var(:,1)==curParamsIN.CF_A_Hz/1e3,4));
    %     SummaryVar(file_var).SR=nanmean(CF_SR_var(((CF_SR_var(:,1)-curParamsIN.CF_A_Hz/1e3)/curParamsIN.CF_A_Hz*1e6)<1,2));  % Don't know why didn't work.
    
    %% Assign SNRenv
    
    temp_psd=load([DataDir '\psd\psd' curFName(4:end)]);
    
    ratio_mat=((temp_psd.PowerMod_SN-temp_psd.PowerMod_N)./temp_psd.PowerMod_N);
    ratio_mat(ratio_mat<0)=eps;
    SummaryVar(file_var).SNRenvAll=20*log10(sqrt(sum(ratio_mat.^2,2)));
    SummaryVar(file_var).SNRenv=mean(SummaryVar(file_var).SNRenvAll);
    
    ratio_mat_nf=((temp_psd.PowerMod_SN_noisefloor-temp_psd.PowerMod_N_noisefloor)./temp_psd.PowerMod_N_noisefloor);
    ratio_mat_nf(ratio_mat_nf<0)=eps;
    SummaryVar(file_var).SNRenvAll_NF=20*log10(sqrt(sum(ratio_mat_nf.^2,2)));
    SummaryVar(file_var).SNRenv_NF=mean(SummaryVar(file_var).SNRenvAll_NF);
    
    %% Plot with inset
    figure(21); clf; hold on; set (gcf, 'Units', 'normalized', 'Position', [.1,.1,.8,.8]);
    set(gcf,'visible','off');
    errorbar(mean(temp_psd.PowerMod_S), std(temp_psd.PowerMod_S),'b','linewidth',lw);
    errorbar(mean(temp_psd.PowerMod_SN), std(temp_psd.PowerMod_SN),'k','linewidth',lw);
    errorbar(mean(temp_psd.PowerMod_N), std(temp_psd.PowerMod_N),'r','linewidth',lw);
    
%     errorbar(mean(temp_psd.PowerMod_S_noisefloor), std(temp_psd.PowerMod_S_noisefloor),'b:','linewidth',lw);
%     errorbar(mean(temp_psd.PowerMod_SN_noisefloor), std(temp_psd.PowerMod_SN_noisefloor),'k:','linewidth',lw);
%     errorbar(mean(temp_psd.PowerMod_N_noisefloor), std(temp_psd.PowerMod_N_noisefloor),'r:','linewidth',lw);
    
    axis tight; set(gca, 'xticklabels', ModFreq); xlabel('Mod Freq'); ylabel('Mod Power'); 
    title(sprintf('%s (SNRenv =%f|| SNRenvNF= %f dB)',strrep(curFName, '_', ':'), SummaryVar(file_var).SNRenv, SummaryVar(file_var).SNRenv_NF));
%     legend('S','N','SN','NF[S]','NF[N]','NF[SN]','location','northwestoutside');
    legend('S','SN','N','location','northwestoutside');
        
    axes('Position',[.7 .7 .2 .2]); box on; hold on;
    errorbar(mean(ratio_mat), std(ratio_mat),'linewidth',lw);
    errorbar(mean(ratio_mat_nf), std(ratio_mat_nf),'linewidth',lw);
    legend('data', 'nf', 'location', 'best'); 
    xlabel('Mod Freq'); ylabel('SNR_e_n_v'); 
    title('Ratio (not dB)');
    axis tight; set(gca, 'xticklabels', ModFreq);
    figName=[DataDir 'envPowereps\SummaryFig_' strrep(curFName(4:end-4),'.','_')];
    saveas(gcf,figName);
    saveas(gcf,figName,'bmp');
    close (21);
    %%
    if isempty(find((temp_psd.PowerMod_SN-temp_psd.PowerMod_SN_noisefloor)./temp_psd.PowerMod_SN_noisefloor <=0, 1))
        %      if isempty(find(sum((temp_psd.PowerMod_SN-temp_psd.PowerMod_SN_noisefloor)./temp_psd.PowerMod_SN_noisefloor,1) <=0, 1))
        SummaryVar(file_var).AboveNFStrict=1;
    else
        SummaryVar(file_var).AboveNFStrict=0;
    end
    
    %%
    temp_params=temp_psd.paramsFile;
    nLinesSN=0;
    nLinesN=0;
    nSpikesSN=0;
    nSpikesN=0;
    spkRateSN=0;
    spkRateN=0;
    for parm_var=1:length(temp_params)
        nLinesSN=nLinesSN+mean(temp_params(parm_var).nLines(2,:));
        nLinesN=nLinesSN+mean(temp_params(parm_var).nLines(3,:));
        nSpikesSN=nSpikesSN+mean(temp_params(parm_var).NumDrivenSpikes(2,:));
        nSpikesN=nSpikesN+mean(temp_params(parm_var).NumDrivenSpikes(3,:));
        spkRateSN=spkRateSN+mean(temp_params(parm_var).AvgRate_sps(2,:));
        spkRateN=nSpikesN+mean(temp_params(parm_var).AvgRate_sps(3,:));
    end
    
    SummaryVar(file_var).nLinesSN=nLinesSN/length(temp_params);
    SummaryVar(file_var).nLinesN=nLinesN/length(temp_params);
    SummaryVar(file_var).nSpikesSN=nSpikesSN/length(temp_params);
    SummaryVar(file_var).nSpikesN=nSpikesN/length(temp_params);
    SummaryVar(file_var).spkRateSN=spkRateSN/length(temp_params);
    SummaryVar(file_var).spkRateN=spkRateN/length(temp_params);
end
SummaryVar=NELfuns.get_revcor(SpikeStimulusData,SummaryVar);

save([DataDir 'SummaryVar.mat'], 'SummaryVar');

%%
tx=extractfield(SummaryVar, 'SNRenv');
ty=extractfield(SummaryVar, 'SNRenvACST');
tz=extractfield(SummaryVar, 'AboveNFStrict')==1;
scatter(tx(tz),ty(tz))
% scatter(tx,ty);
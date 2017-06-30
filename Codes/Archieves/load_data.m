function [spike_data,StimData,StimsFNames]=load_data(DataDir,resultsDir)

%% To return: These are the things that will be used later.
% Sampling Frequency, presentation level, condition_table, spike data corresponding
% to each condition in the condition table, paramsIN.durA_msec=dur_sec*1000;
% paramsIN.durB_msec=dur_sec*1000;
% paramsIN.durC_msec=dur_sec*1000;
% paramsIN.CF_A_Hz=CF_A_kHz*1000;
% paramsIN.CF_B_Hz=CF_B_kHz*1000;
% paramsIN.CF_C_Hz=CF_C_kHz*1000;
% % Need to include CF_A, CF_B, CF_C for more generality
% paramsIN.MAXspikes=2500;
% paramsIN.SNR2use_dB=SNR2use_dB;
% %% *MH July 7 2015 - change to longer delays to allow for lower modulations
% paramsIN.MAXdelay_sec=1;  % sentence duration is 2.7 sec, so 1 sec delays are about 1/3 of stim duration, which is pushing what we can estimate

CurDir=pwd;
addpath(CurDir);

%% load data
cd (DataDir);
allfiles=dir('*SNRenv.mat');

csCell=[];
spike_data=struct([]);
StimData={};
StimsFNames={};
Ncases=0;

allCalibfiles=dir('*calib*');
fprintf('Using %s as calibration file\n', allCalibfiles(end).name);
CalibData=load(allCalibfiles(end).name);
CalibData=CalibData.data.CalibData;

for file_var=1:length(allfiles)
    filename=allfiles(file_var).name;
    PicData=load(filename);
    PicData=PicData.data;
    
    unitNUM=PicData.General.unit;
    trackNUM=PicData.General.track;
    stim_list=PicData.Stimuli.list';
    %     eval(sprintf('x=Unit_%d_%02d;',trackNUM,unitNUM));
    x=load(sprintf('Unit_%d_%02d.mat',trackNUM,unitNUM));
    x=x.data;
    
    paramsIN.durA_msec=PicData.Hardware.Trigger.StmOn;
    paramsIN.durB_msec=PicData.Hardware.Trigger.StmOn;
    paramsIN.durC_msec=PicData.Hardware.Trigger.StmOn;
    paramsIN.MAXspikes=4300;
    
    TypeSNRpolarityMat=zeros(length(stim_list),3);
    for list_var=1:length(stim_list)
        cur_stim_name=stim_list{list_var};
        cur_stim_name=cur_stim_name(strfind(cur_stim_name,'Stim'):end);
        
        if strfind(cur_stim_name,'_N_')
            TypeSNRpolarityMat(list_var,1)=-1;
        elseif strfind(cur_stim_name,'_SN_')
            TypeSNRpolarityMat(list_var,1)=0;
        elseif strfind(cur_stim_name,'_S_')
            TypeSNRpolarityMat(list_var,1)=1;
        end
        
        
        if ~isempty(sscanf(cur_stim_name,'Stim%f*dB_%*s_%*s.wav'))
            TypeSNRpolarityMat(list_var,2)=sscanf(cur_stim_name,'Stim%f*dB_%*s_%*s.wav');
        end
        
        if strcmp(cur_stim_name(end-5),'P')
            TypeSNRpolarityMat(list_var,3)=1;
        else
            TypeSNRpolarityMat(list_var,3)=-1;
        end
    end
    
    SNRs_for_current_unit=unique(TypeSNRpolarityMat(TypeSNRpolarityMat(:,1)==0,2));
    csCell=[csCell;zeros(length(SNRs_for_current_unit),2)]; %#ok<*AGROW>
    csCell(end-length(SNRs_for_current_unit)+1:end,2)=SNRs_for_current_unit';
    csCell(end-length(SNRs_for_current_unit)+1:end,1)= 1e3*x.BFmod; %% Update
    
    
    
    %     paramsIN.level=CalibInterp(csCell(1,1)/1e3,CalibData)-PicData.Stimuli.attens;
    %     %%% If we are interested in energy at CF, else the value below
    %     should be alright.
    paramsIN.level=120-PicData.Stimuli.attens;
    
    for snr_var=1:length(SNRs_for_current_unit)
        Ncases=Ncases+1;
        [spike_data(end+1).SpikeTrains,StimData{end+1},StimsFNames{end+1}]=DataAnal.get_SpikeTrains(PicData,stim_list,TypeSNRpolarityMat,SNRs_for_current_unit(snr_var));
        spike_data(end).nReps=mean(mean(cellfun(@(x) numel(x),spike_data(end).SpikeTrains)));
        spike_data(end).CF=csCell(Ncases,1);
        spike_data(end).SNR=csCell(Ncases,2);
        spike_data(end).Fs=PicData.Stimuli.updateRate_Hz;
        paramsIN.CF_A_Hz=csCell(Ncases,1);
        paramsIN.CF_B_Hz=csCell(Ncases,1);
        paramsIN.CF_C_Hz=csCell(Ncases,1);
        paramsIN.SNR2use_dB=csCell(Ncases,2);
        paramsIN.MAXdelay_sec=floor(paramsIN.durA_msec/1e3);
        paramsIN.plt=1;
        
        paramsIN.FiberType=2;
        spike_data(end).paramsIN=paramsIN;
        spike_data(end).SPL=paramsIN.level;
    end
end

cd (CurDir);
csCell(:,3)=csCell(:,2);
csCell(:,2)=0*csCell(:,2)+1;
csCell=num2cell(csCell);

for i=1:size(csCell,1)
    csCell{i,4}='SSN';
    csCell{i,5}=spike_data(i).paramsIN.level;
end

csLabels = {'CF','sentences','snr','noise','level'}; %#ok<*NASGU>
save([resultsDir 'conditions.mat'],'csCell','csLabels');
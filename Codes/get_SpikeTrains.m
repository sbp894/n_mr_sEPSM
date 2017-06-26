%%% Should be added to path
function [SpikeTrains,StimsData,StimsFNames]=get_SpikeTrains(PicData,stim_list,TypeSNRpolarityMat,curSNR)

SpikeTrains=cell(3,2);
% StimsData=cell(3,2);
all_spikes_times=PicData.spikes{1};
all_files_played=PicData.Line.file;
all_files_played=all_files_played(1:PicData.Stimuli.fully_presented_lines);

% Speech 0dB Positive 
stim_name=stim_list(ismember(TypeSNRpolarityMat,[1 0 1],'rows'));
SpikeTrainsS_plus={};

% [~,DataFolderName,~]=fileparts(pwd);
% UserName=DataFolderName(1:2);
UserName='MH';

[~,stim_name2load,~]=fileparts(stim_name{1});
StimsS_plus_fName=[pwd filesep 'Signals' filesep UserName filesep 'SNRenv' filesep stim_name2load '.wav'];
SpeechStimsS_plus=audioread(StimsS_plus_fName);

for file_var=1:length(all_files_played)
    if strcmp(all_files_played{file_var},stim_name)
        SpikeTrainsS_plus{end+1}=(sort(all_spikes_times((all_spikes_times(:,1)==file_var),2)));
    end
    
end

% Speech 0dB Negative 
stim_name=stim_list(ismember(TypeSNRpolarityMat,[1 0 -1],'rows'));
SpikeTrainsS_minus={};

[~,stim_name2load,~]=fileparts(stim_name{1});
StimsS_minus_fName=[pwd filesep 'Signals' filesep UserName filesep 'SNRenv' filesep stim_name2load '.wav'];
SpeechStimsS_minus=audioread(StimsS_minus_fName);

for file_var=1:length(all_files_played)
    if strcmp(all_files_played{file_var},stim_name)
        SpikeTrainsS_minus{end+1}=(sort(all_spikes_times((all_spikes_times(:,1)==file_var),2)));
    end
    
end

% Noise SNR dB Positive 
stim_name=stim_list(ismember(TypeSNRpolarityMat,[-1 curSNR 1],'rows'));
SpikeTrainsN_plus={};

[~,stim_name2load,~]=fileparts(stim_name{1});
StimsN_plus_fName=[pwd filesep 'Signals' filesep UserName filesep 'SNRenv' filesep stim_name2load '.wav'];
SpeechStimsN_plus=audioread(StimsN_plus_fName);

for file_var=1:length(all_files_played)
    if strcmp(all_files_played{file_var},stim_name)
        SpikeTrainsN_plus{end+1}=(sort(all_spikes_times((all_spikes_times(:,1)==file_var),2)));
    end
    
end

% Noise SNR dB Negative 
stim_name=stim_list(ismember(TypeSNRpolarityMat,[-1 curSNR -1],'rows'));
SpikeTrainsN_minus={};

[~,stim_name2load,~]=fileparts(stim_name{1});
StimsN_minus_fName=[pwd filesep 'Signals' filesep UserName filesep 'SNRenv' filesep stim_name2load '.wav'];
SpeechStimsN_minus=audioread(StimsN_minus_fName);

for file_var=1:length(all_files_played)
    if strcmp(all_files_played{file_var},stim_name)
        SpikeTrainsN_minus{end+1}=(sort(all_spikes_times((all_spikes_times(:,1)==file_var),2)));
    end
    
end

% Speech+Noise SNR dB Positive 
stim_name=stim_list(ismember(TypeSNRpolarityMat,[0 curSNR 1],'rows'));
SpikeTrainsSN_plus={};

[~,stim_name2load,~]=fileparts(stim_name{1});
StimsSN_plus_fName=[pwd filesep 'Signals' filesep UserName filesep 'SNRenv' filesep stim_name2load '.wav'];
SpeechStimsSN_plus=audioread(StimsSN_plus_fName);

for file_var=1:length(all_files_played)
    if strcmp(all_files_played{file_var},stim_name)
        SpikeTrainsSN_plus{end+1}=(sort(all_spikes_times((all_spikes_times(:,1)==file_var),2)));
    end
    
end

% Speech+Noise SNR dB Negative 
stim_name=stim_list(ismember(TypeSNRpolarityMat,[0 curSNR -1],'rows'));
SpikeTrainsSN_minus={};

[~,stim_name2load,~]=fileparts(stim_name{1});
StimsSN_minus_fName=[pwd filesep 'Signals' filesep UserName filesep 'SNRenv' filesep stim_name2load '.wav'];
SpeechStimsSN_minus=audioread(StimsSN_minus_fName);

for file_var=1:length(all_files_played)
    if strcmp(all_files_played{file_var},stim_name)
        SpikeTrainsSN_minus{end+1}=(sort(all_spikes_times((all_spikes_times(:,1)==file_var),2)));
    end
    
end

% SpikeTrains=[SpikeTrainsS_plus,SpikeTrainsS_minus;SpikeTrainsN_plus,SpikeTrainsN_minus;SpikeTrainsSN_plus,SpikeTrainsSN_minus];

SpikeTrains{1,1}=SpikeTrainsS_plus;
SpikeTrains{1,2}=SpikeTrainsS_minus;
SpikeTrains{2,1}=SpikeTrainsN_plus;
SpikeTrains{2,2}=SpikeTrainsN_minus;
SpikeTrains{3,1}=SpikeTrainsSN_plus;
SpikeTrains{3,2}=SpikeTrainsSN_minus;

StimsData={{SpeechStimsS_plus},{SpeechStimsS_minus};{SpeechStimsN_plus},{SpeechStimsN_minus};{SpeechStimsSN_plus},{SpeechStimsSN_minus}};
StimsFNames={{StimsS_plus_fName},{StimsS_minus_fName};{StimsN_plus_fName},{StimsN_minus_fName};{StimsSN_plus_fName},{StimsSN_minus_fName}};
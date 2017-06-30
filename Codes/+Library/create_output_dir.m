function root_dir=create_output_dir(Simulation1DataAnal0,resultDirName, RootOUTPUTDir)

if Simulation1DataAnal0==1
    SwitchDir='Simulation';
elseif Simulation1DataAnal0==0
    SwitchDir='DataAnal';
elseif Simulation1DataAnal0==-1
    SwitchDir='Debug';
end

allfiles=dir([RootOUTPUTDir SwitchDir]);
allfiles=allfiles([allfiles(:).isdir]);

numDIR=0;
for i=1:length(allfiles)
   if strfind(allfiles(i).name,resultDirName)
       numDIR=numDIR+1;
   end
end
numDIR=max(numDIR,1);

root_dir=[RootOUTPUTDir SwitchDir filesep resultDirName '-' num2str(numDIR) filesep]; 

if exist(root_dir,'dir')
    choice = questdlg(sprintf('The result directory- (%s) exists for this data. What should be done?', root_dir),'Already Analyzed Data', 'Create_new_directory','Resume','Replace','Resume');
    switch choice
        case 'Create_new_directory'
                root_dir=[RootOUTPUTDir SwitchDir filesep resultDirName '-' num2str(numDIR+1) filesep];
        case 'Resume'
            % Resumed
        case 'Replace'
            rmdir(root_dir,'s');
    end
    
end

if ~exist(root_dir,'dir')
    mkdir(root_dir);
    mkdir([root_dir 'eps']);
    mkdir([root_dir 'png']);
    mkdir([root_dir 'psd']);
    mkdir([root_dir 'sac']);
    mkdir([root_dir 'sacMET']);
    mkdir([root_dir 'paramsIN']);
    mkdir([root_dir 'envPowerpng']);
    mkdir([root_dir 'envPowereps']);
    mkdir([root_dir 'ModEP']);
    mkdir([root_dir 'SpeechSpikePlots_bmp']);
    mkdir([root_dir 'SpeechSpikePlots_eps']);
end

if ~exist([root_dir 'progress'],'dir')
    mkdir([root_dir 'progress']);
end
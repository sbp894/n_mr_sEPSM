function parse_saved_data_for_SCDCPeaks_Vijay(SACMETDir)

if ~strcmp(SACMETDir(end),filesep)
    SACMETDir=[SACMETDir filesep];
end

mkdir([SACMETDir 'SCDCPeaks']);

allfiles=dir([SACMETDir '*.mat']);
SCDCPeakMat=zeros(length(allfiles),9);

for filevar=1:length(allfiles)
    curfilename=allfiles(filevar).name;
    if ~strcmp(curfilename(1),'.')
        CurSACMet=load([SACMETDir curfilename]);
        CurStrings=strsplit(curfilename,'_');
        SACSCCmetricsAVG=CurSACMet.SACSCCmetrics{end};
        
        %% 
        CFkHz=sscanf(CurStrings{2},'%f');
        dBSPL=str2double(CurStrings{7});
        SNR=str2double(CurStrings{8}(4:end-4));
        
        %%
        SCpeakS=SACSCCmetricsAVG.SCpeaks_A(2);
        SCpeakN=SACSCCmetricsAVG.SCpeaks_B(2);
        SCpeakSN=SACSCCmetricsAVG.SCpeaks_C(1); %% To change
        
        DCpeakS=SACSCCmetricsAVG.DCpeak_A(1);
        DCpeakN=SACSCCmetricsAVG.DCpeak_B(1);
        DCpeakSN=SACSCCmetricsAVG.DCpeak_C(1);        
        
        SCDCPeakMat(filevar,:)=[CFkHz,dBSPL,SNR,SCpeakS,SCpeakN,SCpeakSN,DCpeakS,DCpeakN,DCpeakSN];
        
    end
end

readme={{'cf'},{'dB'},{'SNR'},{'SC_S'},{'SC_N'},{'SC_SN'},{'DC_S'},{'DC_N'},{'DC_SN'}};
save([SACMETDir 'SCDCPeaks' filesep 'SCDCPeakMat.mat'],'SCDCPeakMat');
save([SACMETDir 'SCDCPeaks' filesep 'readme.mat'],'readme');
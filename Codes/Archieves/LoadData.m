function [SpikeTrainsC_plus,SpikeTrainsC_minus,paramsIn]=LoadData(DataDir,curFileName)

CurDir=pwd;
cd(DataDir);
PicData=loadPic(str2double(curFileName(2:5)));

paramsIn.unitNum=PicData.General.unit;
paramsIn.trackNum=PicData.General.track;
paramsIn.trigger=PicData.General.trigger;
paramsIn.nrep=PicData.Stimuli.fully_presented_lines;
paramsIn.T=(PicData.Hardware.Trigger.StmOn+PicData.Hardware.Trigger.StmOff)/1e3;

if paramsIn.unitNum<10
    CFFile=(['Unit_' num2str(paramsIn.trackNum) '_0' num2str(paramsIn.unitNum)]);
else
    CFFile=(['Unit_' num2str(paramsIn.trackNum) '_' num2str(paramsIn.unitNum)]);
end
UnitData=eval('loadPic(%s)',CFFile); %#ok<EVLC>
paramsIn.CF=UnitData.BF*1e3;

spikedata=PicData.spikes{1};

SpikeTrainsC_plus=cell(ceil(paramsIn.nrep/2),1);
SpikeTrainsC_minus=cell(paramsIn.nrep-length(SpikeTrainsC_plus),1);

for i=1:paramsIn.nrep
    if rem(i,2)==1
        SpikeTrainsC_plus{(i+1)/2}=spikedata(spikedata(:,1)==i,2);
    else
        SpikeTrainsC_minus{i/2}=spikedata(spikedata(:,1)==i,2);
    end
end

cd(CurDir);

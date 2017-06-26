function [mrSPikeTrains, paramsIN, PowerModCell]=get_mrSPikeTrains(SpikeTrains, mrWindows, onsetIgnore, paramsIN)
dur_whole=paramsIN.durA_msec/1e3;
mrSPikeTrains=cell(length(mrWindows),1);
PowerModCell=cell(length(mrWindows),1);
mrTimeBoundaries=cell(length(mrWindows),1);


for windowVar=1:length(mrWindows)
    curWindow=mrWindows(windowVar);
    if windowVar==1
        curCellSpk= cell(1, 2); %hardcoded to consider 700ms as another window for 1 Hz modulation
        curCellBnd= cell(1, 2);
        PowerModCell{windowVar}=cell(1, 2);
    else 
        curCellSpk= cell(1, floor((dur_whole-onsetIgnore)/curWindow));
        curCellBnd= cell(1, floor((dur_whole-onsetIgnore)/curWindow));
        PowerModCell{windowVar}=cell(1, floor((dur_whole-onsetIgnore)/curWindow));
    end
    for curLoop=1:length(curCellSpk)
        tempCellSpk=cell(size(SpikeTrains));
        for rowVar=1:size(SpikeTrains,1)
            for colVar=1:size(SpikeTrains,2)
                tempCellSpk{rowVar, colVar}=cellfun(@(x) mrWindowSpktimes(x,onsetIgnore+(curLoop-1)*curWindow, onsetIgnore+curLoop*curWindow), SpikeTrains{rowVar, colVar}, 'Un', 0);
            end
        end
        curCellSpk{curLoop}=tempCellSpk;
        curCellBnd{curLoop}=[onsetIgnore+(curLoop-1)*curWindow, min(onsetIgnore+curLoop*curWindow, dur_whole)];
    end
    mrSPikeTrains{windowVar}=curCellSpk;
    mrTimeBoundaries{windowVar}=curCellBnd;
end
paramsIN.mrTimeBoundaries=mrTimeBoundaries;
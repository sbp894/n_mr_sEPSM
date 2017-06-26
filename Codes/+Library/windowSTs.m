function [SpikeTrain_new,Nspikes,nLines]=windowSTs(SpikeTrain_old,start_sec,end_sec,MAXspikes)
% File: [SpikeTrain_new,Nspikes]=windowSTs(SpikeTrain_old,start_sec,end_sec,MAXspikes)
%
% WINDOWS SpikeTrains given starting and ending time.
% Assumes spikes are in CCC-spiketrain format [cell_array{Nreps}]
SpikeTrain_old=SpikeTrain_old(~cellfun('isempty',SpikeTrain_old));

if isempty(SpikeTrain_old)
    Nspikes=0; SpikeTrain_new={}; nLines=0;
else
    Nspikes=0;
    SpikeTrain_new=cell(size(SpikeTrain_old));
    for i=1:length(SpikeTrain_old)
        if ~isempty(SpikeTrain_old{i})
            spikeINDs=find((SpikeTrain_old{i}>=start_sec) & (SpikeTrain_old{i}<=end_sec));
            if (Nspikes+length(spikeINDs))<=MAXspikes
                Nspikes=Nspikes+length(spikeINDs);
                SpikeTrain_new{i} = SpikeTrain_old{i}(spikeINDs);
                nLines=i;
            else
                break
            end
        end
    end
    
    SpikeTrain_new=SpikeTrain_new(~cellfun('isempty',SpikeTrain_new));
end

return;
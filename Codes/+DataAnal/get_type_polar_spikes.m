% Use cellfun instead of unallocated spktimes
function spktimes=get_type_polar_spikes(SpikeTrains)

spktimes=[];
for i=1:length(SpikeTrains)
    spktimes=[spktimes;SpikeTrains{i}];
end
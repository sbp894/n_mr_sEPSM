function avgRate = calculateFiringRate(SpTimes,anal)

%anal.onsetIgnore = 0.05;
%anal.maxSpikes = 3600;

ST = Library.windowSTs(SpTimes,anal.onsetIgnore,anal.duration,anal.maxSpikes);

numSpikeReps = length(ST);
% Compute AVGrate
kMax=0;
for spikeIdx= 1 : numSpikeReps
    if length(ST{spikeIdx}) > kMax
        kMax = length(ST{spikeIdx});
    end
end

spikeMat=nan(numSpikeReps,kMax); % Setup BIG spike matrix for faster computation

for spikeIdx= 1 : numSpikeReps
    spikeMat(spikeIdx,1 : length(ST{spikeIdx})) = ST{spikeIdx};
end

numSpikes=sum(~isnan(spikeMat(:,:))');  % Count number of real spikes in each line
totalSpikes=sum(numSpikes);
avgRate=totalSpikes/numSpikeReps/(anal.duration-anal.onsetIgnore);

end
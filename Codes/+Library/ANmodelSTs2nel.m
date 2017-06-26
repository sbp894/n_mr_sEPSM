function [NELspikes,nspikes]=ANmodelSTs2nel(sptimes,Nreps)
% File: [NELSpikeTrains,Nspikes]=ANmodelSTs2nel(sptimes)
%
% CONVERT SpikeTrains in ANmodel row format: e.g., ARLO, ZB06/07 to 
% NEL-spiketrain format [rep #,spiketime_sec]

nspikes=length(sptimes);

NELspikes=NaN*ones(nspikes,2);
NELspikes(:,2)=sptimes;
REPendINDs(1)=0;
tempINDs=find(diff(sptimes)<0);  % Allows for no spikes in some reps, by using length(tempINDs) to provide Nreps
Nreps2=length(tempINDs)+1;
if Nreps2<Nreps
    fprintf('***WARNING***NOTE: Nreps2 (%d) < Nreps (%d) - Not all reps had spikes here',Nreps2,Nreps);
end
REPendINDs(2:Nreps2)=tempINDs;
REPendINDs(Nreps2+1)=nspikes;
for REPind=1:Nreps2
	NELspikes(REPendINDs(REPind)+1:REPendINDs(REPind+1),1)=REPind;
end


return;


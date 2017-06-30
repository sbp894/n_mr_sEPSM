% function histbar=get_type_polar_spikes(SpikeTrains)
function spktimes=get_type_polar_spikes(SpikeTrains)

spktimes=[];
for i=1:length(SpikeTrains)
    spktimes=[spktimes;SpikeTrains{i}];
end

% spktimes=sort(spktimes);
% 
% unique_spktimes=unique(spktimes);
% histbar=zeros(1+3*length(unique_spktimes),2);
% 
% for i=1:length(unique_spktimes)
%     histbar(3*i-1,1)=unique_spktimes(i)-1e-5;
%     histbar(3*i  ,1)=unique_spktimes(i);
%     histbar(3*i+1,1)=unique_spktimes(i)+1e-5;
%     histbar(3*i,2)=sum(spktimes==histbar(i,1));
% end

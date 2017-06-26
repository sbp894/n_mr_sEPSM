function histbar=spktimes_to_histbar(spktimes)

spktimes=sort(spktimes);
histbar=zeros(length(unique(spktimes)),2);
histbar=unique(spktimes);

for i=1:length(spktimes)-1
    if spktimes(i)==spktimes(i+1)
        
    end
end
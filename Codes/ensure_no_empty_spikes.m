function outcell=ensure_no_empty_spikes(meanrate_unad,dur_sec,ANmodel_Fs_Hz,Nreps, anal)

count=1;
outcell=cell(Nreps,1);

while count<=Nreps
    SpikeTrain=get_sptimes(meanrate_unad,ANmodel_Fs_Hz,1);
    temp=SpikeTrain{1};
    temp(temp<anal.onsetIgnore)=[];
    temp(temp>dur_sec)=[];
    if ~isempty(temp)
        outcell{count}=temp;
        count=count+1;
    end
end

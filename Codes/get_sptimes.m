function sptimes=get_sptimes(meanrate_unad, Fs, nrep)

sptimes=cell(nrep,1);

for i=1:nrep
    sptimes{i} = Library.SGfast([1/Fs, 1],meanrate_unad)';
end

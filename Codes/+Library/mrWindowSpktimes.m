function spiketimes=mrWindowSpktimes(spiketimes, tStart, tEnd)

spiketimes(spiketimes<tStart)=[];
spiketimes(spiketimes>tEnd)=[];
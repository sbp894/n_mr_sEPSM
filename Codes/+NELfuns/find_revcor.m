function [h,Nspikes]=find_revcor(pin,spike_times,Fs,window_dur)

Fsr=20e3;
pin=resample(pin,Fsr,Fs);

spike_times=spike_times(spike_times>window_dur+0.001);
spike_times=spike_times(spike_times<length(pin)/Fsr);


h=zeros(length(spike_times),window_dur*Fsr+1);
% length(pin)
for i=1:length(spike_times)
    %     if size(h,2)==length((spike_times(i)-window_dur)*Fs:spike_times(i)*Fs)
    %     [(spike_times(i)-window_dur)*Fs,spike_times(i)*Fs]
    
    pinInd=round((spike_times(i)-window_dur)*Fsr:spike_times(i)*Fsr);
    if length(pinInd)==size(h,2)
        h(i,:)=pin(pinInd);
    elseif length(pinInd)>size(h,2)
        h(i,:)=pin(pinInd(1:end-1));
    else 
        h(i,:)=pin([pinInd(1)-1 pinInd]);
    end
    %     else
    %         disp('wait! ');
    %     end
end
h=resample(mean(h)',Fs,Fsr);
Nspikes=length(spike_times);
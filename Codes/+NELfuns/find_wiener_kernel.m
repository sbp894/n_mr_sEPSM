function [k2,Nspikes]=find_wiener_kernel(pin,spike_times,Fs,window_dur,nRep)

spike_times=sort(spike_times);

if size(pin,2)~=1
    pin=pin';
end

spkWindowLow=window_dur+0.001;
spkWindowHi=length(pin)/Fs;

spike_times=spike_times(spike_times>spkWindowLow); % Ignore spikes with latency< 1 ms
spike_times=spike_times(spike_times<spkWindowHi);

Nspikes=length(spike_times);

k2_diag_correction=var(pin)*length(spike_times)/nRep/(length(pin)/Fs)/2/std(pin)^4;


%% Expensive
% k2=zeros(window_dur*Fs+1,window_dur*Fs+1);
% 
% for tau1=1:size(k2,1)
%     for tau2=tau1:size(k2,2)
%         tindClosest=dsearchn(spike_times',tau2/Fs);
%         
%         for tIndVar=tindClosest+1:length(spike_times)
%             ti_ind=round(spike_times(tIndVar)*Fs);
%             k2(tau1,tau2)=k2(tau1,tau2)+pin(ti_ind-tau1+1)*pin(ti_ind-tau2+1);
%         end
%         k2(tau1,tau2)=k2(tau1,tau2)/2/nRep/(length(pin)/Fs)/std(pin)^4;
%         k2(tau2,tau1)=k2(tau1,tau2);
%         if tau1==tau2
%             k2(tau1,tau2)=k2(tau1,tau2)-k2_diag_correction;
%         end
%     end
% end
% k2a=k2;

%% Else
k2=zeros(window_dur*Fs+1);

for tIndVar=1:length(spike_times)
    ti_ind=round(spike_times(tIndVar)*Fs);
    xcur=pin(ti_ind-1:-1:ti_ind-size(k2,1));
    k2=k2+xcur*xcur';
end
k2=k2/2/nRep/(length(pin)/Fs)/std(pin)^4;
for tau=1:length(k2)
    k2(tau,tau)=k2(tau,tau)-k2_diag_correction;
end
% k2b=k2;



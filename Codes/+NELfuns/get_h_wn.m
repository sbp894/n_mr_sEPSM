function h_an_wn_est=get_h_wn(h_an_ssnavg,fName)

[n_ssn,Fs]=audioread(fName);
h_an_ssnavg=h_an_ssnavg-mean(h_an_ssnavg);

NFFT=2056;
lpcP=30;

fftHssn=abs(fft(h_an_ssnavg, NFFT));
freq=fftshift(linspace(-Fs/2, Fs/2, NFFT));

lpc_coefs=lpc(n_ssn,lpcP);
[lpcFx,lpcFreq]=freqz(1,lpc_coefs,NFFT,'whole',Fs);
lpcFreq=lpcFreq-Fs/2;
lpcFreq=fftshift(lpcFreq);
lpcFx=abs(lpcFx);
lpc_ssn=(lpcFx/max(lpcFx))*max(fftHssn);
fraction_tol=.01;
lpc_ssn(lpc_ssn<fraction_tol*max(lpc_ssn))=fraction_tol*max(lpc_ssn);

lpc_multiplier=interp1(lpcFreq,1./lpc_ssn, freq');
lpc_multiplier(isnan(lpc_multiplier))=1; 

h_an_wn_est=real(ifft(fft(h_an_ssnavg, NFFT).*lpc_multiplier, NFFT));
h_an_wn_est=h_an_wn_est(1:length(h_an_ssnavg));
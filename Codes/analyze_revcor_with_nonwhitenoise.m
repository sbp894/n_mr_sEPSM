function [SNRenv,CF]=analyze_revcor_with_nonwhitenoise(StimsFNames,h,attn)

% disp('check');
NFFT=2056;

stimN_name=StimsFNames{2,1}{1};
[stimN, ~] = audioread(stimN_name);
stimN=stimN*10^(-attn/20);

stimSN_name=StimsFNames{3,1}{1};
[stimSN, Fs] = audioread(stimSN_name);
stimSN=stimSN*10^(-attn/20);

stimN_BP_env=abs(hilbert(real(ifft(fft(stimN).*fft(h,length(stimN))))));
stimSN_BP_env=abs(hilbert(real(ifft(fft(stimSN).*fft(h,length(stimSN))))));

fftXcur=fft(stimN_BP_env,NFFT);
fftXcur=fftXcur(1:NFFT/2+1);
psdN=abs(fftXcur).^2/Fs/length(stimN_BP_env);


fftXcur=fft(stimSN_BP_env,NFFT);
fftXcur=fftXcur(1:NFFT/2+1);
psdSN=abs(fftXcur).^2/Fs/length(stimSN_BP_env);

freqvec=linspace(0,Fs/2,NFFT/2+1)';
ModFreq=[1 2 4 8 16 32 64];

PowerModN=nan(length(ModFreq),1);
PowerModSN=nan(length(ModFreq),1);

for modf_var = 1:length(ModFreq)
    PowerModN(modf_var) = Library.MODENERGY(ModFreq(modf_var), freqvec, psdN);
    PowerModSN(modf_var) = Library.MODENERGY(ModFreq(modf_var), freqvec, psdSN);
end

SNRenv=20*log10(sqrt(nansum(((PowerModSN-PowerModN)./PowerModN).^2)));

fftXcur=fft(h,NFFT);
fftXcur=fftXcur(1:NFFT/2+1);
CF=freqvec(fftXcur==max(fftXcur));
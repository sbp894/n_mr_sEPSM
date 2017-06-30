%%
clear;
close all;
clc;

%%
% load('..\SNRenv-OUTPUT\SumVarAll.mat');
[n_ssn,Fs]=audioread('Z:\Users\SP\MATData\SP-2016_10_28-Q285_AN_Normal\Signals\MH\SNRenv\Stim-6dB_N_P.wav');
n_wn=rms(n_ssn)*randn(size(n_ssn));

NFFT=2056*8;
fft_n_ssn=abs(fft(n_ssn, NFFT));
fft_n_wn=abs(fft(n_wn, NFFT));
freq=linspace(-Fs/2, Fs/2, NFFT);

figure;
plot(freq,fftshift(fft_n_ssn));
hold on;
plot(freq,fftshift(fft_n_wn));
xlim([0 10e3]);
xlabel('Freq');
ylabel('Amp');

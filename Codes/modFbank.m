% This function is an implementation of a modulation filterbank similar to 
% the sEPSM-filterbank as presented by Jørgensen & Dau 2011. This 
% implementation consists of a lowpass filter with a cutoff at 1 Hz, in 
% parallel with 8 bandpass filters with octave spacing.  
%
% Inputs: 
%   binWidth       :  FFT bin width in Hz
%   fcs            :  centerfrequencies of the modulation filters
%   qs             :  Q values associated with the filters centered at fcs
%   plotFilter     :  flag to enable (1)/disable (0) plotting of results
% Outputs:
%   wCf            :  Transferfunction
%
% Created by Christoph Scheidiger 2016
% Copyright Christoph Scheidiger 2016

function [wCf] = modFbank(freqs,fcs,qs,plotFilter)
% 
if nargin<4
    % band center frequencies
    fcs=[1 2 4 8 16 32 64 128 256];
end

% Initialize transfer function
TFs = zeros(length(fcs),length(freqs));
Q=qs;
% Calculating frequency-domain transferfunction for each center frequency:
for k = 2:length(fcs)
    TFs(k,1:end) = 1./(1+ (1j*Q(k)*(freqs(1:end)./fcs(k) - fcs(k)./freqs(1:end)))); % p287 Hambley.
end


% squared filter magnitude transfer functions
wCf = (abs(TFs)).^2;
% cutoff frequency of lowpassfilter:
fcut = fcs(1);
% order:
n = 3;

% Lowpass filter squared transfer function:
wCf(1,:) =  1./(1+((2*pi*freqs/(2*pi*fcut)).^(2*n))); % third order butterworth filter TF from: http://en.wikipedia.org/wiki/Butterworth_filter
TFs(1,:) = sqrt(wCf(1,:));


if plotFilter
    figure
    fnts = 14;
    lw = 2;
    % plot(freqs,10*log10(abs(TFs(1,:))),'linewidth',lw), hold on
    semilogx(freqs,wCf,'linewidth',lw)
    title('Squared transfer functions of the filterbank')
    xlabel('Frequency [Hz]','FontSize',fnts)
    ylabel('Filter attenuation [dB]','FontSize',fnts)
    resultDir = '/SCRATCH/csche/20160404tenReps';
    Library.saveFigureAs([resultDir 'doc/sepsm/filterbank.eps'])
    Library.saveFigureAs([resultDir 'doc/sepsm/filterbank.png'])
end

end

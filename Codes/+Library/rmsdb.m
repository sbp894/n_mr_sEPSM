function r=rmsdb(z)
% returns the rms value of a signal in dB
% 20*log10(norm(z)/sqrt(length(z)))

r=20*log10(norm(z)/sqrt(length(z)));
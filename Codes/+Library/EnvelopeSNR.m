function [SNR_env] = EnvelopeSNR(powerSN, powerN)

if ~isnan(powerSN)
    powerEstimate = max(powerSN - powerN,0);  % Assume (as in Jorgensen and Dau, 2011) that PowSN>=PowN, ow/ set to 0.
else
    powerEstimate = NaN;
end
SNR_env = powerEstimate/powerN;

end


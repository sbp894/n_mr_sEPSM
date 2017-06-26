function [results] = calculateNeuralCorrelograms(SpTimes_A_plus,SpTimes_A_minus,anal)

plt = 0;

if isfield(anal,'SNR')
    nSnr = length(anal.SNR);
else
    nSnr = 1;
end

% Analyze Spike Trains
for snr_i = 1: nSnr
    % Calculate SACs
    [x1,SAC_Ap] = Library.SAC(SpTimes_A_plus(snr_i,:),anal.nTrials,anal.binWidth,anal.duration,anal.maxLag,anal.onsetIgnore,anal.fs);
    [x2,SAC_Am] = Library.SAC(SpTimes_A_minus(snr_i,:),anal.nTrials,anal.binWidth,anal.duration,anal.maxLag,anal.onsetIgnore,anal.fs);
    
    results.x= x1;
    results.SAC_p(snr_i,:) = SAC_Ap;
    results.SAC_m(snr_i,:) = SAC_Am;
    
    % Calculate SCCs
    [x1,SCC_Ap_Am] = Library.SCC(SpTimes_A_plus(snr_i,:),SpTimes_A_minus(snr_i,:),anal.nTrials,anal.binWidth,anal.duration,anal.maxLag,anal.onsetIgnore,anal.fs);
    
    results.SCC_p_m(snr_i,:) = SCC_Ap_Am;
    
    if plt
        plot(x1*1e3,SAC_Ap,'-k','LineWidth',2); hold on;
        plot(x2*1e3,SCC_Ap_Am,'-r');
        title('Speech')
        legend('SAC(A+)','SCC(A+,A-)');
        legend boxoff;
    end
    
    %% DIFCORs emphasizing fine structure
    DIFCOR_A = SAC_Ap - SCC_Ap_Am;
    
    maxDIFCOR_A = DIFCOR_A(x1==0); % same as max(DIFCOR_A), but faster
    
    if plt
        plot(x1*1e3,DIFCOR_A,'-k','LineWidth',2);
        legend('DIFCOR','Location','northeast');
        legend boxoff;
    end
    
    results.DIFCOR(snr_i,:) = DIFCOR_A;
    results.maxDIFCOR(snr_i)=maxDIFCOR_A;
    
    %% SUMCORS emphazizing ENV
    SUMCOR_A = mean([SAC_Ap; SCC_Ap_Am]);
    
    maxSUMCOR_A = SUMCOR_A(x1==0); % same as max(SUMCOR_A), but faster
    
    
    results.SUMCOR(snr_i,:) = SUMCOR_A;
    results.maxSUMCOR(snr_i)=maxSUMCOR_A;
end
end
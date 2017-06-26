function [results] = calculateNeuralCorrelograms(SpTimes_A_plus,SpTimes_A_minus,SpTimes_B_plus,SpTimes_B_minus,anal)

plt = anal.plt;
% Analyze Spike Trains

% Calculate SACs
[x1,SAC_Ap] = Library.SAC(SpTimes_A_plus,anal.nTrials,anal.binWidth,anal.duration,anal.maxLag,anal.onsetIgnore,anal.fs);
[x2,SAC_Am] = Library.SAC(SpTimes_A_minus,anal.nTrials,anal.binWidth,anal.duration,anal.maxLag,anal.onsetIgnore,anal.fs);

[x1,SAC_Bp] = Library.SAC(SpTimes_B_plus,anal.nTrials,anal.binWidth,anal.duration,anal.maxLag,anal.onsetIgnore,anal.fs);
[x2,SAC_Bm] = Library.SAC(SpTimes_B_minus,anal.nTrials,anal.binWidth,anal.duration,anal.maxLag,anal.onsetIgnore,anal.fs);

% Calculate SCCs
[x1,SCC_Ap_Am] = Library.SCC(SpTimes_A_plus,SpTimes_A_minus,anal.nTrials,anal.binWidth,anal.duration,anal.maxLag,anal.onsetIgnore,anal.fs);
[x2,SCC_Ap_Bm] = Library.SCC(SpTimes_A_plus,SpTimes_B_minus,anal.nTrials,anal.binWidth,anal.duration,anal.maxLag,anal.onsetIgnore,anal.fs);
[x1,SCC_Bp_Bm] = Library.SCC(SpTimes_B_plus,SpTimes_B_minus,anal.nTrials,anal.binWidth,anal.duration,anal.maxLag,anal.onsetIgnore,anal.fs);
[x1,SCC_Ap_Bp] = Library.SCC(SpTimes_B_plus,SpTimes_B_minus,anal.nTrials,anal.binWidth,anal.duration,anal.maxLag,anal.onsetIgnore,anal.fs);

if plt
    figure
    subplot(331)
    plot(x1*1e3,SAC_Ap,'-k','LineWidth',2); hold on;
    plot(x2*1e3,SCC_Ap_Am,'-r');
    title('Mod PT C1k M10')
    legend('SAC(A+)','SCC(A+,A-)');
    legend boxoff;
    subplot(332)
    plot(x1*1e3,SAC_Bp,'-k','LineWidth',2); hold on;
    plot(x2*1e3,SCC_Bp_Bm,'-r');
    title('Mod PT C1k M100')
    legend('SAC(B+)','SCC(B+,B-)');
    legend boxoff;
    subplot(333)
    plot(x1*1e3,SCC_Ap_Bp,'-k','LineWidth',2); hold on;
    plot(x2*1e3,SCC_Ap_Bm,'-r');
    title('Speech / Noisy Speech')
    legend('SCC(A+,B+)','SCC(A+,B-)');
    legend boxoff;
end

%% DIFCORs emphasizing fine structure
DIFCOR_A = SAC_Ap - SCC_Ap_Am;
DIFCOR_B = SAC_Bp - SCC_Bp_Bm;
DIFCOR_AB = SAC_Ap - SCC_Ap_Bm;

maxDIFCOR_A = DIFCOR_A(x1==0); % same as max(DIFCOR_A), but faster
maxDIFCOR_B = DIFCOR_B(x1==0);
maxDIFCOR_AB = DIFCOR_AB(x1==0);


% Roh_Tanal.fs
rho_tfs = maxDIFCOR_AB / sqrt( maxDIFCOR_A * maxDIFCOR_B );
if plt
    subplot(334)
    plot(x1*1e3,DIFCOR_A,'-k','LineWidth',2);
    ylabel('DIFCOR')
    %legend('DIFCOR','Location','northeast');
    %legend boxoff;
    subplot(335)
    plot(x1*1e3,DIFCOR_B,'-k','LineWidth',2);
    %legend('DIFCOR');
    %legend boxoff;
    subplot(336)
    plot(x1*1e3,DIFCOR_AB,'-k','LineWidth',2);
    text(min(x1*1e3)+0.25,max(DIFCOR_AB),['\rho_{TFS} = ' num2str(rho_tfs)] );
    %legend('DIFCOR');
    %legend boxoff;
end

%% SUMCORS emphazizing ENV
SUMCOR_A = mean([SAC_Ap; SCC_Ap_Am]);
SUMCOR_B = mean([SAC_Bp; SCC_Bp_Bm]);
SUMCOR_AB = mean([SAC_Ap; SCC_Ap_Bm]);

maxSUMCOR_A = SUMCOR_A(x1==0); % same as max(SUMCOR_A), but faster
maxSUMCOR_B = SUMCOR_B(x1==0);
maxSUMCOR_AB = SUMCOR_AB(x1==0);

% Roh_ENV
rho_ENV = (maxSUMCOR_AB-1) / sqrt( (maxSUMCOR_A-1) * (maxSUMCOR_B-1) );

if plt
    subplot(337)
    plot(x1*1e3,SUMCOR_A,'-k','LineWidth',2);
    %legend('SUMCOR','Location','northeast');
    %legend boxoff;
    ylabel('SUMCOR')
    subplot(338)
    plot(x1*1e3,SUMCOR_B,'-k','LineWidth',2);
    xlabel('Delay [ms]')
    %legend('SUMCOR');
    %legend boxoff;
    subplot(339)
    plot(x1*1e3,SUMCOR_AB,'-k','LineWidth',2);
    text(min(x1*1e3)+0.25,max(SUMCOR_AB),['\rho_{ENV} = ' num2str(rho_ENV)] );
    %legend('SUMCOR');
    %legend boxoff;
    
    
    Library.mtit(['Neural Correlations CF = ' num2str(anal.CF) 'Hz'],'yoff',0.02,'xoff',0)
    Library.saveFigureAs( [anal.name 'SacScc_CF' num2str(anal.CF) '.eps'] )
end

results.maxSUMCOR_A=maxSUMCOR_A;
results.maxSUMCOR_B=maxSUMCOR_B;
results.maxSUMCOR_AB=maxSUMCOR_AB;

end
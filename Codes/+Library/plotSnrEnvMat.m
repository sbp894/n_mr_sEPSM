function [ output_args ] = plotSnrEnvMat(  )
%Plots the SNR env matrices from the mat files passed in the argument
%   Plots the SNR env matrices from the mat files passed in the argument
    
    load('snrEnvs_Varsha_10sent_04-11-2015.mat')
    
    minval = min(result.avgMat(:));
    maxval = max(result.avgMat(:));
    
    for snr_i = 1 : length(result.SNRs)
        figure
        imagesc(squeeze(result.avgMat(snr_i,:,:)),[minval,maxval])
        title(['SNRenv matrix, ' num2str(result.SNRs(snr_i)) 'dB SNR'],'FontSize',18)
        set(gca,'YTick',1:length(result.pF))
        set(gca,'YTickLabel',result.pF)
        ylabel('CF')
        set(gca,'XTick',1:length(result.mF))
        set(gca,'XTickLabel',result.mF)
        xlabel('Modulation filter (f_c)')
        h = colorbar;
        ylabel(h, 'SNRenv')
        saveFigureAs(['./doc/snrEnvMatrix_' num2str(result.SNRs(snr_i)) 'dBSnr.png'])
    end
    close all;
end


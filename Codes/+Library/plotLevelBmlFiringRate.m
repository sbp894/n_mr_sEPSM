function plotLevelBmlFiringRate(filename)
    load(filename);
    % bmlData.levels = levels;
    % bmlData.CF = AN.CF;
    % bmlData.avgFR_A = avgFR_A;
    % bmlData.sumcorPeak_A = sumcorPeak_A;
    figure
    xMax = bmlData.levels(end);
    xMin = bmlData.levels(1);
    y1Max = max(bmlData.avgFR_A(:));
    y1Min = min(bmlData.avgFR_A(:));
    y2Max = max(bmlData.sumcorPeak_A(:));
    y2Min = min(bmlData.sumcorPeak_A(:));
    subplot(length( bmlData.CF ), 1,1);
    
   % ha = Library.tight_subplot(length( bmlData.CF ),1,[.01 .03],[.1 .01],[.01 .01])

    for cF_i = 1: length( bmlData.CF )
        %axes(ha(cF_i));
        subplot(length( bmlData.CF ), 1,cF_i);
        [ax,h1,h1]=plotyy( bmlData.levels, bmlData.avgFR_A(cF_i,:), bmlData.levels, bmlData.sumcorPeak_A(cF_i,:));
        set(ax(1),'YLim',[y1Min y1Max])
        set(ax(2),'YLim',[y2Min y2Max])
        set(ax(1),'XLim',[xMin xMax])
        set(ax(2),'XLim',[xMin xMax])
        ylabel([ num2str(bmlData.CF(cF_i)) ' Hz' ])
        if cF_i == length( bmlData.CF )
            xlabel([ 'Level dB' ])
            legend('Firing Rate (sp/s)', 'max SC')
        end
    end
    Library.saveFigureAs([bmlData.name '.eps'])
end
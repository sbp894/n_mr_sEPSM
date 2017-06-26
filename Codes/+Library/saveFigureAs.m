function saveFigureAs( saveAs )
% saves a figure as graphic file with the extensions mentioned in the
% saveAs variable and in the path
if ~isempty(saveAs)
    ext = saveAs(end-2:end);
    switch 1
        case strcmp(ext,'png')
            printDriver =  '-dpng';
        case strcmp(ext,'eps')
            printDriver = '-depsc2';
        case strcmp(ext,'jpg')
            printDriver = '-djpg';
        case strcmp(ext,'pdf')
            printDriver = '-dpdf';
        case strcmp(ext,'iff') % nothing else will end with iff
            printDriver = '-dtiff';
        otherwise
            return
    end
    print(gcf, saveAs, printDriver);
end
end

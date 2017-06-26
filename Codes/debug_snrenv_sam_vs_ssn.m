close all;
lw=2;
for ntypeVar=1:length(noisetypes)
    for cfVar=1:3:length(cfs)
        figure;
        for levelVar=1:length(levels)
            for snrVar=1:length(snrs)
                subplot(4,3,snrVar);
                for sentVar=1:length(sents)
                    
                    errorbar(mean(SimData(ntypeVar, levelVar, snrVar, sentVar, cfVar).SNRenvbrick{1}), std(SimData(ntypeVar, levelVar, snrVar, sentVar, cfVar).SNRenvbrick{1}), '-d','linewidth',lw);
                    hold on;
                    errorbar(mean(SimData(ntypeVar, levelVar, snrVar, sentVar, cfVar).SNRenvbrick{2}), std(SimData(ntypeVar, levelVar, snrVar, sentVar, cfVar).SNRenvbrick{1}), '-d','linewidth',lw);
                    errorbar(mean(SimData(ntypeVar, levelVar, snrVar, sentVar, cfVar).SNRenvbrick{3}), std(SimData(ntypeVar, levelVar, snrVar, sentVar, cfVar).SNRenvbrick{1}), '-d','linewidth',lw);
                    plot(SimData(ntypeVar, levelVar, snrVar, sentVar, cfVar).SNRenv, 'o--','linewidth',lw);
                    hold off;
                    
                    snrBrick=nan; %nansum(mean(SimData(ntypeVar, levelVar, snrVar, sentVar, cfVar).SNRenvbrick).^2);
                    snrFinal=nan; %nansum(SimData(ntypeVar, levelVar, snrVar, sentVar, cfVar).SNRenv.^2);
                    title(sprintf('%s, snr=%i, cf=%1.2f', ... %, snr-brick=%1.2f, snr-from_finalpsd=%1.2f',...
                        noisetypes{ntypeVar}, snrs(snrVar), cfs(cfVar))); %, snrBrick, snrFinal));
                    ylabel('SNRenv'); xlabel('modfreq');
                    set(gca, 'XTick', 1:7, 'xticklabel',  strread(num2str(ModFreqs),'%s'));
                    ylimold=ylim;
                    ylim([0 ylimold(2)]);
                end
            end
        end
        legend('brick-LSR','brick-MSR','brick-HSR', 'final PSD');
    end
end
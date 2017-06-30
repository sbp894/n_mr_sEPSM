function [pModSNOut, pModNOut, SNRenvOut]=mrSNRenvAnalysis(PowerModCell, CompareWithNF, NumSTDsNFtol)

pModSNOut=nan(1,length(PowerModCell));
pModNOut=nan(1,length(PowerModCell));
SNRenvOut=nan(1,length(PowerModCell));

for modVar=1:length(PowerModCell)
    
    modCells=PowerModCell{modVar};
    tempSNmodP=nan(1,length(modCells));
    tempNmodP=nan(1,length(modCells));
    tempSNRenv=nan(1,length(modCells));
    
    for windVar=1:length(modCells)
        
        if ~isempty(modCells{windVar})
            PowerModSN=modCells{windVar}.PowerModSN(:,modVar);
            pModSN=mean(PowerModSN);
            if CompareWithNF
                PowerModSN_noisefloor=modCells{windVar}.PowerModSN_noisefloor(:,modVar);
                zSN=(mean(PowerModSN)-mean(PowerModSN_noisefloor))./(sqrt(std(PowerModSN).*std(PowerModSN_noisefloor)));
                pModSN(zSN<NumSTDsNFtol)=nan;
            end
            
            PowerModN=modCells{windVar}.PowerModN(:,modVar);
            pModN=mean(PowerModN);
            if CompareWithNF
                PowerModN_noisefloor=modCells{windVar}.PowerModN_noisefloor(:,modVar);
                u_nNF=mean(PowerModN_noisefloor);
                std_nNF=std(PowerModN_noisefloor);
                zN=(mean(PowerModN)-u_nNF)./sqrt(std(PowerModN).*std_nNF);
                %pModN(nanInds)=max([pModN(nanInds); u_nNF(nanInds)])+NumSTDsNFtol*std_nNF(nanInds);
                if zN<NumSTDsNFtol
                    pModN=u_nNF+NumSTDsNFtol*std_nNF;
                end
            end
            
            SNRenv=(pModSN-pModN)./pModN;
            SNRenv(SNRenv<0)=nan;
            
            tempSNmodP(windVar)=pModSN;
            tempNmodP(windVar)=pModN;
            tempSNRenv(windVar)=SNRenv;
        end
    end
    
    pModSNOut(modVar)=nansum(tempSNmodP)/length(modCells);
    pModNOut(modVar)=nansum(tempNmodP)/length(modCells);
    SNRenvOut(modVar)=nansum(tempSNRenv)/length(modCells);
end
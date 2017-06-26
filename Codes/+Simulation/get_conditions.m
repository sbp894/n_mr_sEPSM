function cs=get_conditions(A,B,AN,resultsDir)

%% create condition table
nCF = length(AN.CF);
nSentences = A.numberOfSentences;
nSnr = length(B.SNR);
nNoise = length(B.noiseTypes);
nLevel = length(A.level);
nFibertypes=length(AN.fiberType.sr);


cs = Library.cartesianProduct(nCF,nSentences,nLevel,nNoise,nSnr,nFibertypes);
csCell = cell(size(cs));
for i = 1 : size(cs,1)
    csCell{i,1} = AN.CF(cs(i,1));
    csCell{i,2} = A.sentences(cs(i,2));
    csCell{i,3} = A.level;
    csCell{i,4} = B.noiseTypes{cs(i,4)};
    csCell{i,5} = B.SNR(cs(i,5));
    csCell{i,6} = AN.fiberType.sr(cs(i,6));
end

csLabels = {'CF','sentences','snr','noise','level', 'fibertypes'}; %#ok<*NASGU>
nConditions = [nCF, nSentences, nSnr, nNoise, nLevel, nFibertypes];
save([resultsDir 'conditions.mat'],'cs','csCell','csLabels','nConditions');

function [NSCC,NSCCdelays_usec,AVGrates,TOTALspikes] = ShufCrossCorr(SpikeTrains,DELAYbinwidth,Duration)
% File: ShufCrossCorr
% 27Sep2004: M Heinz - updated from Recruitment - NoiseCorr Analysis
%
% Created: 9/3/03 M. Heinz
% Calls: InnerSACmex.dll  (MEX-C file: InnerSACmex.c)
%
% Computes Normalized Shuffled Cross-Correlogram (NSCC) from a set of Spike Trains and Duration

for UNITind=1:2
    NUMspikeREPS{UNITind}=length(SpikeTrains{UNITind});
end

%%%%%%%%%%%%% Setup BIG spike matrix for faster computation
for UNITind=1:2
    Kmax{UNITind}=0;
    for spikeREPind=1:NUMspikeREPS{UNITind}
        if length(SpikeTrains{UNITind}{spikeREPind})>Kmax{UNITind}
            Kmax{UNITind}=length(SpikeTrains{UNITind}{spikeREPind});
        end
    end
end

%%%% Compute AVGrates
for UNITind=1:2
    SpikeMAT{UNITind}=NaN*ones(NUMspikeREPS{UNITind},Kmax{UNITind});
    for REPindREF=1:NUMspikeREPS{UNITind}
        SpikeMAT{UNITind}(REPindREF,1:length(SpikeTrains{UNITind}{REPindREF}))=SpikeTrains{UNITind}{REPindREF};
    end
    NUMspikes{UNITind}=sum(~isnan(SpikeMAT{UNITind}(:,:))',1)';  % Count number of real spikes in each line
    TOTALspikes{UNITind}=sum(NUMspikes{UNITind});
    AVGrates{UNITind}=TOTALspikes{UNITind}/NUMspikeREPS{UNITind}/Duration;
end


[intsMEX,TOTALints] = Library.InnerSCCmex(SpikeMAT{1}',NUMspikes{1},TOTALspikes{1},SpikeMAT{2}',NUMspikes{2},TOTALspikes{2});
% save('haha_sc_gotcha.mat', 'SpikeTrains' ,'DELAYbinwidth', 'Duration', 'SpikeMAT', 'NUMspikes', 'TOTALspikes', 'intsMEX', 'TOTALints');
% fprintf('%i/%i\n', length(intsMEX), TOTALints);
ints=intsMEX(1:min(TOTALints,length(intsMEX)));  % Remove extra ints due to NaN's in SpikeMAT1 matrix
if TOTALints>length(intsMEX)
    warning('something is wrong in ShufCrossCorr\n');
end
clear intsMEX
ints=ints(~isnan(ints))/1e-6;  % Convert into micro-seconds

%%% Normalize SAC such that no temporal correlation = 1
NSCCdelays_usec=0:DELAYbinwidth:Duration;
NSCCdelays_usec=[-fliplr(NSCCdelays_usec(2:end)) NSCCdelays_usec]/1e-6;  % Convert into micro-seconds
SCC=hist(ints,NSCCdelays_usec);
% From Louage et al (2004: J. Neurophysiol) [This is normalized to take out rate effects, such that a NSCC=1 means
% no more coincident than a stationary Poisson Process, above 1 means the two inputs tend to be more coincident than expected]
NSCC=SCC/(NUMspikeREPS{1}*NUMspikeREPS{2}*Duration*AVGrates{1}*AVGrates{2}*DELAYbinwidth);


return;
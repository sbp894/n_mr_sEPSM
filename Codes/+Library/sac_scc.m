function [counts, binCenters] = sac_scc(tSpikeTrains, D, maxLag, w, normalization, spikesToChoose)
% sac_scc - Compute the shuffled autocorrelogram (SAC) or shuffled crosscorrelogram (SAC).
%
% [COUNTS, BINCENTERS] = sac_scc(TSPIKETRAINS, D, MAXLAG, W, NORMALIZATION)
%   Calculate the shuffled autocorrelogram (SAC) or shuffled
%   crosscorrelogram (SCC) of given spike trains TSPIKETRAIN.
%   TSPIKETRAIN is a cell array. Each element corresponds to a spike
%       train, which could be the response to the same stimuli or to a
%       different one. Each column corresponds to a response itself. Each
%       row corresponds to a channel. One or two channels are supported.
%       In the case of a single channel, SAC is calculated. In case of two
%       channels, SCC is calculated.
%   D defines the duration of the whole response.
%   MAXLAG defines the limits of the SAC/SCC, which goes from
%       -MAXLAG to +MAXLAG. However, it is worth noting that the SAC/SCC
%       is symmetric, so we'll only do the computations from
%       0 to +MAXLAG.
%   W is the bin width (coincidence window). A value of 50 us is
%       recommended.
%   NORMALIZATION is a string that defines the type of normalization
%       to be performed on the SAC/SCC. Possible values are:
%       'none'      Perform no normalization at all.
%       'ci'        Perform correlation index normalization.
%       'cr'        Perform coincidence rate normalization.
%   SPIKESTOCHOOSE is a string that defines what spikes to use for the
%       interval calculation (considering the actual peak as a reference).
%       Possible values are:
%       'forward'	Consider only time forward spikes.
%                   This case yields a one-sided (positive) histogram.
%       'backward'  Consider only time backward spikes.
%                   This case yields a one-sided (negative) histogram.
%       'both'      Consider both time forward and time backward spikes.
%                   This case yields a two-sided histogram.
%
%   The output COUNTS contains the corresponding histogram counts.
%   BINCENTERS gives the value of the bin centers. The bin in the middle
%   refers to 0. To generate a plot, something like this can be used:
%
%       bar(binCenters, counts);
%
%   For a complete explanation of SAC/SCC, refer to the following papers:
%
%       Joris, Philip X., et al. "Correlation index: a new metric to
%       quantify temporal coding." Hearing research 216 (2006): 19-30.
%
%       Louage, Dries HG, Marcel van der Heijden, and Philip X. Joris.
%       "Temporal properties of responses to broadband noise in the
%       auditory nerve." Journal of neurophysiology 91.5 (2004): 2051-2065.
%
%       Joris, Philip X., et al. "Dependence of binaural and cochlear
%       "best delays" on characteristic frequency." Auditory Signal
%       Processing. Springer New York, 2005. 477-483.
%
% Created February 24, 2015.
% Arturo Moncada-Torres
%   arturo.moncadatorres@med.kuleuven.be
%   http://www.arturomoncadatorres.com
%
% Updated March 11, 2015.
%   -> Added the option to choose what spikes to use: forward spikes,
%       backward spikes, or both.
% Arturo Moncada-Torres
%   arturo.moncadatorres@med.kuleuven.be
%   http://www.arturomoncadatorres.com
%
% Updated March 13, 2015.
%   -> Corrected bug when nChannels == 2.
% Arturo Moncada-Torres
%   arturo.moncadatorres@med.kuleuven.be
%   http://www.arturomoncadatorres.com


%% Manage inputs.
nChannels = size(tSpikeTrains,1);
switch nChannels
    case {1,2}
        % Do nothing. This is a valid case.
    otherwise
        error([mfilename, ':inputs'],'Invalid number of channels.');
end


%% Calculate important parameters.

% Number of presentations.
for ii = 1:nChannels
    M(ii) = length(tSpikeTrains(ii,:));
end

% Bin centers (including negative ones).
binCenters = -maxLag+w : w : maxLag;


%% Make sure spike trains are in ascending order.
for ii = 1:nChannels
    currM = M(ii);
    for jj = 1:currM
        tSpikeTrains{ii,jj} = sort(tSpikeTrains{ii,jj});
    end
end


%% Calculate actual metric.

for jj = 1:min(M)   % Use min to avoid any possible overflow.
    
    switch nChannels
        case 1
            % Extract current spike train.
            currSpikeTrain = tSpikeTrains{jj};
            
            % Remove current spike train from the rest.
            tSpikeTrains_temp = tSpikeTrains;
            tSpikeTrains_temp(jj) = [];
            
            % Unfold remaining spike trains into a single array.
            tSpikeTrains_temp = [tSpikeTrains_temp{:}];
            
            
        case 2
            % Extract current spike train.
            currSpikeTrain = tSpikeTrains{1,jj};
            
            % Unfold spike trains from second channel into a single array.
            tSpikeTrains_temp = [tSpikeTrains{2,:}];
            
    end % switch nChannels
    
    
    % Go through all the spikes of the current spike train.
    nSpikes = length(currSpikeTrain);
    
    timeIntervals = cell(nSpikes,1);
    for kk = 1:nSpikes
        currSpike = currSpikeTrain(kk);
        
        % Select corresponding time intervals.
        switch lower(spikesToChoose)
            case 'forward'
                tSpikeTrains_chosen = tSpikeTrains_temp(tSpikeTrains_temp >= currSpike);
            case 'backward'
                tSpikeTrains_chosen = tSpikeTrains_temp(tSpikeTrains_temp <= currSpike);
            case 'both'
                tSpikeTrains_chosen = tSpikeTrains_temp;
            otherwise
                error([mfilename, ':inputs'],'Invalid spikesToChoose parameter.');
        end
        
        % Measure the time from the current spike to all forward time
        % spikes of the remaining spike train.
        timeIntervals{kk} = tSpikeTrains_chosen - currSpike;
    end
    timeIntervals = [timeIntervals{:}];
    
    
    % Perform binning.
    if ~exist('counts','var')
        counts = hist(timeIntervals, binCenters);
    else
        counts = counts + hist(timeIntervals, binCenters);
    end

    if mod(jj,25) == 0
        fprintf('\n\t%d of %d repetitions...', jj,currM);
    end
end % for currM
fprintf('\n');


%% Perform normalization.

switch lower(normalization)
    case 'none'
        % Do nothing.
        normFactor = 1;
        
        
    case 'ci'
        % Calculate correlation index normalization factor, as explained
        % in joris2006correlation (p. 22).

        % Calculate average firing rate.
        for ii = 1:nChannels
            for jj = 1:M
                r(ii,jj) = numel(tSpikeTrains{jj}) ./ D;
            end
        end
        r = mean(r,2);

        % Actual normalization factor.
        switch nChannels
            case 1
                normFactor = M .* (M-1) .* (r^2) .* w .* D;
            case 2
                normFactor = prod(M) .* prod(r) .* w .* D;
        end

        
    case 'cr'
        % Calculate coincidence rate normalization, as explained
        % in joris2006correlation (p. 22).

        % Actual normalization factor.
        switch nChannels
            case 1
                normFactor = M * (M-1) * w * D;
            case 2
                normFactor = prod(M) * w * D;
        end
        
 
    otherwise
        error([mfilename, ':inputs'],'Invalid normalization option.');
end

% Actually perform the normalization.
counts = counts / normFactor;
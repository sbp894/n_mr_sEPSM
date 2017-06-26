function [SpTimes_A_Plus, SpTimes_A_Minus, duration] = getSpikeTrains(A, stim_i, AN, cF_i,snr_i,noise_i,level_i)

level = A.level(level_i);
fileName = [A.prefix sprintf('%02d', stim_i) A.fileExtension];
wav = Library.getWavFromPath(A.path,fileName, A.fs);
duration = length(wav)*1/A.fs;

if isfield(A, 'SNR')
    nSnr = length(A.SNR);
    noisySpeech = wav;
    % calculate rms of all sentences
    allSentences = [];
    allNoise = [];
    for sentence_i = 1 : A.numberOfSentences
        fileName = [A.quietPrefix sprintf('%02d', sentence_i) A.fileExtension];
        speech = Library.getWavFromPath(A.path,fileName, A.fs);
        allSentences = [allSentences speech'];
        fileName = [A.noisePrefix sprintf('%02d', sentence_i) A.fileExtension];
        noise = Library.getWavFromPath(A.path,fileName, A.fs);
        allNoise = [allNoise noise'];
    end
    speechMatRms = Library.rms(allSentences);
    sentenceFileLevel = 10*log10(speechMatRms);
    noiseRms = Library.rms(allNoise);
    % get stimulus to calculate spike trains
    fileName = [A.quietPrefix sprintf('%02d', stim_i) A.fileExtension];
    speech = Library.getWavFromPath(A.path,fileName, A.fs);
    fileName = [A.noisePrefix sprintf('%02d', stim_i) A.fileExtension];
    noise = Library.getWavFromPath(A.path,fileName, A.fs);

    % Adjust SNR
    speechRms = Library.rms(speech);
    noiseRms = Library.rms(noise);
    noiseRmsDesired = speechRms/( 10^( A.SNR(snr_i)/20 ) );
    noise = noise .* ( noiseRmsDesired / noiseRms);
    % noiseRmsAdjusted = Library.rms(noise);
    
    % Mix signals
    wav = noise + speech;
    % rmsNoisySpeech = Library.rms(wav);
    % A.SNR(snr_i)
    % p=audioplayer(wav,A.fs);
    % play(p);
    % p=audioplayer(noisySpeech,A.fs);
    % play(p);
    
    % Adjust level
    dBSPL=20*log10(sqrt(mean(wav.^2))/(20e-6));
    wavLevelAdj=wav*10^((level-dBSPL)/20);
    % p=audioplayer(wavLevelAdj,A.fs);
    % play(p);
else
    dBSPL=20*log10(sqrt(mean(wav.^2))/(20e-6));
    wavLevelAdj=wav*10^((level-dBSPL)/20);
end

AN_duration = duration*1; % AN model is run for 1 times the duration of signal
timeVec = 0:1/A.fs:AN_duration*2-1/A.fs; % make time array for AN spike times

    
% A_plus
pin = wavLevelAdj;
%SpTimes_A_Plus = {};
%%%%%% Run peripheral model
vihc = ANModel.model_IHC(pin',AN.CF(cF_i),1,1/A.fs,duration+0.0001,AN.cohc,AN.cihc,AN.species);
for iTrial = 1:AN.nTrials
    %%%%%% Run the auditory nerve model
    [~]=evalc('[~,~,psth] = ANModel.model_Synapse(vihc,AN.CF(cF_i),1,1/A.fs,AN.fiberType.sr,AN.noiseType,AN.implnt)');
    idx = psth==1;
    data = timeVec(idx);
    %%%%%%%%%%% Save data %%%%%%%%%%%
    SpTimes_A_Plus{iTrial} = data;
end

% A_minus
pin = -wavLevelAdj;
%SpTimes_A_Minus = {};
%%%%%% Run peripheral model
vihc = ANModel.model_IHC(pin',AN.CF(cF_i),1,1/A.fs,duration+0.0001,AN.cohc,AN.cihc,AN.species);
for iTrial = 1:AN.nTrials
    %%%%%% Run the auditory nerve model
    [~]=evalc('[~,~,psth] = ANModel.model_Synapse(vihc,AN.CF(cF_i),1,1/A.fs,AN.fiberType.sr,AN.noiseType,AN.implnt)');
    idx = psth==1;
    data = timeVec(idx);
    %%%%%%%%%%% Save data %%%%%%%%%%%
    SpTimes_A_Minus{iTrial} = data;
end



end
function [A,B]=get_speech_params(Fs, ExpControlParams)

if isempty(ExpControlParams)
    A.level = 50; %dB SPL according to email from Varsha
    B.SNR = [-6 0 6];
    B.noiseTypes = {'SSN','SAM'};
    B.noisePrefix = {['noise' filesep 'ssm_simulation_dtu'], ['noise' filesep 'sam_simulation_dtu']};    %     B.noiseTypes = {'SAM'};
    A.sentences=3;
else
    A.level=ExpControlParams.level;
    B.SNR=ExpControlParams.SNR;
    B.noiseTypes=ExpControlParams.noiseTypes;
    B.noisePrefix = ExpControlParams.noisePrefix;
    A.sentences=ExpControlParams.sentences;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%% Stimulus A (speech) %%%%%%%%%%%%%%%%%%%%%
A.path = ['stimuli' filesep];
A.fileExtension = '.wav';
A.prefix = ['quiet' filesep 'dan_sent'];
A.numberOfSentences = length(A.sentences);
A.Fs = Fs;


%% %%%%%%%%%%%%%%%%%%%%%%%%% Stimulus B (noise) %%%%%%%%%%%%%%%%%%%%%

if strcmp(B.noiseTypes{1}, 'SSN')
    A.StimDir=[pwd filesep 'stimuli' filesep 'stimSetStationary' filesep];
%     B.noisePrefix = {['noise' filesep 'ssn_simulation_dtu']};
elseif strcmp(B.noiseTypes{1}, 'SAM')
    A.StimDir=[pwd filesep 'stimuli' filesep 'stimSetFluctuating' filesep];
%     B.noisePrefix = {['noise' filesep 'sam_simulation_dtu']};
end

% B.noiseTypes = {'SSN','SAM'};
B.path = ['stimuli' filesep];
B.fileExtension = '.wav';
% B.noisePrefix = {['noise' filesep 'SN_Varsha_Sentences'];...
%     ['noise' filesep 'SN_Varsha_Sentences_SAM']};
B.prefix = ['speechInNoise' filesep 'speechInNoise'];
B.numberOfSentences = A.numberOfSentences;
B.level = A.level; %dB SPL according to email from Varsha
B.fs = Fs;

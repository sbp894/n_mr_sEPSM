function [A,B]=get_speech_params

%%%%%%%%%%%%%%%%%%%%%%%%% Stimulus A (speech) %%%%%%%%%%%%%%%%%%%%%
A.path = ['stimuli' filesep];
A.fileExtension = '.wav';
A.prefix = ['quiet' filesep 'quiet'];
A.numberOfSentences = 1;
A.level = 50; %dB SPL according to email from Varsha
% A.fs = Fs;

%%%%%%%%%%%%%%%%%%%%%%%%% Stimulus B (noise) %%%%%%%%%%%%%%%%%%%%%
B.noiseTypes = {'SSN','SAM'};
B.path = ['stimuli' filesep];
B.fileExtension = '.wav';
B.noisePrefix = {['noise' filesep 'SN_Varsha_Sentences'];...
    ['noise' filesep 'SN_Varsha_Sentences_SAM']};
B.prefix = ['speechInNoise' filesep 'speechInNoise'];
B.numberOfSentences = A.numberOfSentences;
B.level = A.level; %dB SPL according to email from Varsha
% B.fs = Fs;
% B.SNR = [-6 -3 0 3 6];
function update_progress(progressDir,resultPostfix2,A,B,C,anal,MaxIter)

% update progress dir
a = 0==0; % dummy variable to save in progress mat file
Library.parsave([progressDir filesep resultPostfix2 '.mat'],a);
D = dir([progressDir filesep '*.mat']);
progress = length(D(not([D.isdir])));
resultPostfix1 = sprintf(anal.resultTxt,C.cF_i,C.sentence_i, C.snr_i, B.noiseTypes{C.noise_i}, A.level(C.level_i) );
disp([datestr(now,'HH:MM:SS') ' Progress: ' num2str(progress) '/' num2str(MaxIter) ' (' resultPostfix1 ')']);


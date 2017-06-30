function update_progress(resultsDir,resultPostfix,MaxIter)

% update progress dir
a = 0==0; % dummy variable to save in progress mat file
Library.parsave([resultsDir 'progress' filesep resultPostfix '.mat'],a);
D = dir([resultsDir 'progress' filesep '*.mat']);
progress = length(D(not([D.isdir])));

disp([datestr(now,'HH:MM:SS') ' Progress: ' num2str(progress) '/' num2str(MaxIter) '--' resultPostfix]);


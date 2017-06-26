function fName=show_current_parloop(resultsDir ,resultPostfix)

a = 0==0; % dummy variable to save in progress mat file
updateDir=[resultsDir 'current' filesep];

fName=[updateDir resultPostfix '.mat'];
Library.parsave(fName,a);
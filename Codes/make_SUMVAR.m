%%
clear;
close all;
clc;

%%
ChinIDs=[277, 282, 284, 285];
SumVar=[];

DataDir=[fileparts(pwd) '\OUTPUT\DataAnal\'];
runIndSummary=[0 0 0 0]+1;


for id_var=ChinIDs
    if runIndSummary(ChinIDs==id_var)
        create_summary(id_var);
    end
    allfiles=dir([DataDir '*' num2str(id_var) '*']);
    a=load([DataDir allfiles(end).name '\SummaryVar.mat']);
    SumVar = [SumVar; a.SummaryVar];
end

save([fileparts(pwd) '\OUTPUT\SumVar.mat'],'SumVar');
%% RUN THIS ONCE
% store data in a big matrix
TaskDep.PM=nan(7,5); % 7 days, 5 conditions comparison: 1 task, 2 tasks, 3 tasks,...
TaskDep.M1=nan(7,5);
TaskDep.S1=nan(7,5);
%% RUN THIS FOR THE 5 CONDITIONS
clearvars -except TaskDep
close all
clc

condition='1task/'; % for file directory
conditionNbre=1;

monkey = 'Natalya';

DataDir=(['/Users/Borgo/Google Drive/BSI LOCOMOTION Project/Research Material/Analysis/PROCESSED DATA/' monkey '/dPCAmatricesCtrlTOM/']);

Files=dir([DataDir condition '*.mat']);
%  old stuff for mean across condition and not across days
% if strcmp(condition,'1task/')
%     load([DataDir condition 'varMat1task.mat'])
%     for n = 1:numel(varMat)
%         TaskDep.PM(n,conditionNbre) = varMat(n).PMdep;
%         TaskDep.M1(n,conditionNbre) = varMat(n).M1dep;
%         TaskDep.S1(n,conditionNbre) = varMat(n).S1dep;
%     end
% 
% else
%     for n=1:numel(Files)
%         load([DataDir condition Files(n).name])
%         for day = 1:numel(varMat)
%             tempPM(day)=varMat(day).PMdep;
%             tempM1(day)=varMat(day).M1dep;
%             tempS1(day)=varMat(day).S1dep;
%         end
%             TaskDep.PM(day,conditionNbre)=mean(tempPM);
%             TaskDep.M1(day,conditionNbre)=mean(tempM1);
%             TaskDep.S1(n,conditionNbre)=mean(tempS1);
%     end  
% end


for n=1:numel(Files)
    load([DataDir condition Files(n).name])
    
    for day = 1:numel(varMat)
        tempPM(day,n)=varMat(day).PMdep;
        tempM1(day,n)=varMat(day).M1dep;
        tempS1(day,n)=varMat(day).S1dep;
    end
    
end

TaskDep.PM(:,conditionNbre)=mean(tempPM,2);
TaskDep.M1(:,conditionNbre)=mean(tempM1,2);
TaskDep.S1(:,conditionNbre)=mean(tempS1,2);
%% ONCE THE 5 CONDITIONS RAN, RUN THIS FOR THE PLOT
%Colormap of brain
cmapbrain=[158,55,123;... %[82,163,204]/255; % S1
    238,98,109;... %[46,115,69]/255; % M1
    61,184,224]; %[168,126,229]/255; % PM
cmapbrain=cmapbrain/255;

figure
for condition = 1:5
    plot(condition,mean(TaskDep.PM(:,condition),'omitnan'),'.','MarkerEdgeColor',cmapbrain(3,:),'MarkerFaceColor',cmapbrain(3,:),'MarkerSize',20);
    hold on
    
    plot(condition,mean(TaskDep.M1(:,condition),'omitnan'),'.','MarkerEdgeColor',cmapbrain(2,:),'MarkerFaceColor',cmapbrain(2,:),'MarkerSize',20);
    plot(condition,mean(TaskDep.S1(:,condition),'omitnan'),'.','MarkerEdgeColor',cmapbrain(1,:),'MarkerFaceColor',cmapbrain(1,:),'MarkerSize',20);
   
    errorbar(condition,mean(TaskDep.PM(:,condition),'omitnan'),std(TaskDep.PM(:,condition),'omitnan'),'Color',cmapbrain(3,:))
    errorbar(condition,mean(TaskDep.M1(:,condition),'omitnan'),std(TaskDep.M1(:,condition),'omitnan'),'Color',cmapbrain(2,:))
    errorbar(condition,mean(TaskDep.S1(:,condition),'omitnan'),std(TaskDep.S1(:,condition),'omitnan'),'Color',cmapbrain(1,:))
  
    ylim([0 100])
    xlim([0 6])
end
hold on
plot(mean(TaskDep.PM,1),'Color',cmapbrain(3,:))
plot(mean(TaskDep.M1,1),'Color',cmapbrain(2,:))
plot(mean(TaskDep.S1,1),'Color',cmapbrain(1,:))


title('task-dependent variance with the introduction of more tasks.')
ylabel('neural VAF')
xticks([1 2 3 4 5])
xticklabels({'1 task','2 tasks','3 tasks','4 tasks','5 tasks'})



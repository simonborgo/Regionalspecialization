%% AMAZING CODE FOR 1 TASK CONDITION

clear all
close all
clc

condition='1task/'; % for file directory
conditionNbre=['2';'3';'4';'5'];

monkey = {'Elektra','Natalya'};

for mk=1:2
    %DataDir=(['/Users/Borgo/Google Drive/BSI LOCOMOTION Project/Research Material/Analysis/PROCESSED DATA/' monkey '/dPCAmatricesCtrlTOM/']);
    % store data in a big matrix
    TaskDep.PM=nan(7,4); % 7 days, 4 conditions comparison: 1 task splitted into: 2, 3, 4 and 5
    TaskDep.M1=nan(7,4);
    TaskDep.S1=nan(7,4);

    DataDir=(['/Users/Borgo/Documents/MATLAB/fri_mks/BSI LOCOMOTION Project_end_of_july_2021/Research Material/Analysis/PROCESSED DATA/' monkey{mk} '/dPCAmatricesCtrlTOM/']);
    
    Files=dir([DataDir condition '*.mat']);
    
    
    for cond=1:length(conditionNbre) %loop across the 4 conditions
        
        for n=1:numel(Files) %loop across the varMat files
            %get the correct condition
            idxCond=find(Files(n).name == conditionNbre(cond));
            if ~isempty(idxCond)
                load([DataDir condition Files(n).name])
                
                for day = 1:numel(varMat)
                    tempPM(day,n)=varMat(day).PMdep;
                    tempM1(day,n)=varMat(day).M1dep;
                    tempS1(day,n)=varMat(day).S1dep;
                end
                
            end
        end
        TaskDep.PM(:,cond)=mean(tempPM,2);
        TaskDep.M1(:,cond)=mean(tempM1,2);
        TaskDep.S1(:,cond)=mean(tempS1,2);
        clear tempS1 tempM1 tempPM
    end
    
    %Colormap of brain
    cmapbrain=[158,55,123;... %[82,163,204]/255; % S1
        238,98,109;... %[46,115,69]/255; % M1
        61,184,224]; %[168,126,229]/255; % PM
    cmapbrain=cmapbrain/255;
    cmapPM=[cmapbrain(3,:);cmapbrain(3,:);cmapbrain(3,:);cmapbrain(3,:)];
    cmapM1=[cmapbrain(2,:);cmapbrain(2,:);cmapbrain(2,:);cmapbrain(2,:)];
    cmapS1=[cmapbrain(1,:);cmapbrain(1,:);cmapbrain(1,:);cmapbrain(1,:)];
    
    
    subplot(3,2,mk)
    mybar(TaskDep.PM,'SD',[0 10],'unpaired',{'2','3','4','5'},'VAF by TD',cmapPM);
    if mk==1
        title(monkey(mk))
    elseif mk==2
        title(monkey(mk))
    end
    subplot(3,2,mk+2)
    mybar( TaskDep.M1,'SD',[0 10],'unpaired',{'2','3','4','5'},'VAF by TD',cmapM1);
    subplot(3,2,mk+4)
    mybar( TaskDep.S1,'SD',[0 10],'unpaired',{'2','3','4','5'},'VAF by TD',cmapS1);
    
    
end


%% ANOTHER AMAZING CODE TO CHECK THE TD VAF WITH THE INTRODUCTION OF MORE TASKS

clear all
close all
clc

condition={'2tasks/','3tasks/','4tasks/','5tasks/'};% for file directory
conditionNbre=['2';'3';'4';'5'];
monkey = {'Elektra','Natalya'};


for mk=1:2
    % store data in a big matrix
    TaskDep.PM=nan(7,4); % 7 days, 4 conditions comparison: 1 task splitted into: 2, 3, 4 and 5
    TaskDep.M1=nan(7,4);
    TaskDep.S1=nan(7,4);   
    
    for cond=1:length(conditionNbre)
        
        DataDir=(['/Users/Borgo/Documents/MATLAB/fri_mks/BSI LOCOMOTION Project_end_of_july_2021/Research Material/Analysis/PROCESSED DATA/' monkey{mk} '/dPCAmatricesCtrlTOM/' condition{cond}]);
        Files=dir([DataDir '*.mat']);
        
        
        for n=1:numel(Files) %loop across the varMat files
            load([DataDir Files(n).name])
            
            for day = 1:numel(varMat)
                tempPM(day,n)=varMat(day).PMind;
                tempM1(day,n)=varMat(day).M1ind;
                tempS1(day,n)=varMat(day).S1ind;
            end

        end     
            TaskDep.PM(:,cond)=mean(tempPM,2);
            TaskDep.M1(:,cond)=mean(tempM1,2);
            TaskDep.S1(:,cond)=mean(tempS1,2);
            clear tempS1 tempM1 tempPM
    end
    
   %Colormap of brain
    cmapbrain=[158,55,123;... %[82,163,204]/255; % S1
        238,98,109;... %[46,115,69]/255; % M1
        61,184,224]; %[168,126,229]/255; % PM
    cmapbrain=cmapbrain/255;
    cmapPM=[cmapbrain(3,:);cmapbrain(3,:);cmapbrain(3,:);cmapbrain(3,:)];
    cmapM1=[cmapbrain(2,:);cmapbrain(2,:);cmapbrain(2,:);cmapbrain(2,:)];
    cmapS1=[cmapbrain(1,:);cmapbrain(1,:);cmapbrain(1,:);cmapbrain(1,:)];
    
    
    subplot(3,2,mk)
    mybar(TaskDep.PM,'SD',[0 100],'unpaired',{'2','3','4','5'},'VAF by TD',cmapPM);
    if mk==1
        title(monkey(mk))
    elseif mk==2
        title(monkey(mk))
    end
    subplot(3,2,mk+2)
    mybar( TaskDep.M1,'SD',[0 100],'unpaired',{'2','3','4','5'},'VAF by TD',cmapM1);
    subplot(3,2,mk+4)
    mybar( TaskDep.S1,'SD',[0 100],'unpaired',{'2','3','4','5'},'VAF by TD',cmapS1);
    
end






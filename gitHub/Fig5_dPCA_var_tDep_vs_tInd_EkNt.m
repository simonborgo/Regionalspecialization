close all
clear
load('ColorMapDarkHot.mat')
colormap_Cyber;
% SELECT ALL THE /Users/Borgo/Documents/MATLAB/fri_mks/BSI LOCOMOTION Project_end_of_july_2021/Research Material/Analysis/PROCESSED DATA/Elektra/MultiunitClean/TrainingMatricesSimon
% for EK: TrainingMatricesSimon
% for NT: TrainingMatrices

[Files,PathName,~] = uigetfile('*.mat','MultiSelect','on');
Files = Files';
origPath = cd;
FilesPath = PathName;
% FilesPath = 'D:\Data\Elektra_dPCA\';
saveDir = [PathName 'Analysis' filesep];
% Get Conditions

%% Get trial info
trialCount = 1;
for fileIndex = 1:length(Files)
   
    clear D A B avgTrial
    
    % Load data
    load([FilesPath Files{fileIndex,:}])
    unders = strfind(Files{fileIndex,:},'_');
        
    % Get day info for each file   
    fileDays{fileIndex,:} = Files{fileIndex,:}((unders(3)+1):(unders(4)-1));  
end

% Make the structure with as many days as there are
varMat = struct();
allDays = unique(fileDays);
nDays = size(unique(fileDays),1);

for i = 1:nDays
    varMat(i).Day = allDays{i,1};
end

%% Fill out the info for each day
for fileIndex = 1:length(Files)
    
    % Load data
    load([FilesPath Files{fileIndex,:}])
    unders = strfind(Files{fileIndex,:},'_');
    
     % Get variance explained by taskDep and taskInd
    varTask = Matrices.explVar.totalMarginalizedVar/Matrices.explVar.totalVar*100;
    % TaskDependent first, then TaskIndependent
    
    % Find correct day Index
    for idxDay = 1:numel(varMat)

        Day = Files{fileIndex,:}((unders(3)+1):(unders(4)-1));
        if strcmp(Day,varMat(idxDay).Day)
            break % exit when the day has the same as the variable
        end
    end

    correctDay = idxDay;

     % Get trial Info     
    varMat(correctDay).Monkey = Files{fileIndex,:}((unders(2)+1):(unders(3)-1));
    varMat(correctDay).Day = Files{fileIndex,:}((unders(3)+1):(unders(4)-1)); 
    
    if ~isempty(strfind(Files{fileIndex,:},'PM'))      
        varMat(correctDay).PMdep = varTask(1);
        varMat(correctDay).PMind = varTask(2);
    elseif ~isempty(strfind(Files{fileIndex,:},'M1'))
        varMat(correctDay).M1dep = varTask(1);
        varMat(correctDay).M1ind = varTask(2);
    elseif ~isempty(strfind(Files{fileIndex,:},'S1'))        
        varMat(correctDay).S1dep = varTask(1);
        varMat(correctDay).S1ind = varTask(2);
    end
end

%% Statistical comparison within and across cortices
for i = 1:nDays
    PM(i,:) = [varMat(i).PMdep, varMat(i).PMind];
    M1(i,:) = [varMat(i).M1dep, varMat(i).M1ind];
    S1(i,:) = [varMat(i).S1dep, varMat(i).S1ind];   % uncomment
end
    
figure
comp = [PM(:,1),M1(:,1),S1(:,1)];
comp2 = [PM(:,2),M1(:,2),S1(:,2)];

color = [color_PM(1,:); color_M1(1,:); color_S1(1,:)];

% sub(1) = subplot(1,2,1);
% [p_TaskDep] = mybar(comp,'SD',[0 120],'unpaired',{'PM','M1','S1'},'Task-Dependent',color);
% 
% sub(2) = subplot(1,2,2);
figure
[p_TaskInd] = mybar(comp2,'SD',[],'unpaired',{'PM','M1','S1'},'Task-Independent',color);
ylim([0 120])
set(gcf,'Position',[462   149   344   344],'Renderer','Painters')
set(gca,'TickDir','out')

% figure
% sub2(1) = subplot(1,3,1);
% [p_PM] = mybar(PM,'SD',[],'paired',{'Task-Dep','Task-Ind'},'premotor',color)
% sub2(2) = subplot(1,3,2);
% [p_M1] = mybar(M1,'SD',[],'paired',{'Task-Dep','Task-Ind'},'motor',color)
% sub2(3) = subplot(1,3,3);
% [p_S1] = mybar(S1,'SD',[],'paired',{'Task-Dep','Task-Ind'},'sensory',color)
% linkaxes(sub2,'y');
% set(gcf,'Position',[593   295   609   240],'Renderer','painters')

% Compare cortices
figure
[p_PMM1] = mybar(comp2(:,[1,2]),'SD',[0 1],'paired',{'PM','M1'},'Task-Independent',color)
set(gca,'TickDir','out');
figure
[p_M1S1] = mybar(comp2(:,[2,3]),'SD',[0 1],'paired',{'M1','S1'},'Task-Independent',color)
set(gca,'TickDir','out');
figure
[p_PMS1] = mybar(comp2(:,[1,3]),'SD',[0 1],'paired',{'PM','S1'},'Task-Independent',color)

avgVAF = mean(comp2)
stdVAF = std(comp2)

pause

%% 2 task comparison
% first, run m1_TrainingMatrices twice to get CORR-LADi and STA-OBS matrices
% second, save comp2 variable.
% third, run again this code to get comp2 for the other pair of task
%fourth, concatenate the 2 comp2 and mybar it
close all
compfinal=[comp2;comp2sta]; % I rename comp2 to comp2sta and load comp2 of CORR-LAD

% WARNING: I MODIFIED mybar line 142 to plot in black sta-obs and in grey
% corr-ladd
figure
[p_TaskInd] = mybar(compfinal,'SD',[],'unpaired',{'PM','M1','S1'},'Task-Independent',color);
ylim([0 120])
set(gcf,'Position',[462   149   344   344],'Renderer','Painters')
set(gca,'TickDir','out')
%
figure
[p_PMM1] = mybar(compfinal(:,[1,2]),'SD',[0 1],'paired',{'PM','M1'},'Task-Independent',color)
set(gca,'TickDir','out');
figure
[p_M1S1] = mybar(compfinal(:,[2,3]),'SD',[0 1],'paired',{'M1','S1'},'Task-Independent',color)
set(gca,'TickDir','out');
figure
[p_PMS1] = mybar(compfinal(:,[1,3]),'SD',[0 1],'paired',{'PM','S1'},'Task-Independent',color)

avgVAF = mean(compfinal)
stdVAF = std(compfinal)




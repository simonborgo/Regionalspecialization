close all
clear
load('ColorMapDarkHot.mat')
colormap_Cyber;
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
    end
end

%% Statistical comparison within and across cortices
for i = 1:nDays
    PM(i,:) = [varMat(i).PMdep, varMat(i).PMind];
    M1(i,:) = [varMat(i).M1dep, varMat(i).M1ind];
    %S1(i,:) = [varMat(i).S1dep, varMat(i).S1ind];
end
    
color = [color_PM(1,:); color_M1(1,:)];

figure
comp2 = [PM(:,2),M1(:,2)];
[p_TaskInd] = mybar(comp2,'SD',[],'unpaired',{'PM','M1'},'Task-Independent',color);
ylim([0 90])
set(gcf,'Position',[462   149   344   344],'Renderer','Painters')
set(gca,'TickDir','out')

figure
sub2(1) = subplot(1,2,1);
[p_PM] = mybar(PM,'SD',[],'unpaired',{'Task-Dep','Task-Ind'},'premotor',color);
figure
sub2(2) = subplot(1,2,2);
[p_M1] = mybar(M1,'SD',[],'unpaired',{'Task-Dep','Task-Ind'},'motor',color);

linkaxes(sub2,'y');
set(gcf,'Position',[593   295   609   240],'Renderer','painters')

avgVAF = mean(comp)
stdVAF = std(comp)

%% 
close all
compfinal=[comp2;comp2sta]; % I rename comp2 to comp2sta and load comp2 of CORR-LAD

% WARNING: I MODIFIED mybar line 120 & 143 to plot in black sta-obs and in grey
% corr-ladd
figure
[p_TaskInd] = mybar(compfinal,'SD',[],'unpaired',{'PM','M1'},'Task-Independent',color);
ylim([0 120])
set(gcf,'Position',[462   149   344   344],'Renderer','Painters')
set(gca,'TickDir','out')
%
figure
[p_PMM1] = mybar(compfinal(:,[1,2]),'SD',[0 1],'paired',{'PM','M1'},'Task-Independent',color)
set(gca,'TickDir','out');

avgVAF = mean(compfinal)
stdVAF = std(compfinal)


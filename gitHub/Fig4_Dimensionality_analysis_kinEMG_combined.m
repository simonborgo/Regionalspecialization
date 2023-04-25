%% dimensionality analysis of Natalya
% same as before but automatically loads the Ds and do the stuff
% results saved in the same way as before

%% Dimensionality and his relation to task.
% which are the task that drive the increase in the dimensionality?
% we conduct the same analysis as before but removing one task at a time
%% FOR ELEKTRA MUSCLES
clear all 
close all
Monkey = 'Elektra';
load(['sessions' Monkey])

for Day = 1:size(sessions,1) % length(sessions)
    clearvars -except Day sessions Monkey EMGsalldays
    close all
    
    
    method = 1; % 1 for using the PC coefficients, 2 for using the SCORE, 3 using actual EMGs
    f_timeBased = 0; % Flag to do time-based rather than time-wrapped analysis
    f_norm = 1; % Flag to normalize the correlation by the within-task correlation
    f_plotHist = 0; % Flag to also plot histograms for distribution
    f_dimAnalysis = 1; % Flag to perform a dimensionality analysis on PCA data
    % muscles = [2,4,5,7,8]; % IL, ST, EDL, MG, FHL
    
    muscles = [2:8]; % without left leg and right gluteus (broken electrodes)
    %muscles = [1:10]; % All
    mDim = 15; % Number of dimensions for analysis
    mDimMusc = length(muscles);
    
    %ctxTitles = {'Premotor';'Motor';'Sensory'};
    
    warning('off','stats:pca:ColRankDefX')
    
    %%%%%%%%%%% Multiunit
    
    % sessions = ['20190415';'20190425'];
     dataDir = ['/Users/Borgo/Google Drive/BSI LOCOMOTION Project/Research Material/Analysis/PROCESSED DATA/' Monkey '/MultiUnitClean/'];

%     dataDir = ['/Users/ismaelseanez/Google Drive/Working Folder/G-Lab/MONKEY/'...
%         'BSI LOCOMOTION Project/Research Material/Analysis/PROCESSED DATA/' Monkey '/MultiUnitClean/'];
    % dataDir = ['C:\Users\physio\Google Drive\BSI LOCOMOTION Project\Research Material\Analysis\PROCESSED DATA\Natalya\Multiunit\'];
%           dataDir = 'C:\Users\Phoenix\Google Drive\BSI LOCOMOTION Project\Research Material\Analysis\PROCESSED DATA\Elektra\MultiunitClean\';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%% Single Unit
    % sessions = ['20190425';'20190425'];
    %     dataDir = ['/Users/ismaelseanez/Google Drive/Working Folder/G-Lab/MONKEY/'...
    %         'BSI LOCOMOTION Project/Research Material/Analysis/PROCESSED DATA/Elektra/SingleUnit/'];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    nDays = size(sessions,1);
    Results = struct();

    
    % Load important data and run functions
    load('ColorMapDarkHot.mat')
    colormap_Cyber;
    %color_ctx = [color_PM(1,:);color_M1(1,:);color_S1(1,:)];
    
    Results(Day).Day = sessions(Day,:);
    S = dir([dataDir filesep Monkey '_' sessions(Day,:) '*.mat']);
    Files = cell(numel(S),1);
    for i = 1:numel(S)
        Files{i,:} = S(i).name;
    end
    
    origPath = cd;
    FilesPath = dataDir;
    
    % Get Conditions
    for trialIndex = 1:length(Files)
        unders = strfind(Files{trialIndex,:},'_');
        DataAvg(trialIndex).condition = Files{trialIndex,:}((unders(end-1)+1):(unders(end)-1));
    end
    
    k = 1;
    nanflag = 0;
    for trialIndex = 1:length(Files)
        clear D A B avgTrial
        
        % Load data
        load([FilesPath Files{trialIndex,:}])
        
        for j = 1:numel(D)           
            EMG_cyclic(:,:,j) = D(j).EMG_cyclic(:,muscles);
            DataAllGaits(k).duration = length(D(j).times)/100;
            DataAllGaits(k).FOpercent = D(j).RFO/size(D(j).A,1)*100;
            DataAllGaits(k).EMG_cyclic = D(j).EMG_cyclic(:,muscles);
            DataAllGaits(k).condition = DataAvg(trialIndex).condition;
            
            k = k+1;
        end
        
    end
    
    % Find all gait cycles of each condition
    nCond = 5;
    [nGaits, idx_Cond,condNames] = f_findGaitConditions(DataAllGaits,nCond);
    
    %% Find outliers for each condition
    % Based on gait durations and FO% for each condition
    allOut = zeros(1,length(idx_Cond));
    for i = 1:nCond
        iCond = find(idx_Cond(i,:)==1);
        durations = zeros(1,length(iCond));
        FO = zeros(1,length(iCond));
        for j = 1:length(iCond)
            durations(j) = DataAllGaits(iCond(j)).duration;
            FO(j) = DataAllGaits(iCond(j)).FOpercent;
        end
        
        % Flag Outliers
        gaitOut1 = find(durations>1.5); % Detect steps longer than 1.5 sec
        durations(gaitOut1) = NaN; % And ignore them for the distribution
        FO(gaitOut1) = NaN;
        
        [~,gaitOut2] = outlieriqr(durations); % Find outlier gaits based on time
        [~,gaitOut3] = outlieriqr(FO); % And foot-off location
        
        % Flag them for removal
        allOut(iCond(gaitOut1)) = 1;
        allOut(iCond(gaitOut2)) = 1;
        allOut(iCond(gaitOut3)) = 1;
    end
    
    toRemove = flip(find(allOut==1)); % We flip them because the index will
    % change each time we remove one, so we start from the back
    
    % Remove outlier trials
    for i = 1:length(toRemove)
        DataAllGaits(toRemove(i)) = []; % Now DataAllGaits has clean gaits
    end
    
    % Find all gait cycles of each condition again
    nCond = 5;
    [~, idx_Cond,~] = f_findGaitConditions(DataAllGaits,nCond);
        
    
    
    
    %% Get average behavior for each condition

    nEMG = size(DataAllGaits(1).EMG_cyclic,2);
    lengthEMG = size(DataAllGaits(1).EMG_cyclic,1);
    
    % Get the number of gait cycles for each condition
    for i = 1:nCond
        nGaitsCond(i) = sum(idx_Cond(i,:));
    end
    
    maxGaits = max(nGaitsCond); % Maximum number of gaits
    
    % Create firingRate and EMG arrays
    EMGs = nan(nEMG,nCond,lengthEMG,maxGaits);
    
    % Populate arrays
    countCond = ones(1,nCond);
    for gait = 1:numel(DataAllGaits)
        cond = find(idx_Cond(:,gait)==1);
        if ~isempty(cond)
            EMGs(:,cond,:,countCond(cond)) = DataAllGaits(gait).EMG_cyclic';
            countCond(cond) = countCond(cond)+1;
        end
    end
    
    %save EMGs 
    days={'a','b','c','d','e','f','g'};
    EMGsalldays.(days{Day})= EMGs;
    
 %    clear dim PM M1 S1 A_S1 A_PM B_PM B_S1 B_pred DataAllGaits Neural_cyclic  M1avg PMavg Dim dimension PM_cond M1_cond S1_cond S1avg
 %    close all
end % The end of days                 (heyoooo!)


%% FOR NATALYA KINEMATICS
%% Dimensionality and his relation to task.
% which are the task that drive the increase in the dimensionality?
% we conduct the same analysis as before but removing one task at a time

Monkey = 'Natalya';
load(['sessions' Monkey])


for Day = 1:size(sessions,1) % length(sessions)
    clearvars -except Day sessions Monkey EMGsalldays Kinsalldays muscles
    close all
    
    
    method = 1; % 1 for using the PC coefficients, 2 for using the SCORE, 3 using actual EMGs
    f_timeBased = 0; % Flag to do time-based rather than time-wrapped analysis
    f_norm = 1; % Flag to normalize the correlation by the within-task correlation
    f_plotHist = 0; % Flag to also plot histograms for distribution
    f_dimAnalysis = 1; % Flag to perform a dimensionality analysis on PCA data
    % muscles = [2,4,5,7,8]; % IL, ST, EDL, MG, FHL
    
   load('GoodKinematics.mat');
    %goodKin=[1:18];%All
    
    mDim = 15; % Number of dimensions for analysis
    mDimKin = length(goodKin);
    
    %ctxTitles = {'Premotor';'Motor';'Sensory'};
    
    warning('off','stats:pca:ColRankDefX')
    
    %%%%%%%%%%% Multiunit
    
    % sessions = ['20190415';'20190425'];
     dataDir = ['/Users/Borgo/Google Drive/BSI LOCOMOTION Project/Research Material/Analysis/PROCESSED DATA/' Monkey '/MultiunitClean/'];


    nDays = size(sessions,1);
    Results = struct();

    
    % Load important data and run functions
    load('ColorMapDarkHot.mat')
    colormap_Cyber;
    
    Results(Day).Day = sessions(Day,:);
    S = dir([dataDir filesep Monkey '_' sessions(Day,:) '*.mat']);
    Files = cell(numel(S),1);
    for i = 1:numel(S)
        Files{i,:} = S(i).name;
    end
    
    origPath = cd;
    FilesPath = dataDir;
    
    % Get Conditions
    for trialIndex = 1:length(Files)
        unders = strfind(Files{trialIndex,:},'_');
        DataAvg(trialIndex).condition = Files{trialIndex,:}((unders(end-1)+1):(unders(end)-1));
    end
    
    k = 1;
    nanflag = 0;
    for trialIndex = 1:length(Files)
        clear D A B avgTrial
        
        % Load data
        load([FilesPath Files{trialIndex,:}])
        
        for j = 1:numel(D)
            
            Kin_cyclic(:,:,j) = D(j).Kin_cyclic(:,goodKin);
            DataAllGaits(k).duration = length(D(j).times)/100;
            DataAllGaits(k).FOpercent = D(j).RFO/size(D(j).A,1)*100;
            DataAllGaits(k).Kin_cyclic = D(j).Kin_cyclic(:,goodKin);
            DataAllGaits(k).condition = DataAvg(trialIndex).condition;
            
            k = k+1;
        end
        
    end
    
    % Find all gait cycles of each condition
    nCond = 5;
    [nGaits, idx_Cond,condNames] = f_findGaitConditions(DataAllGaits,nCond);
    
    %% Find outliers for each condition
    % Based on gait durations and FO% for each condition
    allOut = zeros(1,length(idx_Cond));
    for i = 1:nCond
        iCond = find(idx_Cond(i,:)==1);
        durations = zeros(1,length(iCond));
        FO = zeros(1,length(iCond));
        for j = 1:length(iCond)
            durations(j) = DataAllGaits(iCond(j)).duration;
            FO(j) = DataAllGaits(iCond(j)).FOpercent;
        end
        
        % Flag Outliers
        gaitOut1 = find(durations>1.5); % Detect steps longer than 1.5 sec
        durations(gaitOut1) = NaN; % And ignore them for the distribution
        FO(gaitOut1) = NaN;
        
        [~,gaitOut2] = outlieriqr(durations); % Find outlier gaits based on time
        [~,gaitOut3] = outlieriqr(FO); % And foot-off location
        
        % Flag them for removal
        allOut(iCond(gaitOut1)) = 1;
        allOut(iCond(gaitOut2)) = 1;
        allOut(iCond(gaitOut3)) = 1;
    end
    
    toRemove = flip(find(allOut==1)); % We flip them because the index will
    % change each time we remove one, so we start from the back
    
    % Remove outlier trials
    for i = 1:length(toRemove)
        DataAllGaits(toRemove(i)) = []; % Now DataAllGaits has clean gaits
    end
    
    % Find all gait cycles of each condition again
    nCond = 5;
    [~, idx_Cond,~] = f_findGaitConditions(DataAllGaits,nCond);
        
    
    
    
    %% Get average behavior for each condition

    nKin = size(DataAllGaits(1).Kin_cyclic,2);
    lengthKin = size(DataAllGaits(1).Kin_cyclic,1);
    
    % Get the number of gait cycles for each condition
    for i = 1:nCond
        nGaitsCond(i) = sum(idx_Cond(i,:));
    end
    
    maxGaits = max(nGaitsCond); % Maximum number of gaits
    
    % Create KIN arrays
    Kins = nan(nKin,nCond,lengthKin,maxGaits);
    
    % Populate arrays
    countCond = ones(1,nCond);
    for gait = 1:numel(DataAllGaits)
        cond = find(idx_Cond(:,gait)==1);
        if ~isempty(cond)
            
                Kins(:,cond,:,countCond(cond)) = DataAllGaits(gait).Kin_cyclic';
            countCond(cond) = countCond(cond)+1;
        end
    end
    
    %save Kins 
    days={'a','b','c','d','e','f','g'};
    Kinsalldays.(days{Day})= Kins;
end
    
    %% %%%%%%%%%%%% MANIFOLD ANAYSIS %%%%%%%%%%%%%%%%%%%%%%%%%
for Day = 1:size(sessions,1)
    combos = nchoosek(1:nCond,2); % All possible combinations to test
    minSteps = min(countCond-1); % Minimum number of steps on one condition

    
    %get EMG
        EMG = EMGsalldays.(days{Day});
        EMG_avg = mean(EMG(:,:,:,:),4,'omitnan');
        EMG_avg = bsxfun(@minus,EMG_avg,mean(EMG_avg,3));

        EMGavg = [];
        
     %get KIN
        Kin = Kinsalldays.(days{Day});
        Kin_avg = mean(Kin(:,:,:,:),4,'omitnan');
        Kin_avg = bsxfun(@minus,Kin_avg,mean(Kin_avg,3));

        Kinavg = [];
        
      % combine kin and EMGs
        MVT_avg=[EMG_avg; Kin_avg];
        
         MVTavg = [];
         
    %combined space across tasks
        for chan = 1:length(muscles)+length(goodKin)
            MVTfinal = [];
            
            for cond = 1:nCond
               MVTfinal = [MVTfinal, squeeze(MVT_avg(chan,cond,:))'];
                
            end
            MVTavg = [MVTavg; MVTfinal];
            
        end

     %one space per task
        for cond = 1:nCond
            
           
            %MVT_cond = MVT(:,cond,:,1:countCond(cond)-1);
            
            MVTavg_cond = squeeze(MVT_avg(:,cond,:));
            
            dimMVTall(cond) = m_DimensionalityEstimate_Simple(MVTavg_cond');
            

        end

         dimension = [dimMVTall];
         Dim.MVT = m_DimensionalityEstimate_Simple(MVTavg');


    saveDir = [dataDir 'DimEstimatesMVT' filesep];
    if ~isdir(saveDir)
        mkdir(saveDir)
    end
    save([saveDir Monkey '_' sessions(Day,:)],'dimension','Dim')
    
    clear dim PM M1 S1 A_S1 A_PM B_PM B_S1 B_pred DataAllGaits Neural_cyclic  M1avg PMavg Dim dimension PM_cond M1_cond S1_cond S1avg
    close all
end % The end of days                 (heyoooo!)

%a_DimensionalityEstimates_Histograms_kinEMG


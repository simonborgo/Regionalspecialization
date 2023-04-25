%% dimensionality analysis of Natalya
% same as before but automatically loads the Ds and do the stuff
% results saved in the same way as before

%% Dimensionality and his relation to task.
% which are the task that drive the increase in the dimensionality?
% we conduct the same analysis as before but removing one task at a time
clear all
close all
Monkey = 'Kara';

load(['sessions' Monkey])
% sessions = sessionsObs;

for Day = 1:size(sessions,1)
    close all
      
    ctxTitles = {'Premotor';'Motor';'Sensory'};
    
    warning('off','stats:pca:ColRankDefX')
    
    %%%%%%%%%%% Multiunit
    
%     dataDir = ['/Users/ismaelseanez/Google Drive/Working Folder/G-Lab/MONKEY/'...
%         'BSI LOCOMOTION Project/Research Material/Analysis/PROCESSED DATA/' monkey '/MultiUnitClean/'];
    
    dataDir = ['/Users/Borgo/Google Drive/BSI LOCOMOTION Project/Research Material/Analysis/PROCESSED DATA/' Monkey '/MultiunitClean/'];
    
%     dataDir = ['C:\Users\physio\Google Drive\BSI LOCOMOTION Project\Research Material\Analysis\PROCESSED DATA\' ...
%         monkey '\MultiunitClean\'];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    nDays = size(sessions,1);
    
    % Load important data and run functions
    load('ColorMapDarkHot.mat')
    colormap_Cyber;
    color_ctx = [color_PM(1,:);color_M1(1,:);color_S1(1,:)];
    
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
            
            Neural_cyclic(:,:,j) = D(j).Neural_cyclic;
            %         EMG_cyclic(:,:,j) = D(j).EMG_cyclic;
            DataAllGaits(k).duration = length(D(j).times)/100;
            DataAllGaits(k).FOpercent = D(j).RFO/size(D(j).A,1)*100;
            DataAllGaits(k).Neural_cyclic = D(j).Neural_cyclic;
            %         DataAllGaits(k).EMG_cyclic = D(j).EMG_cyclic;
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
   % nCond = 2;
    [~, idx_Cond,~] = f_findGaitConditions(DataAllGaits,nCond);
    
    %% Get average behavior for each condition
    nChan = size(DataAllGaits(1).Neural_cyclic,2);
    lengthNeural = size(DataAllGaits(1).Neural_cyclic,1);
    
    % Get the number of gait cycles for each condition
    for i = 1:nCond
        nGaitsCond(i) = sum(idx_Cond(i,:));
    end
    
    maxGaits = max(nGaitsCond); % Maximum number of gaits
    
    % Create firingRate and EMG arrays
    firingRatesAll = nan(nChan,nCond,lengthNeural,maxGaits);
    % EMGs = nan(nEMG,nCond,lengthEMG,maxGaits);
    
    % Populate arrays
    countCond = ones(1,nCond);
    for gait = 1:numel(DataAllGaits)
        cond = find(idx_Cond(:,gait)==1);
        if ~isempty(cond)
            
            firingRatesAll(:,cond,:,countCond(cond)) = DataAllGaits(gait).Neural_cyclic';
            %     EMGs(:,cond,:,countCond(cond)) = DataAllGaits(gait).EMG_cyclic';
            countCond(cond) = countCond(cond)+1;
        end
    end
    
    % Get the average for each condition
    firingRatesAvg = mean(firingRatesAll,4,'omitnan');
    firingRatesSTE = std(firingRatesAll,0,4,'omitnan')./sqrt(size(firingRatesAll,4));
    % EMGsAvg = mean(EMGs,4,'omitnan');
    % EMGsSTE = std(EMGs,0,4,'omitnan')./sqrt(size(EMGs,4));
    
    % Get the normalized FR
    firingRatesNorm = bsxfun(@rdivide,firingRatesAvg,max(firingRatesAvg,[],3));
    % firingRatesNorm = firingRatesAvg;
    
    % EMGsNorm = bsxfun(@rdivide,EMGsAvg,max(EMGsAvg,[],3));
    
    % Get the peak location for each condition for sorting
    if nChan == 96 % Number of channels 
    nPM = 48; % Number of neurons on each cortex
    nM1 = 48;
    elseif nChan <96
        nPM = Data.nCellPM; % Number of neurons on each cortex
        nM1 = Data.nCellM1;
    end
    
    %% Average setup
    % other parameters needed
    removeCond = find(nGaitsCond==0); % Find conditions with no gaits
    nCondGood = length(nGaitsCond)-length(removeCond);
    goodCond = 1:nCond;
    goodCond(removeCond) = [];
    countCond(removeCond) = [];
    nCond = length(goodCond);
    
    firingRatesAll(:,removeCond,:,:) = [];
    firingRatesAvg(:,removeCond,:) = [];
    
    % Populate trial array
    trialLength = size(DataAllGaits(1).Neural_cyclic,1);
    trialNumAll = nan(nChan,nCond);
    for i = 1:nCond
        trialNumAll(:,i) = countCond(i)-1;
    end
    

    %% %%%%%%%%%%%% MANIFOLD ANAYSIS %%%%%%%%%%%%%%%%%%%%%%%%%
    
    combos = nchoosek(1:nCond,2); % All possible combinations to test
    minSteps = min(countCond-1); % Minimum number of steps on one condition

    neurons_PM = 1:nPM;
    neurons_M1 = nPM+1:nPM+nM1;
        
    % Firing Rates for motor and premotor
        PM = firingRatesAll(neurons_PM,:,:,:);
        PM_avg = mean(firingRatesAll(neurons_PM,:,:,:),4,'omitnan');
        M1 = firingRatesAll(neurons_M1,:,:,:);
        M1_avg = mean(firingRatesAll(neurons_M1,:,:,:),4,'omitnan');

    % Subtract the mean of each condition for each channel
    PM = bsxfun(@minus,PM,mean(PM_avg,3));
    PM_avg = bsxfun(@minus,PM_avg,mean(PM_avg,3));
    M1 = bsxfun(@minus,M1,mean(M1_avg,3));
    M1_avg = bsxfun(@minus,M1_avg,mean(M1_avg,3));
        
    PMavg = [];
    M1avg = [];
    
    for chan = 1:length(neurons_PM)
        PMfinal = [];
        for cond = 1:nCond
            PMfinal = [PMfinal, squeeze(PM_avg(chan,cond,:))'];
        end
        PMavg = [PMavg; PMfinal];
    end
    
    for chan = 1:length(neurons_M1)
        M1final = [];
        for cond = 1:nCond
            M1final = [M1final, squeeze(M1_avg(chan,cond,:))'];
        end
        M1avg = [M1avg; M1final];
    end
        
    % Communication space between cortices per task
        for cond = 1:nCond
            
            PM_cond = PM(:,cond,:,1:countCond(cond)-1);
            PMavg_cond = squeeze(PM_avg(:,cond,:));
            dimPMall(cond) = m_DimensionalityEstimate_Simple(PMavg_cond');
            
            M1_cond = M1(:,cond,:,1:countCond(cond)-1);
            M1avg_cond = squeeze(M1_avg(:,cond,:));
            dimM1all(cond) = m_DimensionalityEstimate_Simple(M1avg_cond');

        end
    dimension = [dimPMall; dimM1all];
    
    Dim.PM = m_DimensionalityEstimate_Simple(PMavg');
    Dim.M1 = m_DimensionalityEstimate_Simple(M1avg');
    
    saveDir = [dataDir 'DimEstimates' filesep];
    if ~isdir(saveDir)
        mkdir(saveDir)
    end
    save([saveDir Monkey '_' sessions(Day,:)],'dimension','Dim')
    
    clear dim PM M1 S1 A_S1 A_PM B_PM B_S1 B_pred DataAllGaits Neural_cyclic  M1avg PMavg Dim dimension PM_cond M1_cond S1_cond S1avg
    close all
end % The end of days                 (heyoooo!)
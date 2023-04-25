%% dimensionality analysis of Natalya
% same as before but automatically loads the Ds and do the stuff
% results saved in the same way as before
clear all

% WARNING: CHANGE DIRECTORY IN a_DimensionalityEstimates_Histograms

%% Dimensionality and his relation to task.
% which are the task that drive the increase in the dimensionality?
% we conduct the same analysis as before but removing one task at a time
Monkey = 'Natalya';
load(['sessions' Monkey])

for Day = 1:size(sessions,1) % length(sessions)
    clearvars -except Day sessions Monkey
    close all
    
    
    method = 1; % 1 for using the PC coefficients, 2 for using the SCORE, 3 using actual EMGs
    f_timeBased = 0; % Flag to do time-based rather than time-wrapped analysis
    f_norm = 1; % Flag to normalize the correlation by the within-task correlation
    f_plotHist = 0; % Flag to also plot histograms for distribution
    f_dimAnalysis = 0; % Flag to perform a dimensionality analysis on PCA data
    % muscles = [2,4,5,7,8]; % IL, ST, EDL, MG, FHL
    
    muscles = [1:10]; % All
    mDim = 15; % Number of dimensions for analysis
    mDimMusc = length(muscles);
    
    ctxTitles = {'Premotor';'Motor';'Sensory'};
    
    warning('off','stats:pca:ColRankDefX')
    
    %%%%%%%%%%% Multiunit
    
    % sessions = ['20190415';'20190425'];
    dataDir = ['/Users/Borgo/Documents/MATLAB/fri_mks/BSI LOCOMOTION Project_end_of_july_2021/Research Material/Analysis/PROCESSED DATA/' Monkey '/MultiunitClean/'];

%     dataDir = ['//Users/Borgo/Google Drive/BSI LOCOMOTION Project/Research Material/Analysis/PROCESSED DATA/' Monkey '/MultiunitClean/'];
%     dataDir = ['/Users/ismaelseanez/Google Drive/Working Folder/G-Lab/MONKEY/'...
%         'BSI LOCOMOTION Project/Research Material/Analysis/PROCESSED DATA/' Monkey '/MultiunitClean/'];
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
    all_prinAngles = zeros(nDays,3,7,7,mDim);
    all_prinAnglesSurr = zeros(nDays,3,7,7,mDim);
    EMG_prinAngles = zeros(nDays,7,7,mDimMusc);
    EMG_prinAnglesSurr = zeros(nDays,7,7,mDimMusc);
    
    % Load important data and run functions
    load('ColorMapDarkHot.mat')
    colormap_Cyber;
    color_ctx = [color_PM(1,:);color_M1(1,:);color_S1(1,:)];
    
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
    nCond = 5;
    [~, idx_Cond,~] = f_findGaitConditions(DataAllGaits,nCond);
    
    %% Get average behavior for each condition
    nChan = size(DataAllGaits(1).Neural_cyclic,2);
    % if nChan>nPM+nM1
    %     f_sensory = 1; % Check if there's sensory data
    % else
    %     f_sensory = 0;
    % end
    lengthNeural = size(DataAllGaits(1).Neural_cyclic,1);
    % nEMG = size(DataAllGaits(1).EMG_cyclic,2);
    % lengthEMG = size(DataAllGaits(1).EMG_cyclic,1);
    
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
    
    % Demean the data
%     firingRatesAll = bsxfun(@minus,firingRatesAll,mean(firingRatesAll,3,'omitnan'));
%     FR_Avg = bsxfun(@minus,firingRatesAvg,mean(firingRatesAvg,3,'omitnan'));
    % Get the average for each condition
    %     firingRatesAvg = mean(firingRatesAll,4,'omitnan');
    %     firingRatesSTE = std(firingRatesAll,0,4,'omitnan')./sqrt(size(firingRatesAll,4));
    %     % EMGsAvg = mean(EMGs,4,'omitnan');
    %     % EMGsSTE = std(EMGs,0,4,'omitnan')./sqrt(size(EMGs,4));
    %
    %     % Get the normalized FR
    %     firingRatesNorm = bsxfun(@rdivide,firingRatesAvg,max(firingRatesAvg,[],3));
    %     % firingRatesNorm = firingRatesAvg;
    
    % EMGsNorm = bsxfun(@rdivide,EMGsAvg,max(EMGsAvg,[],3));
    
    % Get the peak location for each condition for sorting
    if strcmp(Monkey,'Elektra')
        if nChan == 160 % Number of channels
            nPM = 48; % Number of neurons on each cortex
            nM1 = 48;
            nS1 = 64;
        elseif nChan <160
            nPM = Data.nCellPM; % Number of neurons on each cortex
            nM1 = Data.nCellM1;
            nS1 = Data.nCellS1;
        end
    else
        if nChan == 160 % Number of channels
            nPM = 48; % Number of neurons on each cortex
            nM1 = 48;
            nS1 = 64;
        elseif nChan <160
            nPM = Data.nCellPM; % Number of neurons on each cortex
            nM1 = Data.nCellM1;
            nS1 = Data.nCellS1;
        elseif nChan == 192
            nS1 = 32;
            nM1 = 64;
            nPE = 32;
            nPM = 64;
        end
    end
    
    %     nChan = nChan - length(IdxElectOut);
    %% Bad channel deletion
    if strcmp(Monkey,'Elektra')
       
        Session = Files{1,1}(9:16);
        MkName = Files{1,1}(1:7);

        H.PMM1_badch = load(['badChan_' MkName '_M1PM_' Session '.mat']);
        H.S1_badch = load(['badChan_' MkName '_S1_' Session '.mat']);
        H.PMM1_badch.channels(H.PMM1_badch.channels == 0) = []; 
        H.S1_badch.channels(H.S1_badch.channels == 0) = []; 

        % Get number of bad guys
        nBadPM = length(find(H.PMM1_badch.channels<=48));
        nBadM1 = length(find(H.PMM1_badch.channels>48));
        nBadS1 = length(H.S1_badch.channels);

        % Remove them
        firingRatesAll(H.S1_badch.channels+96,:,:,:) = [];
        firingRatesAll(H.PMM1_badch.channels,:,:,:) = [];
%         firingRatesNorm(H.S1_badch.channels+96,:,:) = [];
%         firingRatesNorm(H.PMM1_badch.channels,:,:) = [];

        nPM = nPM - nBadPM;
        nM1 = nM1 - nBadM1;
        nS1 = nS1 - nBadS1;

        nChan = size(firingRatesAll,1);
    elseif strcmp(Monkey,'Natalya')
        Session = Files{1,1}(9:16);
        MkName = Files{1,1}(1:7);

        H.S1M1_badch = load(['badChan_' MkName '_S1M1_' Session '.mat']);
        H.PM_badch = load(['badChan_' MkName '_PM_' Session '.mat']);
        H.S1M1_badch.channels(H.S1M1_badch.channels == 0) = []; 
        H.PM_badch.channels(H.PM_badch.channels == 0) = []; 

        % Get number of bad guys
        nBadS1 = length(find(H.S1M1_badch.channels<=32));
        nBadM1 = length(find(H.S1M1_badch.channels>32));
        nBadPM = length(H.PM_badch.channels);

        % Remove them
        firingRatesAll(H.PM_badch.channels+nS1+nM1+nPE,:,:,:) = [];
        firingRatesAll(H.S1M1_badch.channels,:,:,:) = [];
%         firingRatesNorm(H.PM_badch.channels+nS1+nM1+nPE,:,:) = [];
%         firingRatesNorm(H.S1M1_badch.channels,:,:) = [];

        nPM = nPM - nBadPM;
        nM1 = nM1 - nBadM1;
        nS1 = nS1 - nBadS1;

        nChan = size(firingRatesAll,1);
    end
    
    %% Average setup
    % other parameters needed
    removeCond = find(nGaitsCond==0); % Find conditions with no gaits
    nCondGood = length(nGaitsCond)-length(removeCond);
    goodCond = 1:nCond;
    goodCond(removeCond) = [];
    
    % Populate trial array
    trialLength = size(DataAllGaits(1).Neural_cyclic,1);
    trialNumAll = nan(nChan,nCond);
    for i = 1:nCond
        trialNumAll(:,i) = countCond(i)-1;
    end
    
    
    %% %%%%%%%%%%%% MANIFOLD ANAYSIS %%%%%%%%%%%%%%%%%%%%%%%%%
    
    combos = nchoosek(1:nCond,2); % All possible combinations to test
    minSteps = min(countCond-1); % Minimum number of steps on one condition
    
    if strcmp(Monkey,'Natalya')
        neurons_S1 = 1:nS1;
        neurons_M1 = nS1+1:nS1+nM1;
        %    neurons_PE = nS1+nM1+1:nS1+nM1+nPE;
        neurons_PM = nS1+nM1+nPE+1:nS1+nM1+nPE+nPM;
    else
        neurons_PM = 1:nPM;
        neurons_M1 = nPM+1:nPM+nM1;
        neurons_S1 = nPM+nM1+1:nPM+nM1+nS1;
    end

        % Firing Rates for motor and premotor
        PM = firingRatesAll(neurons_PM,:,:,:);
        PM_avg = mean(firingRatesAll(neurons_PM,:,:,:),4,'omitnan');
        M1 = firingRatesAll(neurons_M1,:,:,:);
        M1_avg = mean(firingRatesAll(neurons_M1,:,:,:),4,'omitnan');
        S1 = firingRatesAll(neurons_S1,:,:,:);
        S1_avg = mean(firingRatesAll(neurons_S1,:,:,:),4,'omitnan');


        %Inspect all trials within individual conditions
%         f_inspect = 0;
%         if f_inspect == 1
%             figure
%             condInspect = 1;
%             a = squeeze(PM(:,condInspect,:,1:countCond(condInspect)-1));
%             for i = 1:countCond(condInspect)-1
%                 if i <= length(color_Centered)/2
%                     colIdx = i;
%                 else
%                     colIdx = i-length(color_Centered)/2;
%                 end
%                 plot([100*i-99:100*i],squeeze(a(:,:,i)),'color',color_Centered(colIdx*2,:))
%                 hold on
%             end
%             set(gcf,'Position',[60         152        1192         420],'Renderer','Painters')
%         end

        % Subtract the mean of each condition for each channel
        PM = bsxfun(@minus,PM,mean(PM_avg,3));
        PM_avg = bsxfun(@minus,PM_avg,mean(PM_avg,3));
        M1 = bsxfun(@minus,M1,mean(M1_avg,3));
        M1_avg = bsxfun(@minus,M1_avg,mean(M1_avg,3));
        S1 = bsxfun(@minus,S1,mean(S1_avg,3));
        S1_avg = bsxfun(@minus,S1_avg,mean(S1_avg,3));

        PMavg = [];
        M1avg = [];
        S1avg = [];
        %    PEavg = [];

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

        for chan = 1:length(neurons_S1)
            S1final = [];
            for cond = 1:nCond
                S1final = [S1final, squeeze(S1_avg(chan,cond,:))'];
            end
            S1avg = [S1avg; S1final];
        end
        
        % Communication space between cortices per task
        for cond = 1:nCond
            
            PM_cond = PM(:,cond,:,1:countCond(cond)-1);
            PMavg_cond = squeeze(PM_avg(:,cond,:));
            dimPMall(cond) = m_DimensionalityEstimate_Simple(PMavg_cond');
            
            M1_cond = M1(:,cond,:,1:countCond(cond)-1);
            M1avg_cond = squeeze(M1_avg(:,cond,:));
            dimM1all(cond) = m_DimensionalityEstimate_Simple(M1avg_cond');
            
            S1_cond = S1(:,cond,:,1:countCond(cond)-1);
            S1avg_cond = squeeze(S1_avg(:,cond,:));
            dimS1all(cond) = m_DimensionalityEstimate_Simple(S1avg_cond');            

        end
        dimension = [dimPMall; dimM1all; dimS1all];
        
        Dim.PM = m_DimensionalityEstimate_Simple(PMavg');
        Dim.M1 = m_DimensionalityEstimate_Simple(M1avg');
        Dim.S1 = m_DimensionalityEstimate_Simple(S1avg');

    saveDir = [dataDir 'DimEstimates' filesep];
    if ~isdir(saveDir)
        mkdir(saveDir)
    end
    save([saveDir Monkey '_' sessions(Day,:)],'dimension','Dim')
    
    clear dim PM M1 S1 A_S1 A_PM B_PM B_S1 B_pred DataAllGaits Neural_cyclic  M1avg PMavg Dim dimension PM_cond M1_cond S1_cond S1avg
    close all
end % The end of days                 (heyoooo!)

a_DimensionalityEstimates_Histograms
% a_Subspace_BetweenTasks
% How closely aligned are neural manifolds in the sensorimotor cortex for a
% variety of locomotor tasks?
% If a cortex generate movement through flexible combinations of stable,
% well-defined neural modes, the neural manifolds corresponding to
% different tasks should be similarly oriented.
% Hypothesis: The manifolds in M1 and S1 will be very closely overlapped
% for locomotor tasks that similar in complexity. In contrast, PM might
% recruit neurons in more task-specific combinations that are not part of
% stable neural modes, so the neural manifolds will not be similarly
% oriented.

% LIMITED AMOUNT OF DIM (according to dim analysis)

% INSTRUCTIONS:
% Run once for each day. The day variable will change each time
% When you're done with all days, run the plot part
clear
close all

rng('shuffle', 'twister') % randomize the seed

method = 1; % 1 for using the PC coefficients, 2 for using the SCORE, 3 using actual EMGs
f_timeBased = 0; % Flag to do time-based rather than time-wrapped analysis
f_norm = 1; % Flag to normalize the correlation by the within-task correlation
f_plotHist = 0; % Flag to also plot histograms for distribution
f_dimAnalysis = 0; % Flag to perform a dimensionality analysis on PCA data
% muscles = [2,4,5,7,8]; % IL, ST, EDL, MG, FHL
muscles = [1:10]; % All
mDim = 15; % Number of dimensions for analysis
mDimMusc = length(muscles);
surrogate_type = 'surrogate-TC'; % TME parameters
nSurrogates = 2;
ctxTitles = {'Premotor';'Motor';'Sensory'};

warning('off','stats:pca:ColRankDefX')

Monkey = 'Elektra';
%%%%%%%%%%% Multiunit
load(['sessions' Monkey '.mat']);
% sessions = ['20190415';'20190425'];
dataDir = ['/Users/ismaelseanez/Google Drive/Working Folder/G-Lab/MONKEY/'...
    'BSI LOCOMOTION Project/Research Material/Analysis/PROCESSED DATA/' Monkey '/MultiUnitClean/'];
dim_directory = ['/Users/ismaelseanez/Google Drive/Working Folder/G-Lab/MONKEY/'...
    'BSI LOCOMOTION Project/Research Material/Analysis/PROCESSED DATA/' Monkey '/MultiUnitClean/DimEstimates/'];

% dataDir = ['C:\Users\physio\Google Drive\BSI LOCOMOTION Project\Research Material\Analysis\PROCESSED DATA\Natalya\MultiUnit\'];
% dim_directory = ['C:\Users\physio\Google Drive\BSI LOCOMOTION Project\Research Material\Analysis\Nicolo\DimensionalityEstimates\'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% Single Unit
% sessions = ['20190425';'20190425'];
% dataDir = ['/Users/ismaelseanez/Google Drive/Working Folder/G-Lab/MONKEY/'...
%     'BSI LOCOMOTION Project/Research Material/Analysis/PROCESSED DATA/Elektra/SingleUnit/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nDays = size(sessions,1);
Results = struct();
all_prinAngles = zeros(nDays,3,5,5,mDim);
all_prinAnglesSurr = zeros(nDays,3,5,5,mDim);
EMG_prinAngles = zeros(nDays,5,5,mDimMusc);
EMG_prinAnglesSurr = zeros(nDays,5,5,mDimMusc);

cond_def = [1:4,5];

% Load important data and run functions
load('ColorMapDarkHot.mat')
load('hope.mat')
colormap_Cyber;
color_ctx = [color_PM(1,:);color_M1(1,:);color_S1(1,:)];

%% for loop for all the days
for Day = 1:nDays
    clear DataAllGaits Neural_cyclic EMG_cyclic
    
    load([ dim_directory Monkey '_' sessions(Day,:)])
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
        
        % Flag outliers
        %     gaitOut1 = find(durations>1.5); % Detect steps longer than 1.5 sec and
        %
        %     [~,gaitOut] = outlieriqr(durations); % Find outlier gaits based on time
        %     [~,gaitOut2] = outlieriqr(FO);
        %     gaitOut3 = find(durations>1.5); % Also brute force and remove longer than 1.5 sec steps
        %     allOut(iCond(gaitOut)) = 1;
        %     allOut(iCond(gaitOut2)) = 1;
        %     allOut(iCond(gaitOut3)) = 1;
        
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
        firingRatesAll(:,cond,:,countCond(cond)) = DataAllGaits(gait).Neural_cyclic';
        %     EMGs(:,cond,:,countCond(cond)) = DataAllGaits(gait).EMG_cyclic';
        countCond(cond) = countCond(cond)+1;
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
        firingRatesAvg(H.S1_badch.channels+96,:,:) = [];
        firingRatesAvg(H.PMM1_badch.channels,:,:) = [];

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
        firingRatesAvg(H.PM_badch.channels+nS1+nM1+nPE,:,:) = [];
        firingRatesAvg(H.S1M1_badch.channels,:,:) = [];

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
    
    
    % Subtract the mean of each condition for each channel
    FR_All =  firingRatesAll; 
    FR_Avg =  firingRatesAvg; 
    FR_All = bsxfun(@minus,FR_All,mean(FR_Avg,3));
    FR_Avg = bsxfun(@minus,FR_Avg,mean(FR_Avg,3));
    
    %% %%%%%%%%%%%% MANIFOLD ANAYSIS %%%%%%%%%%%%%%%%%%%%%%%%%
    
    combos = nchoosek(1:nCond,2); % All possible combinations to test
    minSteps = min(countCond-1); % Minimum number of steps on one condition
    
    for ctx = 1:3 % 1 for Premotor, 2 for Motor, 3 for Sensory
        
        if strcmp(Monkey,'Elektra')
            if ctx == 1
                neurons = 1:nPM;
                shortname = 'PM';
            elseif ctx == 2
                neurons = nPM+1:nPM+nM1;
                shortname = 'M1';
            elseif ctx == 3
                neurons = nPM+nM1+1:nPM+nM1+nS1;
                shortname = 'S1';
            end
            
        else
            
            if ctx == 1
                neurons = nS1+nM1+nPE+1:nS1+nM1+nPE+nPM;
            elseif ctx == 2
                neurons = nS1+1:nS1+nM1;
            elseif ctx == 3
                neurons = 1:nS1;
            end
        end
        
        for comb = 1:size(combos,1)
            
            % GET TASK DATA
            cond1 = combos(comb,1);
            cond2 = combos(comb,2);
            
            % Get data for each conditon
            %
            if method == 1 || method == 3
                X1 = squeeze(FR_All(neurons,cond1,:,1:countCond(cond1)-1));
                X2 = squeeze(FR_All(neurons,cond2,:,1:countCond(cond2)-1));
            elseif method == 2
                X1 = squeeze(FR_All(neurons,cond1,:,1:minSteps));
                X2 = squeeze(FR_All(neurons,cond2,:,1:minSteps));
            end
            % DO PCA FOR ACTUAL DATA AND COMPUTE PRINCIPAL ANGLES
            x1 = reshape(X1,size(X1,1),[]);
            x2 = reshape(X2,size(X2,1),[]);
            
            dim1 = max([2,dimension(ctx,cond_def(cond1))]);
            dim2 = max([2,dimension(ctx,cond_def(cond2))]);
            
            if method == 1 || method == 3
                [ca1,wa1] = pca(x1');
                [ca2,wa2] = pca(x2');
                prinAngles = rad2deg(principal_angles(ca1(:,1:dim1),ca2(:,1:dim2)));
            elseif method == 2
                [ca1,wa1] = pca(x1');
                [ca2,wa2] = pca(x2');
                prinAngles = rad2deg(principal_angles(wa1(:,1:dim1),wa2(:,1:dim2)));
            end
                            
            prinAngles(end+1:mDim)  = nan;
            all_prinAngles(Day,ctx,cond1,cond2,:) = prinAngles;
            all_prinAngles(Day,ctx,cond2,cond1,:) = prinAngles;
            
            fprintf('Day %d cortex %d combination %d of %d \n', Day, ctx, comb, size(combos,1))
        end % End of all possible combinations
        
        
    end % End of for each cortex
   
    
    
end % The end of days                 (heyoooo!)


%% Step 2: Calculate principal angles

prinAngles_Avg = squeeze(mean(all_prinAngles,1,'omitnan'));
all_prinAangles_perDayMean = squeeze(mean(squeeze(mean(all_prinAngles,3,'omitnan')),3,'omitnan'));



%% mybar 
for ii = 1:2
    figure
    pValue = mybar(squeeze(all_prinAangles_perDayMean(:,:,ii)),'SD',[],'unpaired',{'PM','M1','S1'},...
    ['PA component# ' mat2str(ii)],...
    color_ctx);
    set(gca,'TickDir','out')
end

%% Stats
figure
p1_PMM1 = mybar(squeeze(all_prinAangles_perDayMean(:,[2,1],1)),'SD',[],'paired',{'M1','PM'},...
    'Ignore figure (use for stats)',...
    color_ctx)
p1_PMS1 = mybar(squeeze(all_prinAangles_perDayMean(:,[3,1],1)),'SD',[],'paired',{'S1','PM'},...
    'Ignore figure (use for stats)',...
    color_ctx)
p1_M1S1 = mybar(squeeze(all_prinAangles_perDayMean(:,[2,3],1)),'SD',[],'paired',{'M1','S1'},...
    'Ignore figure (use for stats)',...
    color_ctx)
p2_PMM1 = mybar(squeeze(all_prinAangles_perDayMean(:,[2,1],2)),'SD',[],'paired',{'M1','PM'},...
    'Ignore figure (use for stats)',...
    color_ctx)
p2_PMS1 = mybar(squeeze(all_prinAangles_perDayMean(:,[3,1],2)),'SD',[],'paired',{'S1','PM'},...
    'Ignore figure (use for stats)',...
    color_ctx)
p2_M1S1 = mybar(squeeze(all_prinAangles_perDayMean(:,[2,3],2)),'SD',[],'paired',{'M1','S1'},...
    'Ignore figure (use for stats)',...
    color_ctx)
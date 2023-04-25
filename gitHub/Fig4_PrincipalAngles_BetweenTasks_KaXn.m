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

% INSTRUCTIONS:
% Run once for each day. The day variable will change each time
% When you're done with all days, run the plot part
clear
close all

monkey = 'Kara';

load(['sessions' monkey])
% sessions = sessionsCorr;


rng('shuffle', 'twister') % randomize the seed

method = 1; % 1 for using the PC coefficients, 2 for using the SCORE, 3 using actual EMGs
f_timeBased = 0; % Flag to do time-based rather than time-wrapped analysis
f_norm = 1; % Flag to normalize the correlation by the within-task correlation
f_plotHist = 0; % Flag to also plot histograms for distribution
f_dimAnalysis = 0; % Flag to perform a dimensionality analysis on PCA data
% muscles = [2,4,5,7,8]; % IL, ST, EDL, MG, FHL
muscles = [1:10]; % All
mDim = 10; % Number of dimensions for analysis
mDimMusc = length(muscles);
surrogate_type = 'surrogate-TC'; % TME parameters
nSurrogates = 2;
ctxTitles = {'Premotor';'Motor';'Sensory'};

warning('off','stats:pca:ColRankDefX')

%%%%%%%%%%% Multiunit
    dataDir = ['/Users/ismaelseanez/Google Drive/Working Folder/G-Lab/MONKEY/'...
        'BSI LOCOMOTION Project/Research Material/Analysis/PROCESSED DATA/' monkey '/MultiUnitClean/'];
    dim_directory = [dataDir 'DimEstimates/'];
% dataDir = ['C:\Users\physio\Google Drive\BSI LOCOMOTION Project\Research Material\Analysis\PROCESSED DATA\Elektra\MultiUnit\'];
% dim_directory = ['C:\Users\physio\Google Drive\BSI LOCOMOTION Project\Research Material\Analysis\Nicolo\DimensionalityEstimates\'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% Single Unit
% sessions = ['20190425';'20190425'];
% dataDir = ['/Users/ismaelseanez/Google Drive/Working Folder/G-Lab/MONKEY/'...
%     'BSI LOCOMOTION Project/Research Material/Analysis/PROCESSED DATA/Elektra/SingleUnit/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nDays = size(sessions,1);
Results = struct();
all_prinAngles = zeros(nDays,2,2,2,mDim);
all_prinAnglesSurr = zeros(nDays,2,2,2,mDim);
EMG_prinAngles = zeros(nDays,2,2,mDimMusc);
EMG_prinAnglesSurr = zeros(nDays,2,2,mDimMusc);

cond_def = [1,2];

% Load important data and run functions
load('ColorMapDarkHot.mat')
colormap_Cyber;
color_ctx = [color_PM(1,:);color_M1(1,:);color_S1(1,:)];

for Day = 1:nDays

    clear DataAllGaits Neural_cyclic EMG_cyclic
    load([ dim_directory monkey '_' sessions(Day,:)])
    Results(Day).Day = sessions(Day,:);
    S = dir([dataDir filesep monkey '_' sessions(Day,:) '*.mat']);
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
            EMG_cyclic(:,:,j) = D(j).EMG_cyclic;
            DataAllGaits(k).duration = length(D(j).EMG)/100;
            DataAllGaits(k).FOpercent = D(j).RFO/size(D(j).A,1)*100;
            DataAllGaits(k).Neural_cyclic = D(j).Neural_cyclic;
            DataAllGaits(k).EMG_cyclic = D(j).EMG_cyclic;
            DataAllGaits(k).condition = DataAvg(trialIndex).condition;
            DataAllGaits(k).Neural_time = D(j).A;
            DataAllGaits(k).EMG_time = D(j).EMG';
            
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
    [~, idx_Cond,~] = f_findGaitConditions(DataAllGaits,nCond);
    
    % Check if it's time-based approach and correct data
    if f_timeBased == 1
        
        stepLen = zeros(1,numel(DataAllGaits));
        for i = 1:numel(DataAllGaits)
            stepLen(i) = size(DataAllGaits(i).EMG_time,1);
        end
        shortStep = min(stepLen); % shortest step in whole session
        % Cut data to the shortest duration
        for i = 1:numel(DataAllGaits)
            DataAllGaits(i).Neural_cyclic = DataAllGaits(i).Neural_time(1:shortStep,:);
            DataAllGaits(i).EMG_cyclic = DataAllGaits(i).EMG_time(1:shortStep,:);
        end
    end
    
    %% Get average behavior for each condition
    nChan = size(DataAllGaits(1).Neural_cyclic,2);
    lengthNeural = size(DataAllGaits(1).Neural_cyclic,1);
    nEMG = size(DataAllGaits(1).EMG_cyclic,2);
    lengthEMG = size(DataAllGaits(1).EMG_cyclic,1);
    
    % Get the number of gait cycles for each condition
    for i = 1:nCond
        nGaitsCond(i) = sum(idx_Cond(i,:));
    end
    
    maxGaits = max(nGaitsCond); % Maximum number of gaits
    
    % Create firingRate and EMG arrays
    firingRatesAll = nan(nChan,nCond,lengthNeural,maxGaits);
    EMGs = nan(nEMG,nCond,lengthEMG,maxGaits);
    
    % Populate arrays
    countCond = ones(1,nCond);
    for gait = 1:numel(DataAllGaits)
        cond = find(idx_Cond(:,gait)==1);
        firingRatesAll(:,cond,:,countCond(cond)) = DataAllGaits(gait).Neural_cyclic';
        EMGs(:,cond,:,countCond(cond)) = DataAllGaits(gait).EMG_cyclic';
        countCond(cond) = countCond(cond)+1;
    end
    
    % Get the average for each condition
    firingRatesAvg = mean(firingRatesAll,4,'omitnan');
    firingRatesSTE = std(firingRatesAll,0,4,'omitnan')./sqrt(size(firingRatesAll,4));
    EMGsAvg = mean(EMGs,4,'omitnan');
    EMGsSTE = std(EMGs,0,4,'omitnan')./sqrt(size(EMGs,4));
    
    % Get the normalized FR
    firingRatesNorm = bsxfun(@rdivide,firingRatesAvg,max(firingRatesAvg,[],3));
    % firingRatesNorm = firingRatesAvg;
    
    EMGsNorm = bsxfun(@rdivide,EMGsAvg,max(EMGsAvg,[],3));
    
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
    
    % Subtract the mean of each condition for each channel
    FR_All =  firingRatesAll; 
    FR_Avg =  firingRatesAvg; 
    FR_All = bsxfun(@minus,FR_All,mean(FR_Avg,3));
    FR_Avg = bsxfun(@minus,FR_Avg,mean(FR_Avg,3));
    
    %% %%%%%%%%%%%% MANIFOLD ANAYSIS %%%%%%%%%%%%%%%%%%%%%%%%%
    
    combos = nchoosek(1:nCond,2); % All possible combinations to test
    minSteps = min(countCond-1); % Minimum number of steps on one condition
   
    for ctx = 1:2 % 1 for Premotor, 2 for Motor, 3 for Sensory
        
        if ctx == 1
            neurons = 1:nPM;
        elseif ctx == 2
            neurons = nPM+1:nPM+nM1;
        elseif ctx == 3
            neurons = nPM+nM1+1:nPM+nM1+nS1;
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

prinAngles_Avg = squeeze(mean(all_prinAngles,1));
prinAnglesSurr_Avg = squeeze(mean(all_prinAnglesSurr,1));
EMGprinAngles_Avg = squeeze(mean(EMG_prinAngles,1));
EMGprinAnglesSurr_Avg = squeeze(mean(EMG_prinAnglesSurr,1));

all_prinAangles_perDayMean = squeeze(mean(squeeze(mean(all_prinAngles,3,'omitnan')),3,'omitnan'));

prinAnglesSurr_Avg(prinAnglesSurr_Avg==0) = NaN;
sigLine = squeeze(min(prinAnglesSurr_Avg,[],2,'omitnan'));
sigLine = squeeze(min(sigLine,[],2,'omitnan'));


%%

for ii = 1:2
    figure
    pValue = mybar(squeeze(all_prinAangles_perDayMean(:,:,ii)),'SD',[],'unpaired',{'PM','M1'},...
    ['PA component #' mat2str(ii)],...
    color_ctx);
    if ii == 1
        ylim([0 30])
    elseif ii == 2
        ylim([0 60])
    end
    set(gca,'TickDir','out')
    set(gcf,'Position',[461   159   348   391],'Renderer','Painters')
end

figure
p1_PMM1 = mybar(squeeze(all_prinAangles_perDayMean(:,[1,2],1)),'SD',[],'paired',{'M1','PM'},...
    'PA first component',...
    color_ctx)

figure
p2_PMM1 = mybar(squeeze(all_prinAangles_perDayMean(:,[1,2],2)),'SD',[],'paired',{'M1','PM'},...
    'PA first component',...
    color_ctx)


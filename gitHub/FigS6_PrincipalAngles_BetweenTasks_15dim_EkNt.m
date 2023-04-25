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
nSurrogates = 1000;
ctxTitles = {'Premotor';'Motor';'Sensory'};

warning('off','stats:pca:ColRankDefX')

%%%%%%%%%%% Multiunit
monkey = 'Natalya'
load(['sessions' monkey '.mat'])
% sessions = ['20190415';'20190425'];
% dataDir = ['/Users/ismaelseanez/Google Drive/Working Folder/G-Lab/MONKEY/'...
%     'BSI LOCOMOTION Project/Research Material/Analysis/PROCESSED DATA/Elektra/MultiUnitClean/'];
dataDir = ['D:\Google Drive\BSI LOCOMOTION Project\Research Material\Analysis\PROCESSED DATA\Natalya\MultiUnitClean\'];
% dataDir = ['C:\Users\physio\Documents\MATLAB\Neural_PCA\PROCESSED DATA\Elektra\MultiUnit'];

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

% Load important data and run functions
load('ColorMapDarkHot.mat')
colormap_Cyber;
color_ctx = [color_PM(1,:);color_M1(1,:);color_S1(1,:)];

for Day = 1:nDays
    clear DataAllGaits Neural_cyclic EMG_cyclic
    
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
    if strcmp(monkey,'Elektra')
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
    if strcmp(monkey,'Elektra')
       
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
    elseif strcmp(monkey,'Natalya')
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
        
%         nS1 = 32;
%         nM1 = 64;
%         nPE = 32;
%         nPM = 64;

        if strcmp(monkey,'Elektra')
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
            
            if method == 1 || method == 3
                [ca1,wa1] = pca(x1');
                [ca2,wa2] = pca(x2');
                prinAngles = rad2deg(principal_angles(ca1(:,1:mDim),ca2(:,1:mDim)));
            elseif method == 2
                [ca1,wa1] = pca(x1');
                [ca2,wa2] = pca(x2');
                prinAngles = rad2deg(principal_angles(wa1(:,1:mDim),wa2(:,1:mDim)));
            end
            
            % Reshape matrix to: time x unit x steps
            X1 = permute(X1,[2 1 3]);
            X2 = permute(X2,[2 1 3]);
            
            % GET NOISE FEATURES FOR ORIGINAL DATASET
            [sigmaT1, sigmaN1, sigmaC1, M1] = extractFeatures_M(X1);
            [sigmaT2, sigmaN2, sigmaC2, M2] = extractFeatures_M(X2);
            
            % Define surrogate parameters
            params1.margCov{1} = sigmaT1;
            params1.margCov{2} = [];
            params1.margCov{3} = sigmaC1;
            params1.meanTensor = M1.TC;
            
            params2.margCov{1} = sigmaT2;
            params2.margCov{2} = [];
            params2.margCov{3} = sigmaC2;
            params2.meanTensor = M2.TC;
            
            % Fit the maximum entropy distribution
            maxEntropy1 = fitMaxEntropy(params1);
            maxEntropy2 = fitMaxEntropy(params2);
            
            % COMPUTE PRINCIPAL ANGLES BETWEEN DATA SETS
            prinAngles_surr = zeros(nSurrogates,mDim);
            
            for it = 1:nSurrogates
                
                % Generate surrogate datasets
                surr1 = sampleTME(maxEntropy1);
                surr2 = sampleTME(maxEntropy2);
                
                % PCA FOR EACH DATASET
                surr1 = permute(surr1,[2 1 3]);
                surr2 = permute(surr2,[2 1 3]);
                
                surr1 = reshape(surr1,size(surr1,1),[]);
                surr2 = reshape(surr2,size(surr2,1),[]);
                
                if method == 1 || method == 3
                    [c1, w1,~] = pca(surr1');
                    [c2, w2,~] = pca(surr2');
                    prinAngles_surr(it,:) = rad2deg(principal_angles(c1(:,1:mDim),c2(:,1:mDim)));
                elseif method == 2
                    [c1, w1,~] = pca(surr1');
                    [c2, w2,~] = pca(surr2');
                    prinAngles_surr(it,:) = rad2deg(principal_angles(w1(:,1:mDim),w2(:,1:mDim)));
                end
                
            end
            
            % Take the TME surrogate bottom x percentile
            TME_PA = prctile(prinAngles_surr,0.1);
            
            all_prinAngles(Day,ctx,cond1,cond2,:) = prinAngles;
            all_prinAngles(Day,ctx,cond2,cond1,:) = prinAngles;
            all_prinAnglesSurr(Day,ctx,cond1,cond2,:) = TME_PA;
            all_prinAnglesSurr(Day,ctx,cond2,cond1,:) = TME_PA;
            
            fprintf('Day %d cortex %d combination %d of %d \n', Day, ctx, comb, size(combos,1))
        end % End of all possible combinations
        
        
    end % End of for each cortex
   
    
    
end % The end of days                 (heyoooo!)

%% Look at dimensionality of each data acording to ta randomly chosen VAF threshold
% Look for maxDim, to account for at least 80% of the total neural variance
% across all tasks and monkeys.(The results do not really depend on this
% 80% arbitrary choice)
if f_dimAnalysis == 1
    clear VAF VAFmean VAFstd dimCross dimMax
    PCthresh = 80; % threshold for PC
    pcMax = 20;
    for Day = 1:nDays
        for ctx = 1:3
            for cond = 1:nCond
                VAF(Day,ctx,cond,:) = cumsum(Results(Day).Ctx(ctx).Cond(cond).explained(1:pcMax));
            end
        end
    end
    
    VAFmean = squeeze(mean(VAF,1));
    VAFstd = squeeze(std(VAF,1));
    
    % Find the first dimension above the threshold for each ctx and condition
    for ctx = 1:3
        for cond = 1:nCond
            dimCross(ctx,cond) = find(VAFmean(ctx,cond,:)>PCthresh,1,'first');
        end
    end
    dimMax = max(max(dimCross)); % Maximum dimension needed to cross the threshold
    
    figure
    for ctx = 1:3
        if ctx == 1
            color = color_PM(1,:);
        elseif ctx == 2
            color = color_M1(1,:);
        elseif ctx == 3
            color = color_S1(1,:);
        end
        
        for cond = 1:nCond
            subplot(2,nCond,cond);
            shadedErrorBar(1:pcMax,squeeze(VAFmean(ctx,cond,:)),squeeze(VAFstd(ctx,cond,:)),{'color',color});
            hold on
            plot([1 pcMax],[PCthresh PCthresh],'k--')
            plot([dimCross(ctx,cond) dimCross(ctx,cond)],[0 VAFmean(ctx,cond,dimCross(ctx,cond))],...
                'color',color,'LineWidth',1)
            
            if ctx==3
                xlabel('Neural mode')
                h = gca;
                h.TickDir = 'out';
                h.YTick = [0 25 50 75 100];
                h.XTick = [0 pcMax/4 pcMax/4*2 pcMax/4*3 pcMax];
                h.Box = 'off';
                if cond == 1
                    ylabel('Neural VAF (%)')
                end
            end
            
        end
    end
    
% %     % Add EMGs
% %     clear VAF VAFmean VAFstd dimCross
% %     for Day = 1:nDays
% %         for ctx = 4
% %             for cond = 1:nCond
% %                 VAF(Day,cond,:) = cumsum(Results(Day).Ctx(ctx).Cond(cond).explained);
% %             end
% %         end
% %     end
% %     
%     VAFmean = squeeze(mean(VAF,1));
%     VAFstd = squeeze(std(VAF,1));
%     for cond = 1:nCond
%         dimCross(cond) = find(VAFmean(cond,:)>PCthresh,1,'first');
%     end
%     
%     % dimMax = max(max(dimCross)); % Maximum dimension needed to cross the threshold
%     
%     
%     ctx = 4;
%     pcMax = size(VAF,3);
%     for cond = 1:nCond
%         color = color_EMG(2,:);
%         subplot(2,nCond,cond+nCond);
%         shadedErrorBar(1:pcMax,squeeze(VAFmean(cond,:)),squeeze(VAFstd(cond,:)),{'color',color});
%         hold on
%         plot([1 pcMax],[PCthresh PCthresh],'k--')
%         plot([dimCross(cond) dimCross(cond)],[0 VAFmean(cond,dimCross(cond))],...
%             'color',color,'LineWidth',1)
%         
%         xlabel('EMG mode')
%         ylim([30 100])
%         h = gca;
%         h.TickDir = 'out';
%         h.XTick = [0 2 4 6 8 10];
%         h.YTick = [0 25 50 75 100];
%         h.Box = 'off';
%         if cond == 1
%             ylabel('EMG VAF (%)')
%         end
%         
%         
%     end
    
    set(gcf,'Position',[7    351    1273   300],'Renderer','Painters')
    
end

%% Step 2: Calculate principal angles

prinAngles_Avg = squeeze(mean(all_prinAngles,1));
prinAnglesSurr_Avg = squeeze(mean(all_prinAnglesSurr,1));
% EMGprinAngles_Avg = squeeze(mean(EMG_prinAngles,1));
% EMGprinAnglesSurr_Avg = squeeze(mean(EMG_prinAnglesSurr,1));

prinAnglesSurr_Avg(prinAnglesSurr_Avg==0) = NaN;
sigLine = squeeze(min(prinAnglesSurr_Avg,[],2,'omitnan'));
sigLine = squeeze(min(sigLine,[],2,'omitnan'));

%% Plot all principal angle comparisons for each cortex together
colorJet = jet;
figure
for ctx = 1:3
    for comb = 1:size(combos,1)
        
        % GET TASK DATA
        cond1 = combos(comb,1);
        cond2 = combos(comb,2);
        
        subplot(2,3,ctx);
        
        plot(1:mDim,squeeze(prinAngles_Avg(ctx,cond1,cond2,1:mDim)),...
            'LineWidth',1,'color',colorJet(floor(comb*length(colorJet)/size(combos,1)),:))
        hold on
        
    end
    
    % Add minimum significance line for motor cortex (so consistent line)
    plot(1:mDim,sigLine(ctx,:),'--','LineWidth',1,'color',[0.7 0.7 0.7])
    title(ctxTitles{ctx,1})
    if ctx == 1
        ylabel('Principal angle (deg)')
    end
    h = gca;
    h.TickDir = 'out';
    xlabel('Neural mode')
end

% set(gcf,'Position',[ 282   238   796   196],'Renderer','painters')


% Plot histogram of how many dimensions are sifnificant for each comparison

compNames        = cell(1,size(combos,1));
for p = 1:size(combos,1)
    compNames{p} = [condNames{combos(p,1),:} ' vs ' ...
        condNames{combos(p,2),:} ];
end

for ctx = 1:3
    nSig = zeros(size(combos,1),nDays);
    for Day = 1:nDays
        howSig = zeros(size(combos,1),mDim);
        for comb = 1:size(combos,1)
            % GET TASK DATA
            cond1 = combos(comb,1);
            cond2 = combos(comb,2);
            
            a = squeeze(all_prinAngles(Day,ctx,cond1,cond2,:));
            b = sigLine(ctx,:)';
            nSig(comb,Day) = length(find(a<b)); % how many are below significance
            howSig(comb,:) = b-a; % how significant are they
        end
        
        % Do a histogram for every combination
        edges = 0.5:mDim+0.5;
        for comb = 1:size(combos,1)
            [allCounts(comb,:)] = histcounts(nSig(comb,:),edges);
            allCounts(comb,:) = allCounts(comb,:)/sum(allCounts(comb,:)); % Normalize count
        end
    end
    
    allCounts = allCounts/sum(sum(allCounts))*100; % Normalize by total number of counts
    %     [~,sortIdx] = sort(color_combos(:,1)); % Group same colors together
    
    sub(ctx) = subplot(2,3,ctx+3);
    b = bar(allCounts','stack','EdgeColor','none');
    title(ctxTitles{ctx,1})
    colormap jet
    xlim([0 mDim+1])
    %     ylim([0 55])
    set(gca,'Tickdir','out'), box off
    xlabel('Number of similar modes')
    if ctx == 1
        ylabel('Task comparisons (%)')
    end
end
legend(compNames,'Location','NorthWest','FontSize',8), legend boxoff

linkaxes(sub,'y');
set(gcf,'Position',[1  1  1280   704],'Renderer','painters')
% set(gcf,'Position',[1          41        1920         964],'Renderer','painters')
%% Plot principal angles for each combination for each cortex
% figure
% for ctx = 1:3
%     for cond1 = 1:nCond
%         for cond2 = 1:nCond
%             if cond1 ~= cond2
%                 subplot(3,nCond,nCond*ctx-nCond+cond1)
%                 plot(1:mDim,squeeze(prinAngles_Avg(ctx,cond1,cond2,:)),...
%                     'LineWidth',0.5,'color',color_Cyber(cond2,:))
%                 hold on
%                 plot(1:mDim,squeeze(prinAnglesSurr_Avg(ctx,cond1,cond2,:)),...
%                     ':','LineWidth',0.5,'color','k')
%             end
%         end
%
%         ylim([0 90])
%         h = gca;
%         h.TickDir = 'out';
%         h.XTick = [0 mDim/4 mDim*2/4 mDim*3/4 mDim];
%         h.YTick = [0 30  60 90];
%         h.Box = 'off';
%         if cond1 == 1
%             ylabel('Principal angles (deg)')
%
%         end
%         if ctx == 1
%             title(condNames{cond1,:})
%         end
%         if ctx == 3
%             xlabel('Neural mode')
%         end
%
%     end
% end
%
%
% set(gcf,'Position',[ 21  131 1228  477],'Renderer','painters')

%% Plot "effect size" in principal angles for each combination for each cortex
% figure
% k = 1;
% for ctx = 1:3
%     for cond1 = 1:nCond
%         for cond2 = 1:nCond
%             if cond1 ~= cond2
%                 subplot(3,nCond,nCond*ctx-nCond+cond1)
%                 a = squeeze(prinAngles_Avg(ctx,cond1,cond2,:)); % angle between conditions
%                 b = squeeze(prinAnglesSurr_Avg(ctx,cond1,cond2,:)); % significance line
%                 
%                 plot(1:mDim,b-a,...
%                     'LineWidth',0.5,'color',color_Cyber(cond2,:))
%                 hold on
%             end
%         end
%         
%         
%         %         ylim([-30 30])
%         plot([1 mDim],[0 0],'k--')
%         h = gca;
%         h.TickDir = 'out';
%         h.XTick = [0 mDim/4 mDim*2/4 mDim*3/4 mDim];
%         %         h.YTick = [0 30  60 90];
%         h.Box = 'off';
%         if cond1 == 1
%             ylabel('Significant similarity (deg)')
%             
%         end
%         if ctx == 1
%             title(condNames{cond1,:})
%         end
%         if ctx == 3
%             xlabel('Neural mode')
%         end
%         
%     end
% end
% 
% 
% set(gcf,'Position',[ 21  131 1228  477],'Renderer','painters')

%% Plot principal angles for each combination for EMGs
% figure
% k = 1;
% for cond1 = 1:nCond
%     for cond2 = 1:nCond
%         if cond1 ~= cond2
%             subplot(1,nCond,nCond-nCond+cond1)
%             plot(1:mDimMusc,squeeze(EMGprinAngles_Avg(cond1,cond2,:)),...
%                 'LineWidth',0.5,'color',color_Cyber(cond2,:))
%             hold on
%             plot(1:mDimMusc,squeeze(EMGprinAnglesSurr_Avg(cond1,cond2,:)),...
%                 ':','LineWidth',0.5,'color','k')
%         end
%     end
%     
%     %         ylim([0 90])
%     h = gca;
%     h.TickDir = 'out';
%     h.XTick = [0 mDim/4 mDim*2/4 mDim*3/4 mDim];
%     h.YTick = [0 30  60 90];
%     h.Box = 'off';
%     if cond1 == 1
%         ylabel('Principal angles (deg)')
%         
%     end
%     if ctx == 1
%         title(condNames{cond1,:})
%     end
%     if ctx == 3
%         xlabel('EMG mode')
%     end
%     
% end
% 
% 
% set(gcf,'Position',[21 455 1228 153],'Renderer','painters')

%% Plot "effect size" in principal angles for each combination for each cortex
% figure
% k = 1;
% 
% for cond1 = 1:nCond
%     for cond2 = 1:nCond
%         if cond1 ~= cond2
%             subplot(1,nCond,nCond-nCond+cond1)
%             a = squeeze(EMGprinAngles_Avg(cond1,cond2,:)); % angle between conditions
%             b = squeeze(EMGprinAnglesSurr_Avg(cond1,cond2,:)); % significance line
%             
%             plot(1:mDimMusc,b-a,...
%                 'LineWidth',0.5,'color',color_Cyber(cond2,:))
%             hold on
%         end
%     end
%     
%     
%     %         ylim([-30 30])
%     plot([1 mDim],[0 0],'k--')
%     h = gca;
%     h.TickDir = 'out';
%     h.XTick = [0 mDim/4 mDim*2/4 mDim*3/4 mDim];
%     %         h.YTick = [0 30  60 90];
%     h.Box = 'off';
%     if cond1 == 1
%         ylabel('Significant similarity (deg)')
%         
%     end
%     if ctx == 1
%         title(condNames{cond1,:})
%     end
%     if ctx == 3
%         xlabel('EMG mode')
%     end
%     
% end
% 
% 
% 
% set(gcf,'Position',[ 38          46        1228         197],'Renderer','painters')

%% Plot histogram for 3 cortices together
figure
for ctx = [3,2,1]
    nSig = zeros(size(combos,1),nDays);
    for Day = 1:nDays
        howSig = zeros(size(combos,1),mDim);
        for comb = 1:size(combos,1)
            % GET TASK DATA
            cond1 = combos(comb,1);
            cond2 = combos(comb,2);
            
            a = squeeze(all_prinAngles(Day,ctx,cond1,cond2,:));
            b = sigLine(ctx,:)';
%             b = squeeze(all_prinAnglesSurr(Day,ctx,cond1,cond2,:));
            nSig(comb,Day) = length(find(a<b)); % how many are below significance
            howSig(comb,:) = b-a; % how significant are they
        end
    end
    
    %    AllComps(ctx,:) = reshape(nSig,1,[]);
    AllComps(ctx,:) = mean(nSig,2);
    h(ctx) = histogram(AllComps(ctx,:),edges,'Normalization','probability','EdgeColor','none');
    h(ctx).FaceColor = color_ctx(ctx,:);
    h(ctx).FaceAlpha = 0.6;
    hold on
end

figure
pVal = mybar(AllComps','SD',[],'unpaired',{'PM','M1','S1'},...
    'Number of significant differences',...
    color_ctx)
%     b = bar(allCounts','stack','EdgeColor','none');
%     title(ctxTitles{ctx,1})
%     colormap jet
%     xlim([0 mDim+1])
% %     ylim([0 55])
%     set(gca,'Tickdir','out'), box off
%     xlabel('Number of similar modes')
%     if ctx == 1
%         ylabel('Task comparisons (%)')
%     end
%
% legend(compNames,'Location','NorthWest','FontSize',8), legend boxoff
%
% linkaxes(sub,'y');
% set(gcf,'Position',[1  1  1280   704],'Renderer','painters')

%% Plot average and std of each cortex
figure
for ctx = 1:3
    for comb = 1:size(combos,1)
        
        % GET TASK DATA
        cond1 = combos(comb,1);
        cond2 = combos(comb,2);
        
        prinTask(ctx,comb,:) = squeeze(prinAngles_Avg(ctx,cond1,cond2,1:mDim));
        
        plot(squeeze(mean(prinTask(ctx,:,:))),'color',color_ctx(ctx,:))
        hold on
    end
end
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
monkey = 'Kara'
load(['sessions' monkey '7sessions.mat'])
% sessions = ['20190415';'20190425'];
% dataDir = ['/Users/ismaelseanez/Google Drive/Working Folder/G-Lab/MONKEY/'...
%     'BSI LOCOMOTION Project/Research Material/Analysis/PROCESSED DATA/Elektra/MultiUnitClean/'];
dataDir = ['D:\Google Drive\BSI LOCOMOTION Project\Research Material\Analysis\PROCESSED DATA\Kara\MultiUnitClean\'];
% dataDir = ['D:\Google Drive\BSI LOCOMOTION Project\Research Material\Analysis\PROCESSED DATA\Kara\MultiUnit\'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sessions = sessionsCorr; % Run once for this then for Obs
sessions = sessionsObs; 

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
            
            % Run only if conditions are not empty
            if (countCond(cond1)-1 ~= 0) && (countCond(cond2)-1 ~= 0)
                
                
                % Get data for each conditon
                
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
            else
                continue
            end
        end % End of all possible combinations
        
        
    end % End of for each cortex   
    
    
end % The end of days                 (heyoooo!)

%% Remove conditions that weren't even tried
condNone = flip(find(countCond<=1));

for i = 1:length(condNone)
    removCond = condNone(i);
    
    all_prinAngles(:,:,:,removCond,:) = [];
    all_prinAnglesSurr(:,:,:,removCond,:) = [];
    EMG_prinAngles(:,:,removCond,:) = [];
    EMG_prinAnglesSurr(:,:,removCond,:) = [];
    
    all_prinAngles(:,:,removCond,:,:) = [];
    all_prinAnglesSurr(:,:,removCond,:,:) = [];
    EMG_prinAngles(:,removCond,:,:) = [];
    EMG_prinAnglesSurr(:,removCond,:,:) = [];
end

nCond = 2;
combos = nchoosek(1:nCond,2);
%% Step 2: Calculate principal angles

prinAngles_Avg = squeeze(mean(all_prinAngles,1));
prinAnglesSurr_Avg = squeeze(mean(all_prinAnglesSurr,1));
EMGprinAngles_Avg = squeeze(mean(EMG_prinAngles,1));
EMGprinAnglesSurr_Avg = squeeze(mean(EMG_prinAnglesSurr,1));

prinAnglesSurr_Avg(prinAnglesSurr_Avg==0) = NaN;
sigLine = squeeze(min(prinAnglesSurr_Avg,[],2,'omitnan'));
sigLine = squeeze(min(sigLine,[],2,'omitnan'));

%% Plot all principal angle comparisons for each cortex together
colorJet = jet;
figure
for ctx = 1:2
    for comb = 1:size(combos,1)
        
        % GET TASK DATA
        cond1 = combos(comb,1);
        cond2 = combos(comb,2);
        
        % Run only if conditions are not empty
            
            subplot(2,3,ctx);

            plot(1:mDim,squeeze(prinAngles_Avg(ctx,cond1,cond2,1:mDim)),...
                'LineWidth',1,'color',colorJet(floor(comb*length(colorJet)/size(combos,1)),:))
            hold on

    end
    
    % Add minimum significance line for motor cortex (so consistent line)
    plot(1:mDim,sigLine(2,:),'--','LineWidth',1,'color',[0.7 0.7 0.7])
    title(ctxTitles{ctx,1})
    if ctx == 1
        ylabel('Principal angle (deg)')
    end
    ylim([0 90])
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
            b = sigLine(2,:)';
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



%% Plot principal angles for each combination for EMGs
figure
k = 1;
for cond1 = 1:nCond
    for cond2 = 1:nCond
        if cond1 ~= cond2
            subplot(1,nCond,nCond-nCond+cond1)
            plot(1:mDimMusc,squeeze(EMGprinAngles_Avg(cond1,cond2,:)),...
                'LineWidth',0.5,'color',color_Cyber(cond2,:))
            hold on
            plot(1:mDimMusc,squeeze(EMGprinAnglesSurr_Avg(cond1,cond2,:)),...
                ':','LineWidth',0.5,'color','k')
        end
    end
    
    %         ylim([0 90])
    h = gca;
    h.TickDir = 'out';
    h.XTick = [0 mDim/4 mDim*2/4 mDim*3/4 mDim];
    h.YTick = [0 30  60 90];
    h.Box = 'off';
    if cond1 == 1
        ylabel('Principal angles (deg)')
        
    end
    if ctx == 1
        title(condNames{cond1,:})
    end
    if ctx == 3
        xlabel('EMG mode')
    end
    
end


set(gcf,'Position',[21 455 1228 153],'Renderer','painters')

%% Plot "effect size" in principal angles for each combination for each cortex
figure
k = 1;

for cond1 = 1:nCond
    for cond2 = 1:nCond
        if cond1 ~= cond2
            subplot(1,nCond,nCond-nCond+cond1)
            a = squeeze(EMGprinAngles_Avg(cond1,cond2,:)); % angle between conditions
            b = squeeze(EMGprinAnglesSurr_Avg(cond1,cond2,:)); % significance line
            
            plot(1:mDimMusc,b-a,...
                'LineWidth',0.5,'color',color_Cyber(cond2,:))
            hold on
        end
    end
    
    
    %         ylim([-30 30])
    plot([1 mDim],[0 0],'k--')
    h = gca;
    h.TickDir = 'out';
    h.XTick = [0 mDim/4 mDim*2/4 mDim*3/4 mDim];
    %         h.YTick = [0 30  60 90];
    h.Box = 'off';
    if cond1 == 1
        ylabel('Significant similarity (deg)')
        
    end
    if ctx == 1
        title(condNames{cond1,:})
    end
    if ctx == 3
        xlabel('EMG mode')
    end
    
end



set(gcf,'Position',[ 38          46        1228         197],'Renderer','painters')

%% Plot histogram for 2 cortices together
figure
for ctx = 1:2
    nSig = zeros(size(combos,1),nDays);
    for Day = 1:nDays
        howSig = zeros(size(combos,1),mDim);
        for comb = 1:size(combos,1)
            % GET TASK DATA
            cond1 = combos(comb,1);
            cond2 = combos(comb,2);
            
            a = squeeze(all_prinAngles(Day,ctx,cond1,cond2,:));
            b = squeeze(all_prinAnglesSurr(Day,ctx,cond1,cond2,:));
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
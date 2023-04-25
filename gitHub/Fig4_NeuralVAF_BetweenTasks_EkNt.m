% a_NeuralVAF_BetweenTasks
% How closely aligned are neural manifolds in the sensorimotor cortex for a
% variety of locomotor tasks?
% This code will project the data from one task onto the PC manifold of
% another task, then compute the neural variance accounted for by this new
% manifold.
% Finally, it will compare the VAF of that manifold, to the VAF if the
% same-task data had been projected onto its original 12-D manifold

clear
close all
%%
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
% load('sessionsElektra.mat')
% dataDir = ['/Users/ismaelseanez/Google Drive/Working Folder/G-Lab/MONKEY/'...
%     'BSI LOCOMOTION Project/Research Material/Analysis/PROCESSED DATA/Elektra/MultiUnit/'];
% dataDir = ['D:\Google Drive\BSI LOCOMOTION Project\Research Material\Analysis\PROCESSED DATA\Elektra\MultiUnit\'];
% 

monkey = 'Elektra'
load(['sessions' monkey '.mat'])
% dataDir = ['/Users/ismaelseanez/Google Drive/Working Folder/G-Lab/MONKEY/'...
%     'BSI LOCOMOTION Project/Research Material/Analysis/PROCESSED DATA/Natalya/MultiUnit/'];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% Single Unit
% sessions = ['20190425';'20190425'];
% dataDir = ['/Users/ismaelseanez/Google Drive/Working Folder/G-Lab/MONKEY/'...
%     'BSI LOCOMOTION Project/Research Material/Analysis/PROCESSED DATA/Elektra/SingleUnit/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%Dimensionality analysis for VAF
% dataDir = ['C:\Users\physio\Google Drive\BSI LOCOMOTION Project\'...
%     'Research Material\Analysis\PROCESSED DATA\Natalya\Multiunit\'];
% dim_directory = ['C:\Users\physio\Google Drive\BSI LOCOMOTION Project\Research Material\'...
%     'Analysis\Nicolo\DimensionalityEstimates\'];

% dataDir = ['/Users/ismaelseanez/Google Drive/Working Folder/G-Lab/MONKEY/'...
%     'BSI LOCOMOTION Project/Research Material/Analysis/PROCESSED DATA/' monkey '/MultiUnitClean/'];


dataDir = ['/Users/Borgo/Documents/MATLAB/fri_mks/BSI LOCOMOTION Project_end_of_july_2021/Research Material/Analysis/PROCESSED DATA/' monkey '/MultiunitClean/'];

dim_directory = ['/Users/Borgo/Documents/MATLAB/fri_mks/BSI LOCOMOTION Project_end_of_july_2021/Research Material/Analysis/PROCESSED DATA/' monkey '/MultiUnitClean/DimEstimates/'];

% dataDir = ['C:\Users\physio\Google Drive\Working Folder\G-Lab\MONKEY\BSI LOCOMOTION Project\Research Material\Analysis\PROCESSED DATA\Elektra\MultiUnit\'];
% dim_directory = ['C:\Users\physio\Google Drive\Working Folder\G-Lab\MONKEY\BSI LOCOMOTION Project\Research Material\Analysis\PROCESSED DATA\Elektra\MultiUnit\DimEstimates\'];

cond_def = [1:4,5];


nDays = size(sessions,1);
Results = struct();
normVAF = zeros(nDays,3,10,2);


% Load important data and run functions
load('ColorMapDarkHot.mat')
colormap_Cyber;
color_ctx = [color_PM(1,:);color_M1(1,:);color_S1(1,:)];

for Day = 1:nDays

    load([ dim_directory monkey '_' sessions(Day,:)])
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
        DataAllGaits(k).duration = length(D(j).times)/100;
        DataAllGaits(k).FOpercent = D(j).RFO/size(D(j).A,1)*100;
        DataAllGaits(k).Neural_cyclic = D(j).Neural_cyclic;
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
    %     nEMG = size(DataAllGaits(1).EMG_cyclic,2);
    %     lengthEMG = size(DataAllGaits(1).EMG_cyclic,1);
    
    % Get the number of gait cycles for each condition
    for i = 1:nCond
        nGaitsCond(i) = sum(idx_Cond(i,:));
    end
    
    maxGaits = max(nGaitsCond); % Maximum number of gaits
    
    % Create firingRate and EMG arrays
    firingRatesAll = nan(nChan,nCond,lengthNeural,maxGaits);
    %     EMGs = nan(nEMG,nCond,lengthEMG,maxGaits);
    
    % Populate arrays
    countCond = ones(1,nCond);
    for gait = 1:numel(DataAllGaits)
        cond = find(idx_Cond(:,gait)==1);
        firingRatesAll(:,cond,:,countCond(cond)) = DataAllGaits(gait).Neural_cyclic';
        %         EMGs(:,cond,:,countCond(cond)) = DataAllGaits(gait).EMG_cyclic';
        countCond(cond) = countCond(cond)+1;
    end
    
    % Get the average for each condition
    firingRatesAvg = mean(firingRatesAll,4,'omitnan');
    firingRatesSTE = std(firingRatesAll,0,4,'omitnan')./sqrt(size(firingRatesAll,4));
    %     EMGsAvg = mean(EMGs,4,'omitnan');
    %     EMGsSTE = std(EMGs,0,4,'omitnan')./sqrt(size(EMGs,4));
    
    % Get the normalized FR
    firingRatesNorm = bsxfun(@rdivide,firingRatesAvg,max(firingRatesAvg,[],3));
    % firingRatesNorm = firingRatesAvg;
    
    %     EMGsNorm = bsxfun(@rdivide,EMGsAvg,max(EMGsAvg,[],3));
    
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
%     EMG_All = bsxfun(@minus,EMGs,mean(EMGs,3,'omitnan'));
%     EMG_Avg = bsxfun(@minus,EMGsAvg,mean(EMGsAvg,3,'omitnan'));
%     
    %% %%%%%%%%%%%% MANIFOLD ANAYSIS %%%%%%%%%%%%%%%%%%%%%%%%%
    
    combos = nchoosek(1:nCond,2); % All possible combinations to test
    minSteps = min(countCond-1); % Minimum number of steps on one condition
   
    for ctx = 1:3 % 1 for Premotor, 2 for Motor, 3 for Sensory
        
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
        
        for cond1 = 1:nCond
            
            % Source activity
            X1 = squeeze(FR_All(neurons,cond1,:,1:countCond(cond1)-1));
            x1 = reshape(X1,size(X1,1),[]);
            dim1 = max([2,dimension(ctx,cond_def(cond1))]);
            
            for cond2 = 1:nCond
                
                dim2 = max([2,dimension(ctx,cond_def(cond2))]);
                
                % Target Encoding Matrix
                X2 = squeeze(FR_All(neurons,cond2,:,1:countCond(cond2)-1));

                x2 = reshape(X2,size(X2,1),[]);
                
                [C2,W2] = pca(x2');
                
                % Compute VAF by source activity on target encoding
                min2dims = min([dim1,dim2]);
                VAF(Day,ctx,cond1,cond2) = a_Cross_Vaf_Computation(x1,C2(:,1:min2dims)');
            end
        end
                        
    end % End of for each cortex
   

    
end % The end of days                 (heyoooo!)


%% normalize according to intratask 
if f_norm == 1
for i = 1:nCond
    VAFnorm(:,:,i,:) = VAF(:,:,i,:) ./ VAF(:,:,i,i);
end

% Mean intra-task VAF for each cortex 
meanVAF = squeeze(mean(VAF,1));
for ctx = 1:3
    meanVAFctx(ctx) = mean(diag(squeeze(meanVAF(ctx,:,:))));
end
else
     VAFnorm(:,:,i,:) = VAF(:,:,i,:);
    
end
colormap_Cyber
%% plots histogram then we will see
figure
cortex_data = {};
for ctx = 1:3
    Gen = [];
    for cond1 = 1:nCond
        for cond2 = 1:nCond
            if cond1 ~=cond2
                Gen = [Gen; VAFnorm(:,ctx,cond1,cond2)];
            end             
        end
    end
    
    
%     sub(ctx) = subplot(1,3,ctx);
    nb = 0:0.1:1;
    %nb = calcnbins(Gen, 'middle');
    h1 = histogram(Gen,nb,'Normalization','probability');
    h1.FaceColor = color_ctx(ctx,:);
    h1.EdgeColor = 'none';
    h1.FaceAlpha = 0.6;
    hold on
    
    axis('square')
    
   % plot mean and std
   yMean = 0.8;
    plot(mean(Gen,'omitnan'),yMean,'o','MarkerFaceColor',color_ctx(ctx,:),'MarkerEdgeColor','none')
    hold on
    err = std(reshape(Gen,1,[]),'omitnan');
    plot([mean(Gen,'omitnan')-err mean(Gen,'omitnan')+err],[yMean yMean],'-','LineWidth',1,'color',color_ctx(ctx,:))
    
    ylabel('Task comparisons (%)')
    xlim([0 1])
    ylim([0 1])
    xticks([0 .25 .5 .75 1])
    xticks([0 .25 .5 .75 1])
    set(gca,'TickDir','out');
    
end

%% Compute mean generalization performance of each day
normVAF = zeros(nDays,3);
for ctx = 1:3
    for Day = 1:nDays
        Gen = [];
        for cond1 = 1:5
            for cond2 = 1:5
                if cond1 ~=cond2
                    Gen = [Gen; VAFnorm(Day,ctx,cond1,cond2)];
                end             
            end
        end
        normVAF(Day,ctx) = mean(Gen);
    end
end

color = [color_PM(1,:); color_M1(1,:); color_S1(1,:)];

figure
[p] = mybar(normVAF,'SD',[0 1],'unpaired',{'PM','M1','S1'},'Neural VAF',color)
set(gca,'TickDir','out');

%% stats

figure
[p_PMM1] = mybar(normVAF(:,[1,2]),'SD',[0 1],'paired',{'PM','M1'},'Task-Independent',color)
set(gca,'TickDir','out');
[p_M1S1] = mybar(normVAF(:,[2,3]),'SD',[0 1],'paired',{'M1','S1'},'Task-Independent',color)
set(gca,'TickDir','out');
[p_PMS1] = mybar(normVAF(:,[1,3]),'SD',[0 1],'paired',{'PM','S1'},'Task-Independent',color)


avgVAF = mean(normVAF)
stdVAF = std(normVAF)


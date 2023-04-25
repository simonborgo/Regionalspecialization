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

monkey = 'Kara';

load(['sessions' monkey])

f_timeBased = 0; % Flag to do time-based rather than time-wrapped analysis
f_norm = 1; % Flag to normalize the correlation by the within-task correlation
f_plotHist = 0; % Flag to also plot histograms for distribution
f_dimAnalysis = 0; % Flag to perform a dimensionality analysis on PCA data
% muscles = [2,4,5,7,8]; % IL, ST, EDL, MG, FHL
muscles = [1:10]; % All
mDim = 10; % Number of dimensions for analysis
mDimMusc = length(muscles);
ctxTitles = {'Premotor';'Motor';'Sensory'};

warning('off','stats:pca:ColRankDefX')

%%%%%%%%%%% Multiunit
% load('sessionsElektra.mat')
% dataDir = ['/Users/ismaelseanez/Google Drive/Working Folder/G-Lab/MONKEY/'...
%     'BSI LOCOMOTION Project/Research Material/Analysis/PROCESSED DATA/Elektra/MultiUnit/'];
% dataDir = ['D:\Google Drive\BSI LOCOMOTION Project\Research Material\Analysis\PROCESSED DATA\Elektra\MultiUnit\'];
% 

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

dataDir = ['C:\Users\physio\Google Drive\BSI LOCOMOTION Project\Research Material\Analysis\PROCESSED DATA\Kara\MultiUnitClean\'];

dim_directory = ['C:\Users\physio\Google Drive\BSI LOCOMOTION Project\Research Material\Analysis\PROCESSED DATA\Kara\MultiUnitClean\DimEstimates\'];
% dim_directory = ['/Users/ismaelseanez/Google Drive/Working Folder/G-Lab/MONKEY/'...
%     'BSI LOCOMOTION Project/Research Material/Analysis/PROCESSED DATA/' monkey '/MultiUnitClean/DimEstimates/'];


cond_def = [1,2];


% dataDir = ['C:\Users\physio\Google Drive\BSI LOCOMOTION Project\Research Material\Analysis\PROCESSED DATA\' ...
%         monkey '\Multiunit\'];
% dim_directory = ['C:\Users\physio\Google Drive\BSI LOCOMOTION Project\Research Material\Analysis\PROCESSED DATA\' ...
%         monkey '\Multiunit\DimEstimates\'];
    
nDays = size(sessions,1);
Results = struct();
VAF = nan(nDays,2,4,4); % days, ctx, nCond, nCond

condActual = zeros(nDays,2);
for i = 1:nDays
    if dayType(i) == 1
        condActual(i,:) = [1,2];
    elseif dayType(i) == 2
        condActual(i,:) = [3,4];
    end
end

% Load important data and run functions
load('ColorMapDarkHot.mat')
colormap_Cyber;
color_ctx = [color_PM(1,:);color_M1(1,:);color_S1(1,:)];

for Day = 1:nDays

    load([ dim_directory monkey '_' sessions(Day,:)])
    clear DataAllGaits Neural_cyclic EMG_cyclic DataAvg
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
%     
    %% %%%%%%%%%%%% MANIFOLD ANAYSIS %%%%%%%%%%%%%%%%%%%%%%%%%
    
    combos = nchoosek(1:nCond,2); % All possible combinations to test
    minSteps = min(countCond-1); % Minimum number of steps on one condition
   
    for ctx = 1:2 % 1 for Premotor, 2 for Motor, 3 for Sensory
        
        
            if ctx == 1
                neurons = 1:nPM;
                shortname = 'PM';
            elseif ctx == 2
                neurons = nPM+1:nPM+nM1;
                shortname = 'M1';
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
nCond = 4;
for i = 1:nCond
    VAFnorm(:,:,i,:) = VAF(:,:,i,:) ./ VAF(:,:,i,i);
end

% Mean intra-task VAF for each cortex 
meanVAF = squeeze(mean(VAF,1,'omitnan'));
for ctx = 1:2
    meanVAFctx(ctx) = mean(diag(squeeze(meanVAF(ctx,:,:))),'omitnan');
end

colormap_Cyber
%% plots histogram then we will see
figure
cortex_data = {};
for ctx = 1:2
    Gen = [];
    for cond1 = 1:nCond
        for cond2 = 1:nCond
            if cond1 ~=cond2
                Gen = [Gen; VAFnorm(:,ctx,cond1,cond2)];
            end             
        end
    end
    
    Gen(isnan(Gen)) = [];
    sub(ctx) = subplot(1,3,ctx);
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
normVAF = zeros(nDays,2);
for ctx = 1:2
    for Day = 1:nDays
        Gen = [];
        for cond1 = 1:nCond
            for cond2 = 1:nCond
                if cond1 ~=cond2
                    Gen = [Gen; VAFnorm(Day,ctx,cond1,cond2)];
                end             
            end
        end
        normVAF(Day,ctx) = mean(Gen,'omitnan');
    end
end

color = [color_PM(1,:); color_M1(1,:); color_S1(1,:)];

figure
[p] = mybar(normVAF,'SD',[0 1],'unpaired',{'PM','M1'},'VAF across tasks',color);
ylim([0 1])
set(gca,'TickDir','out');
set(gcf,'Position',[461   159   348   391],'Renderer','Painters')
%% stats

figure
[p_PMM1] = mybar(normVAF(:,[1,2]),'SD',[0 1],'paired',{'PM','M1'},'VAF across tasks',color)
set(gca,'TickDir','out');
set(gcf,'Position',[461   159   348   391],'Renderer','Painters')


avgVAF = mean(normVAF)
stdVAF = std(normVAF)


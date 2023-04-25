% a_Corr_BetweenTasks
% This script will analyze the correlation of individual
% channels/neurons/emgs across different tasks
% INSTRUCTIONS:
% Run once for each day. The day variable will change each time
% When you're done with all days, run the plot part
clear
% close all
f_stat_between_each_combinaition = 1;
f_norm = 0; % Flag to normalize the correlation by the within-task correlation
f_kin = 0;
f_singleUnit = 1;
f_plotHist = 0; % Flag to also plot histograms for distribution
muscles = [2,7,5,8]; % IL, EDL, MG, FHL

%%%%%%%%%%% Multiunit
% load('sessionsElektra.mat')
% dataDir = ['/Users/ismaelseanez/Google Drive/Working Folder/G-Lab/MONKEY/'...
%     'BSI LOCOMOTION Project/Research Material/Analysis/PROCESSED DATA/Elektra/MultiUnit/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% Multiunit with Kinematics
% sessions = ['20190415';'20190425';'20190426'];
% load('GoodKinematics.mat')
% dataDir = ['/Users/ismaelseanez/Google Drive/Working Folder/G-Lab/MONKEY/'...
%     'BSI LOCOMOTION Project/Research Material/Analysis/PROCESSED DATA/Elektra/MultiunitWithKin/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% Single Unit
% sessions = ['20190425';'20190426';'20190430';'20190506';'20190514';'20190522';'20190531';'20190604'];
% dataDir = ['/Users/ismaelseanez/Google Drive/Working Folder/G-Lab/MONKEY/'...
%     'BSI LOCOMOTION Project/Research Material/Analysis/PROCESSED DATA/Elektra/SingleUnit/'];
% mk = 1;
% f_singleUnit = 1;
% f_kin = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% Monkey specification %%%%%%%%%%%%%%%%
Monkey = 'Elektra';
load(['sessions' Monkey '.mat']);
% dataDir = ['/Users/ismaelseanez/Google Drive/Working Folder/G-Lab/MONKEY/'...
%     'BSI LOCOMOTION Project/Research Material/Analysis/PROCESSED DATA/' Monkey '/SingleUnitClean/'];


%dataDir = ['/Users/Borgo/Google Drive/BSI LOCOMOTION Project/Research Material/Analysis/PROCESSED DATA/' Monkey '/SingleUnitClean/'];
 dataDir = ['/Users/Borgo/Documents/MATLAB/fri_mks/BSI LOCOMOTION Project_end_of_july_2021/Research Material/Analysis/PROCESSED DATA/' Monkey '/SingleUnitClean/'];

if strcmp(Monkey,'Elektra')
    mk = 1;
elseif strcmp(Monkey,'Natalya')
    mk = 2;
end
f_singleUnit = 1;
f_kin = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nDays = size(sessions,1);
Results = struct();

for Day = 1:nDays
    load('ColorMapDarkHot.mat')
    colormap_Cyber;
    clear DataAllGaits Neural_cyclic EMG_cyclic Kin_cyclic
    

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
        
        Results(Day).avgChan = mean(Data.SpikeRate,1);
        for j = 1:numel(D)
            
            Neural_cyclic(:,:,j) = D(j).Neural_cyclic;
%             EMG_cyclic(:,:,j) = D(j).EMG_cyclic;
             DataAllGaits(k).duration = length(D(j).times)/100;
            DataAllGaits(k).FOpercent = D(j).RFO/size(D(j).A,1)*100;
            DataAllGaits(k).Neural_cyclic = D(j).Neural_cyclic;
%             DataAllGaits(k).EMG_cyclic = D(j).EMG_cyclic;
            DataAllGaits(k).condition = DataAvg(trialIndex).condition;
            
            if f_kin == 1
                Kin_cyclic(:,:,j) = D(j).Kin_cyclic;
                DataAllGaits(k).Kin_cyclic = D(j).Kin_cyclic;
            end
            
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
%     
    % Get the normalized FR
    firingRatesNorm = bsxfun(@rdivide,firingRatesAvg,max(firingRatesAvg,[],3));
    % firingRatesNorm = firingRatesAvg;
    
%     EMGsNorm = bsxfun(@rdivide,EMGsAvg,max(EMGsAvg,[],3));
    
    % Get the peak location for each condition for sorting
    if nChan == 160 % Number of channels
        nPM = 48; % Number of neurons on each cortex
        nM1 = 48;
        nS1 = 64;
    elseif nChan <160
        nPM = Data.nCellPM; % Number of neurons on each cortex
        nM1 = Data.nCellM1;
        nS1 = Data.nCellS1;
    end
    
    % Process kinematics
    if f_kin == 1
        nKin = size(DataAllGaits(1).Kin_cyclic,2);
        lengthKin = size(DataAllGaits(1).Kin_cyclic,1);
        KinAll = nan(nKin,nCond,lengthKin,maxGaits);
        % Populate arrays
        countCond = ones(1,nCond);
        for gait = 1:numel(DataAllGaits)
            cond = find(idx_Cond(:,gait)==1);
            KinAll(:,cond,:,countCond(cond)) = DataAllGaits(gait).Kin_cyclic';
            countCond(cond) = countCond(cond)+1;
        end
        
        % Smooth the kinematics with Savitzky-Golay filter
        for i = 1:size(KinAll,1)
            for cond = 1:nCond
                for j = 1:countCond(cond)-1
                    a = squeeze(KinAll(i,cond,:,j));
                    a2 = zeros(300,1);
                    a2(1:100) = a(1); % Add the first and last values to the ends
                    a2(101:200) = a;
                    a2(201:300) = a(end);
                    
                    b = smooth(a2-mean(a,'omitnan'),15,'sgolay');
                    KinAll(i,cond,:,j) = b(101:200);
                end
            end
        end
        
        % Remove trials based on the outlier on the means of the hip angle
        for cond = 1:nCond
            av = zeros(1,countCond(cond)-1);
            numNan = zeros(1,countCond(cond)-1);
            for j = 1:countCond(cond)-1
                av(j) = mean(KinAll(1,cond,:,j),'omitnan');
                numNan(j) = length(find(isnan(KinAll(1,cond,:,j))));
            end
            [~,id2remove] = outlieriqr(av(~isnan(av)));
            id2remove2 = find(numNan>20); % Remove the ones that are more than 30% incomplete
            id2removeAll = unique([id2remove id2remove2]);
            
            % Remove them from kinematics, EMG, and FR
            KinAll(:,cond,:,id2removeAll) = NaN;
            EMG_All(:,cond,:,id2removeAll) = NaN;
            FR_All(:,cond,:,id2removeAll) = NaN;
            countCond(cond) = countCond(cond) - length(id2removeAll); % Tell countCond that you removed trials
        end
        
        Kin_Avg = trimmean(KinAll,5,4); % mean without outliers in time profile (omits nans)
        Kin_STE = std(KinAll,0,4,'omitnan')./sqrt(size(KinAll,4));
        Kin_Norm = bsxfun(@rdivide,Kin_Avg,max(Kin_Avg,[],3));
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
    
    % Do not demean the data
    FR_All = firingRatesAll; % bsxfun(@minus,firingRatesAll,mean(firingRatesAll,3,'omitnan'));
    FR_Avg = firingRatesAvg; % bsxfun(@minus,firingRatesAvg,mean(firingRatesAvg,3,'omitnan'));
%     EMG_All = EMGs; % bsxfun(@minus,EMGs,mean(EMGs,3,'omitnan'));
%     EMG_Avg = EMGsAvg; % bsxfun(@minus,EMGsAvg,mean(EMGsAvg,3,'omitnan'));
    
    %% Calculate correlation between average firing rate and modulation depth between tasks
    meanFR = squeeze(mean(FR_All,3,'omitnan'));
    modDepth = squeeze(bsxfun(@minus,max(FR_All,[],3,'omitnan'),min(FR_All,[],3,'omitnan')));
    maxFR = squeeze(max(FR_All,[],3,'omitnan'));
    
    if f_singleUnit == 1
        % Save the average result for each neuron
        Results(Day).meanFR = mean(meanFR,3,'omitnan');
        Results(Day).modDepth = mean(modDepth,3,'omitnan');
        Results(Day).maxFR = mean(maxFR,3,'omitnan');
        Results(Day).nChans = [nPM,nM1,nS1];
        
    end
end % End of days                 (heyoooo!)

%% Concatenate all neurons per cortex
% Get correlations for all days of cortex channels into one matrix
meanFR_PM = [];
modDepth_PM = [];
maxFR_PM = [];
meanFR_M1 = [];
modDepth_M1 = [];
maxFR_M1 = [];
meanFR_S1 = [];
modDepth_S1 = [];
maxFR_S1 = [];

for ctx = 1:3 % 1 for Premotor, 2 for Motor, 3 for Sensory
    
    for Day = 1:nDays
        
        nPM = Results(Day).nChans(1);
        nM1 = Results(Day).nChans(2);
        nS1 = Results(Day).nChans(3);
        
        if ctx == 1
            chans = 1:nPM;
            color = color_PM(1,:);
            colorConf = color_PMw;
            
            meanFR_PM = [ meanFR_PM; Results(Day).meanFR(chans,:)];
            modDepth_PM = [ modDepth_PM; Results(Day).modDepth(chans,:)];
            maxFR_PM = [ maxFR_PM; Results(Day).maxFR(chans,:)];
            
        elseif ctx == 2
            chans = nPM+1:nPM+nM1;
            color = color_M1(1,:);
            colorConf = color_M1w;
            
            meanFR_M1 = [ meanFR_M1; Results(Day).meanFR(chans,:)];
            modDepth_M1 = [ modDepth_M1; Results(Day).modDepth(chans,:)];
            maxFR_M1 = [ maxFR_M1; Results(Day).maxFR(chans,:)];
            
        elseif ctx == 3
            chans = nPM+nM1+1:nPM+nM1+nS1;
            color = color_S1(1,:);
            colorConf = color_S1w;
            
            meanFR_S1 = [ meanFR_S1; Results(Day).meanFR(chans,:)];
            modDepth_S1 = [ modDepth_S1; Results(Day).modDepth(chans,:)];
            maxFR_S1 = [ maxFR_S1; Results(Day).maxFR(chans,:)];
        end
        
    end
end

% And then put all cortices together... to make rest of code easier
meanFRall = [meanFR_PM; meanFR_M1; meanFR_S1];
modDepthall = [modDepth_PM; modDepth_M1; modDepth_S1];
maxFRall = [maxFR_PM; maxFR_M1; maxFR_S1];
nPM = size(meanFR_PM,1);
nM1 = size(meanFR_M1,1);
nS1 = size(meanFR_S1,1);

%% Scatterplot for each
for  ctx = 1:3
    figure(ctx)
    if ctx == 1
        chans = 1:nPM;
        color = color_PM(1,:);
        colorConf = color_PMw;
    elseif ctx == 2
        chans = nPM+1:nPM+nM1;
        color = color_M1(1,:);
        colorConf = color_M1w;
    elseif ctx == 3
        chans = nPM+nM1+1:nPM+nM1+nS1;
        color = color_S1(1,:);
        colorConf = color_S1w;
    end
    k = 1;
    for cond2 = 1:nCond
        for cond1 = 1:nCond
            subplot(nCond,nCond,k)
            x = meanFRall(chans,cond1); 
            y = meanFRall(chans,cond2);
%             x = modDepthall(chans,cond1);
%             y = modDepthall(chans,cond2);
            if mk == 1 % Elektra
                scatter(x,y,5,'filled','o','MarkerFaceAlpha',.7,'MarkerEdgeAlpha',.7,'MarkerFaceColor',[0.5 0.5 0.5]','MarkerEdgeColor','none');
            elseif mk == 2 % Natalya
                scatter(x,y,5,'filled','o','MarkerFaceAlpha',.7,'MarkerEdgeAlpha',.7,'MarkerFaceColor',color,'MarkerEdgeColor','none');
            end
            m = x\y;
            x1 = 1:150;
            y1 = m*x1;
            hold on
            if mk == 1 % Elektra
                plot(x1,y1,'color',[0.5 0.5 0.5])
            elseif mk == 2 % Natalya
                plot(x1,y1,'color',color)
            end
            plot([0 150],[0 150],'--','color','k')
            k = k+1;
            ylim([0 100])
            xlim([0 100])
            set(gca,'TickDir','out')
            
            % Add regresion info
            yCalc1 = m*x;
            Rsq1 = 1 - sum((y - yCalc1).^2)/sum((y - mean(y)).^2);
            if mk == 1
                text(10,80,mat2str(round(Rsq1,1)),'Color',[0.5 0.5 0.5])
                text(10,120,mat2str(round(m,1)),'Color',[0.5 0.5 0.5])
            elseif mk == 2
                text(75,20,mat2str(round(Rsq1,1)),'Color',color)
                text(75,50,mat2str(round(m,1)),'Color',color)
            end
            
            xlabel(condNames{cond1,:})
            ylabel(condNames{cond2,:})
            
        end
    end
    set(gcf,'Position',[176     6   888   689],'Renderer','Painters')
end


%% Plot average firing rate and modulation depth for each cortex and each condition
% figure
% for  ctx = 1:3
%     if ctx == 1
%         chans = 1:nPM;
%         color = color_PM(1,:);
%         colorConf = color_PMw;
%     elseif ctx == 2
%         chans = nPM+1:nPM+nM1;
%         color = color_M1(1,:);
%         colorConf = color_M1w;
%     elseif ctx == 3
%         chans = nPM+nM1+1:nPM+nM1+nS1;
%         color = color_S1(1,:);
%         colorConf = color_S1w;
%     end
% 
%     subplot(1,3,ctx)
%     mybar(modDepthall(chans,:),'SD',[],'unpaired',condNames','ModDepth',color_Cyber(1:nCond,:));
%     ylim([0 80])
%     set(gca,'TickDir','out')
% end
%         
% set(gcf,'Position',[143   348   848   213],'Renderer','Painters')

%% Plot histogram for firing rate and modulation depth for each cortex and each condition
% figure
% for  ctx = 1:3
%     if ctx == 1
%         chans = 1:nPM;
%         color = color_PM(1,:);
%         colorConf = color_PMw;
%     elseif ctx == 2
%         chans = nPM+1:nPM+nM1;
%         color = color_M1(1,:);
%         colorConf = color_M1w;
%     elseif ctx == 3
%         chans = nPM+nM1+1:nPM+nM1+nS1;
%         color = color_S1(1,:);
%         colorConf = color_S1w;
%     end
% 
%     subplot(1,3,ctx)
%     
%     binrange = 0:10:70;
%     a = maxFRall(chans,:);
%     a(a>70) = 70;
%     
%     allCounts = [];
%     for cond1 = 1:nCond
%         [allCounts(cond1,:)] = histc(a(:,cond1),binrange);
%     end
%        
%     allCounts = allCounts/sum(sum(allCounts))*100;
%     b1 = bar(binrange,allCounts','stack','EdgeColor','none');
%     colormap(color_Cyber(1:nCond,:));   
%     hold on
%         
%     set(gca,'TickDir','out')
% end
%         
% set(gcf,'Position',[143   348   848   213],'Renderer','Painters')

%% Plot average firing rate and modulation depth across days
newOrder = [5,1,2,3,4];
for i = 1:length(newOrder)
    newNames{i,:} = condNames{newOrder(i),:};
end
figure
for  ctx = 1:3
    
    subplot(1,3,ctx)
    
    for Day = 1:nDays
        nPM = Results(Day).nChans(1);
        nM1 = Results(Day).nChans(2);
        nS1 = Results(Day).nChans(3);
        
        if ctx == 1
            chans = 1:nPM;
            color = color_PM(1,:);
            colorConf = color_PMw;
        elseif ctx == 2
            chans = nPM+1:nPM+nM1;
            color = color_M1(1,:);
            colorConf = color_M1w;
        elseif ctx == 3
            chans = nPM+nM1+1:nPM+nM1+nS1;
            color = color_S1(1,:);
            colorConf = color_S1w;
        end
        
        meanFRdays(Day,:,ctx) = mean(Results(Day).meanFR(chans,:));
    end
    toPlot = squeeze(meanFRdays(:,newOrder,ctx));
    mybar(toPlot,'SD',[],'unpaired',newNames','MeanFR',color_Cyber(newOrder,:));
    ylim([0 50])
    set(gca,'TickDir','out')
    
    % Wilcoxon paired test between all different combinaison
    if f_stat_between_each_combinaition == 1
        cortex={'PMd','M1','S1'};
        stufftocompare = nchoosek(1:5,2);
        for pairofcomp = 1:size(stufftocompare,1)
            compNames{1,pairofcomp} = [newNames{stufftocompare(pairofcomp,1)} '_' newNames{stufftocompare(pairofcomp,2)}];
            pvaluesMeanFR.(cortex{ctx}).(compNames{pairofcomp})=[];
            toPlotPaired = [];
            toPlotPaired = toPlot(:,stufftocompare(pairofcomp,:));
            
            pvaluesMeanFR.(cortex{ctx}).(compNames{pairofcomp})=mybar(toPlotPaired,'SD',[],'paired',newNames','MeanFR',color_Cyber(newOrder,:));
            ylim([0 50])
            set(gca,'TickDir','out')
            close
        end
    end
    
end

set(gcf,'Position',[6         300        862   213],'Renderer','Painters')


figure
for  ctx = 1:3
    
    subplot(1,3,ctx)
    
    for Day = 1:nDays
        nPM = Results(Day).nChans(1);
        nM1 = Results(Day).nChans(2);
        nS1 = Results(Day).nChans(3);
        
        if ctx == 1
            chans = 1:nPM;
            color = color_PM(1,:);
            colorConf = color_PMw;
        elseif ctx == 2
            chans = nPM+1:nPM+nM1;
            color = color_M1(1,:);
            colorConf = color_M1w;
        elseif ctx == 3
            chans = nPM+nM1+1:nPM+nM1+nS1;
            color = color_S1(1,:);
            colorConf = color_S1w;
        end
        
        modDepthdays(Day,:,ctx) = mean(Results(Day).modDepth(chans,:));
    end
    toPlot = squeeze(modDepthdays(:,newOrder,ctx));
    mybar(toPlot,'SD',[],'unpaired',newNames','ModDepth',color_Cyber(newOrder,:))
    ylim([0 70])
    set(gca,'TickDir','out')
    
    % Wilcoxon paired test between all different combinaison
    if f_stat_between_each_combinaition == 1
        cortex={'PMd','M1','S1'};
        stufftocompare = nchoosek(1:5,2);
        for pairofcomp = 1:size(stufftocompare,1)
            compNames{1,pairofcomp} = [newNames{stufftocompare(pairofcomp,1)} '_' newNames{stufftocompare(pairofcomp,2)}];
            pvaluesModDepth.(cortex{ctx}).(compNames{pairofcomp})=[];
            toPlotPaired = [];
            toPlotPaired = toPlot(:,stufftocompare(pairofcomp,:));
            
            pvaluesModDepth.(cortex{ctx}).(compNames{pairofcomp})=mybar(toPlotPaired,'SD',[],'paired',newNames','MeanFR',color_Cyber(newOrder,:));
            ylim([0 50])
            set(gca,'TickDir','out')
            close
        end
        
    end
  
end

set(gcf,'Position',[6         10       862   213],'Renderer','Painters')
if f_stat_between_each_combinaition == 1
    close
end
%% Comparison bewteen cortex

modDepthtoPlotCtx = squeeze(mean(modDepthdays,2));
meanFRdaystoPlotCtx = squeeze(mean(meanFRdays,2));

%Colormap of brain
cmapbrain=[61,184,224;... %[82,163,204]/255; % S1
    238,98,109;... %[46,115,69]/255; % M1    
    158,55,123]; %[168,126,229]/255; % PM
cmapbrain=cmapbrain/255;

subplot(1,2,1)
mybar(modDepthtoPlotCtx,'SD',[],'unpaired',{'PMd','M1','S1'},'modDepth',cmapbrain)
ylim([0 60])
set(gca,'TickDir','out')



subplot(1,2,2)
mybar(meanFRdaystoPlotCtx,'SD',[],'unpaired',{'PMd','M1','S1'},'meanFR',cmapbrain)
ylim([0 40])
set(gca,'TickDir','out')
set(gcf,'Position',[6         10       862   213],'Renderer','Painters')
%% stat brain
figure
mybar([modDepthtoPlotCtx(:,1) modDepthtoPlotCtx(:,2)] ,'SD',[],'paired',{'PMd','M1'},'modDepth',[])
figure
mybar([modDepthtoPlotCtx(:,1) modDepthtoPlotCtx(:,3)] ,'SD',[],'paired',{'PMd','S1'},'modDepth',[])
figure
mybar([modDepthtoPlotCtx(:,2) modDepthtoPlotCtx(:,3)] ,'SD',[],'paired',{'M1','S1'},'modDepth',[])
%%
figure
mybar([meanFRdaystoPlotCtx(:,1) meanFRdaystoPlotCtx(:,2)] ,'SD',[],'paired',{'PMd','M1'},'meanFR',[])
figure
mybar([meanFRdaystoPlotCtx(:,1) meanFRdaystoPlotCtx(:,3)] ,'SD',[],'paired',{'PMd','S1'},'meanFR',[])
figure
mybar([meanFRdaystoPlotCtx(:,2) meanFRdaystoPlotCtx(:,3)] ,'SD',[],'paired',{'M1','S1'},'meanFR',[])










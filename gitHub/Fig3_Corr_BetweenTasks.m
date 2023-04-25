% a_Corr_BetweenTasks
% This script will analyze the correlation of individual
% channels/neurons/emgs across different tasks
% INSTRUCTIONS:
% Run once for each day. The day variable will change each time
% When you're done with all days, run the plot part
clear
% close all
f_norm = 0; % Flag to normalize the correlation by the within-task correlation
f_emg = 0;
f_kin = 0;
f_singleUnit = 1;
f_plotHist = 0; % Flag to also plot histograms for distribution
muscles = [2,3,6,7,4,5,8]; % 1=GLU; 2=IL, 3=RF, 4=ST, 5=MG, 6=TA, 7=EDL, 8=FHL

load('ColorMapDarkHot.mat')
colormap_Cyber;

color_ctx = [color_PM(1,:);color_M1(1,:);color_S1(1,:)];

%%%%%%%%%%% Single Unit %%%%%%%%%%%%%%%%%
% load('sessionsElektra.mat')
% dataDir = ['/Users/ismaelseanez/Google Drive/Working Folder/G-Lab/MONKEY/'...
%     'BSI LOCOMOTION Project/Research Material/Analysis/PROCESSED DATA/Elektra/SingleUnit/'];
% f_emg = 1;
% f_singleUnit = 1;
f_hope='no'; % 'yes' or 'no'. If yes it will remove bad electrodes (outside the brain)
if regexp(f_hope,'yes')
    load('hope.mat') %load IdxElectOut variable
end

Monkey = 'Natalya';
load(['sessions' Monkey '.mat'])
%sessions(5:end,:) = [];
% dataDir = ['/Users/ismaelseanez/Google Drive/Working Folder/G-Lab/MONKEY/'...
%     'BSI LOCOMOTION Project/Research Material/Analysis/PROCESSED DATA/' Monkey '/SingleUnitClean/'];
% dataDir = ['/Users/Borgo/Google Drive/BSI LOCOMOTION Project/Research Material/Analysis/PROCESSED DATA/' Monkey '/MultiUnitClean/'];
 dataDir = ['/Users/Borgo/Documents/MATLAB/fri_mks/BSI LOCOMOTION Project_end_of_july_2021/Research Material/Analysis/PROCESSED DATA/' Monkey '/SingleUnitClean/'];

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% Multiunit with Kinematics
% sessions = ['20190415';'20190425';'20190426'];
% load('GoodKinematics.mat')
% dataDir = ['/Users/ismaelseanez/Google Drive/Working Folder/G-Lab/MONKEY/'...
%     'BSI LOCOMOTION Project/Research Material/Analysis/PROCESSED DATA/Elektra/MultiunitWithKin/'];
% dataDir = ['/Users/Borgo/Google Drive/BSI LOCOMOTION Project/Research Material/Analysis/PROCESSED DATA/Elektra/SingleUnit/'];

% f_kin = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% Multi Unit
% load('sessionsElektra.mat')
% dataDir = ['/Users/ismaelseanez/Google Drive/Working Folder/G-Lab/MONKEY/'...
%     'BSI LOCOMOTION Project/Research Material/Analysis/PROCESSED DATA/Elektra/MultiUnit/'];
% f_emg = 1;
% f_singleUnit = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nDays = size(sessions,1);
Results = struct();

for Day =  1:nDays
    clear DataAllGaits Neural_cyclic EMG_cyclic Kin_cyclic
    
        % filter parameters for kinematics
    sampFq=100; % sample frequency of the kinematic data
    half_rate = sampFq/2;
    fre_filt = 6; % Hz
    filter_order = 4;
    cutoff = fre_filt/half_rate;
    [b,a]= butter(filter_order, cutoff);
    
    
    
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
        
        % remove bad electrodes in Natalya
        if regexp(f_hope,'yes')
            cellIDtokickOut=[];
            cnt=1;
            for electout=1:size(IdxElectOut,2)
                tempidx=find(Data.SpikeEvent.M1(1,:)==IdxElectOut(electout));
                if ~isempty(tempidx)
                    cellIDtokickOut(cnt)=Data.SpikeEvent.M1(3,tempidx(1));
                    cnt=cnt+1;
                end
                clear tempidx
            end         
        end
        
        
        
        for j = 1:numel(D)
            
            Neural_cyclic(:,:,j) = D(j).Neural_cyclic;
            DataAllGaits(k).duration = length(D(j).times)/100;
            DataAllGaits(k).FOpercent = D(j).RFO/size(D(j).A,1)*100;
            DataAllGaits(k).Neural_cyclic = D(j).Neural_cyclic;
            DataAllGaits(k).condition = DataAvg(trialIndex).condition;
            
            if f_emg == 1
                EMG_cyclic(:,:,j) = D(j).EMG_cyclic;
                DataAllGaits(k).EMG_cyclic = D(j).EMG_cyclic;
            end
            
            if f_kin == 1
                Kin_cyclic(:,:,j) = filtfilt(b,a,D(j).Kin_cyclic);
                DataAllGaits(k).Kin_cyclic = filtfilt(b,a,D(j).Kin_cyclic);
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
    nCond = 5;
    [~, idx_Cond,~] = f_findGaitConditions(DataAllGaits,nCond);
    
    %% Get average behavior for each condition
    nChan = size(DataAllGaits(1).Neural_cyclic,2);
    lengthNeural = size(DataAllGaits(1).Neural_cyclic,1);
    if f_emg == 1
        nEMG = size(DataAllGaits(1).EMG_cyclic,2);
        lengthEMG = size(DataAllGaits(1).EMG_cyclic,1);
    end
    
    % Get the number of gait cycles for each condition
    for i = 1:nCond
        nGaitsCond(i) = sum(idx_Cond(i,:));
    end
    
    maxGaits = max(nGaitsCond); % Maximum number of gaits
    
    % Create firingRate and EMG arrays
    firingRatesAll = nan(nChan,nCond,lengthNeural,maxGaits);
    
    if f_emg == 1
        EMGs = nan(nEMG,nCond,lengthEMG,maxGaits);
    end
    
    % Populate arrays
    countCond = ones(1,nCond);
    for gait = 1:numel(DataAllGaits)
        cond = find(idx_Cond(:,gait)==1);
        firingRatesAll(:,cond,:,countCond(cond)) = DataAllGaits(gait).Neural_cyclic';
        if f_emg == 1
            EMGs(:,cond,:,countCond(cond)) = DataAllGaits(gait).EMG_cyclic';
        end
        countCond(cond) = countCond(cond)+1;
    end
    
    % Get the average for each condition
    firingRatesAvg = mean(firingRatesAll,4,'omitnan');
    firingRatesSTE = std(firingRatesAll,0,4,'omitnan')./sqrt(size(firingRatesAll,4));
    
    
    % Get the normalized FR
    firingRatesNorm = bsxfun(@rdivide,firingRatesAvg,max(firingRatesAvg,[],3));
    % firingRatesNorm = firingRatesAvg;
    
    if f_emg == 1
        EMGsAvg = mean(EMGs,4,'omitnan');
        EMGsSTE = std(EMGs,0,4,'omitnan')./sqrt(size(EMGs,4));
        EMGsNorm = bsxfun(@rdivide,EMGsAvg,max(EMGsAvg,[],3));
    end
    
    % Get the peak location for each condition for sorting
    if nChan == 160 % Number of channels
        nPM = 48; % Number of neurons on each cortex
        nM1 = 48;
        nS1 = 64;
    elseif nChan <160
        nPM = Data.nCellPM; % Number of neurons on each cortex
        if regexp(f_hope,'yes')
            nM1 = Data.nCellM1 - length(cellIDtokickOut);
            nChan = nChan - length(cellIDtokickOut);
        else
            nM1 = Data.nCellM1;
        end
        nS1 = Data.nCellS1;
    elseif nChan == 192
        nS1 = 32;
        nM1 = 64;
        nPE = 32;
        nPM = 64;
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
    
    % Demean the data
%     FR_All = bsxfun(@minus,firingRatesAll,mean(firingRatesAll,3,'omitnan'));
%     FR_Avg = bsxfun(@minus,firingRatesAvg,mean(firingRatesAvg,3,'omitnan'));
    FR_All = firingRatesAll;
    FR_Avg = firingRatesAvg;
    
    if f_emg == 1
        EMG_All = bsxfun(@minus,EMGs,mean(EMGs,3,'omitnan'));
        EMG_Avg = bsxfun(@minus,EMGsAvg,mean(EMGsAvg,3,'omitnan'));
    end
    
    %% Remove cells from bad electrodes in Natalya
    if regexp(f_hope,'yes')
        cellIDtokickOut2=cellIDtokickOut+nPM;
        FR_Avg(cellIDtokickOut2,:,:)=[];     
    end
    
    %% Calculate corr coeff for each neuron/channel across tasks
    corr_All = zeros(nChan,nCond,nCond);
    for cond1 = 1:nCond
        
        % To split data in two sets randomly
        nTrials = countCond(cond1)-1;
        temp = 1:nTrials;
        randIdx = temp(randperm(length(temp))); % randomize order
        cutoff = floor(nTrials/2); % cut halfway
        
        for cond2 = 1:nCond
            for chan = 1:nChan
                if cond1 ~= cond2
                    corr_All(chan,cond1,cond2) = corr(squeeze(FR_Avg(chan,cond1,:)),...
                        squeeze(FR_Avg(chan,cond2,:)));
                    
                    % When it's the same condition, look at the average intra-task
                    % correlation of each channel
                elseif cond1 == cond2
                    
                    % Get the behavior for all those random steps
                    a = FR_All(chan,cond1,:,randIdx(1:cutoff));
                    b = FR_All(chan,cond1,:,randIdx(cutoff+1:end));
                    
                    % And the average
                    amean = squeeze(mean(a,4,'omitnan'));
                    bmean = squeeze(mean(b,4,'omitnan'));
                    
                    corr_All(chan,cond1,cond2)  = corr(amean,bmean);
                end
                
                
            end
        end
    end
    
    %% Calculate corr coeff for each EMG across tasks
    if f_emg == 1
        corr_AllEmg = zeros(nEMG,nCond,nCond);
        for cond1 = 1:nCond
            
            % To split data in two sets randomly
            nTrials = countCond(cond1)-1;
            temp = 1:nTrials;
            randIdx = temp(randperm(length(temp))); % randomize order
            cutoff = floor(nTrials/2); % cut halfway
            
            for cond2 = 1:nCond
                for chan = 1:nEMG
                    if cond1 ~= cond2
                        corr_AllEmg(chan,cond1,cond2) = corr(squeeze(EMG_Avg(chan,cond1,:)),...
                            squeeze(EMG_Avg(chan,cond2,:)));
                        
                        % When it's the same condition, look at the average intra-task
                        % correlation of each channel
                    elseif cond1 == cond2
                        
                        % Get the behavior for all those random steps
                        a = EMG_All(chan,cond1,:,randIdx(1:cutoff));
                        b = EMG_All(chan,cond1,:,randIdx(cutoff+1:end));
                        
                        % And the average
                        amean = squeeze(mean(a,4,'omitnan'));
                        bmean = squeeze(mean(b,4,'omitnan'));
                        
                        corr_AllEmg(chan,cond1,cond2)  = corr(amean,bmean);
                    end
                end
            end
        end
    end
    %% Calculate corr coeff for each Kin across tasks
    if f_kin == 1
        corr_AllKin = zeros(nKin,nCond,nCond);
        for cond1 = 1:nCond
            
            % Get new number of trials
            
            % To split data in two sets randomly
            nTrials = countCond(cond1)-1;
            temp = 1:nTrials;
            randIdx = temp(randperm(length(temp))); % randomize order
            cutoff = floor(nTrials/2); % cut halfway
            
            for cond2 = 1:nCond
                for chan = 1:nKin
                    if cond1 ~= cond2
                        corr_AllKin(chan,cond1,cond2) = corr(squeeze(Kin_Avg(chan,cond1,:)),...
                            squeeze(Kin_Avg(chan,cond2,:)));
                        
                        % When it's the same condition, look at the average intra-task
                        % correlation of each channel
                    elseif cond1 == cond2
                        
                        % Get the behavior for all those random steps
                        a = KinAll(chan,cond1,:,randIdx(1:cutoff));
                        b = KinAll(chan,cond1,:,randIdx(cutoff+1:end));
                        
                        % And the average
                        amean = squeeze(mean(a,4,'omitnan'));
                        bmean = squeeze(mean(b,4,'omitnan'));
                        
                        corr_AllKin(chan,cond1,cond2)  = corr(amean,bmean);
                    end
                end
            end
        end
    end
    %%%%%%%%%%%% POPULATION ANAYSIS %%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Calculate corr coeff for each feature in FR population
    corr_popuFR = zeros(nChan,nCond,nCond);
    for ctx = 1:3 % 1 for Premotor, 2 for Motor, 3 for Sensory
        
        if ctx == 1
            chansCtx = 1:nPM;
        elseif ctx == 2
            chansCtx = nPM+1:nPM+nM1;
        elseif ctx == 3
            chansCtx = nPM+nM1+1:nPM+nM1+nS1;
        end
        neurons = chansCtx;
        corr_popuCtx = zeros(length(neurons),nCond,nCond);
        [coeff_FR,pcaScore_FR] = m_PCA(FR_Avg(neurons,:,:),1:length(neurons));
        
        for cond1 = 1:nCond
            
            % To split data in two sets randomly
            nTrials = countCond(cond1)-1;
            temp = 1:nTrials;
            randIdx = temp(randperm(length(temp))); % randomize order
            cutoff = floor(nTrials/2); % cut halfway
            
            a = squeeze(FR_All(neurons,cond1,:,randIdx(1:cutoff)));
            b = squeeze(FR_All(neurons,cond1,:,randIdx(cutoff+1:end)));
            
            aAvg = mean(a,3);
            bAvg = mean(b,3);
            
            % reconstructed population response
            Za = aAvg' * coeff_FR;
            Zb = bAvg' * coeff_FR;
            
            for cond2 = 1:nCond
                for chan = 1:size(pcaScore_FR,2)
                    if cond1 ~= cond2
                        
                        corr_popuCtx(chan,cond1,cond2) = corr(squeeze(pcaScore_FR(:,chan,cond1)),...
                            squeeze(pcaScore_FR(:,chan,cond2)));
                        
                        % When it's the same condition, look at the average intra-task
                        % correlation of each channel
                    elseif cond1 == cond2
                        
                        % Get the behavior for all those random steps
                        a = Za(:,chan);
                        b = Zb(:,chan);
                        
                        corr_popuCtx(chan,cond1,cond2)  = corr(a,b);
                    end
                end
            end
        end
        corr_popuFR(neurons,:,:) = corr_popuCtx;
    end
    
    
    %% Calculate corr coeff for each feature in EMG population
    if f_emg == 1
        corr_popuEmg = zeros(length(muscles),nCond,nCond);
        [coeff_emg,pcaScore_EMG] = m_PCA(EMG_Avg,muscles);
        
        for cond1 = 1:nCond
            
            % To split data in two sets randomly
            nTrials = countCond(cond1)-1;
            temp = 1:nTrials;
            randIdx = temp(randperm(length(temp))); % randomize order
            cutoff = floor(nTrials/2); % cut halfway
            
            a = squeeze(EMG_All(muscles,cond1,:,randIdx(1:cutoff)));
            b = squeeze(EMG_All(muscles,cond1,:,randIdx(cutoff+1:end)));
            
            aAvg = mean(a,3);
            bAvg = mean(b,3);
            
            % reconstructed population response
            Za = aAvg' * coeff_emg;
            Zb = bAvg' * coeff_emg;
            
            for cond2 = 1:nCond
                for chan = 1:size(pcaScore_EMG,2)
                    if cond1 ~= cond2
                        
                        corr_popuEmg(chan,cond1,cond2) = corr(squeeze(pcaScore_EMG(:,chan,cond1)),...
                            squeeze(pcaScore_EMG(:,chan,cond2)));
                        
                        % When it's the same condition, look at the average intra-task
                        % correlation of each channel
                    elseif cond1 == cond2
                        
                        % Get the behavior for all those random steps
                        a = Za(:,chan);
                        b = Zb(:,chan);
                        
                        corr_popuEmg(chan,cond1,cond2)  = corr(a,b);
                    end
                end
            end
        end
    end
    
    if f_singleUnit == 1
        % Save results
        Results(Day).corr_All = corr_All;
        Results(Day).nChans = [nPM,nM1,nS1];
        
        if f_emg == 1
            Results(Day).corr_AllEmg = corr_AllEmg;
        end
        if f_kin ==1
            Results(Day).corr_AllKin = corr_AllKin;
        end
        
    else
        Results(Day).corr_All = corr_All;
        Results(Day).nChans = [nPM,nM1,nS1];
        if ~exist('corr_AllDays','var')
            corr_AllDays = zeros(nChan,nCond,nCond,nDays);
            corr_AllDaysPopuFR = zeros(size(FR_Avg,1),nCond,nCond,nDays);
            corr_AllDays(:,:,:,Day) = corr_All;
            corr_AllDaysPopuFR(:,:,:,Day) = corr_popuFR;
            
            if f_emg == 1
                corr_AllDaysEmg = zeros(nEMG,nCond,nCond,nDays);
                corr_AllDaysPopuEMG = zeros(length(muscles),nCond,nCond,nDays);
                corr_AllDaysEmg(:,:,:,Day) = corr_AllEmg;
                corr_AllDaysPopuEMG(:,:,:,Day) = corr_popuEmg;
            end
            
            if f_kin == 1
                corr_AllDaysKin = zeros(nKin,nCond,nCond,nDays);
                corr_AllDaysKin(:,:,:,Day) = corr_AllKin;
            end
            
        else
            corr_AllDays(:,:,:,Day) = corr_All;
            corr_AllDaysPopuFR(:,:,:,Day) = corr_popuFR;
            
            if f_emg == 1
                corr_AllDaysEmg(:,:,:,Day) = corr_AllEmg;
                corr_AllDaysPopuEMG(:,:,:,Day) = corr_popuEmg;
            end
            
            if f_kin == 1
                corr_AllDaysKin(:,:,:,Day) = corr_AllKin;
            end
        end
    end
end % End of days                 (heyoooo!)

%% Plot histogram for each cortex but keeping conditions separate
if f_singleUnit == 1
    figure
    combos = nchoosek(1:nCond,2); % unique combinations
    nCombos = length(combos);
    if f_norm == 1
        binrange = -2:0.4:2;
    else
        binrange = -1:0.2:1;
    end
    for ctx = 1:3 % 1 for Premotor, 2 for Motor, 3 for Sensory
        
        clear corr_AlDays
        
        % Get correlations for all days of cortex channels into one matrix
        corr_AllDays = [];
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
            
            corrAvgChans(:,:,Day) = squeeze(mean(Results(Day).corr_All(chans,:,:)));
            
        end
        
        allCounts = zeros(nCombos,length(binrange));
        corrHist = [];
        for j = 1:nCombos
            cond1 = combos(j,1);
            cond2 = combos(j,2);
            
            if f_norm == 1
                a = squeeze(corrAvgChans(cond1,cond2,:)); % comparison correlation
                
                b = squeeze(corrAvgChans(cond1,cond1,:)); % correlation within task 1
                c = squeeze(corrAvgChans(cond2,cond2,:)); % correlation within task 2
                d = mean([b;c],'omitnan'); % Get the average correlation of both tasks being compared
                a = a./d; % divide over average
            else
                a = squeeze(corrAvgChans(cond1,cond2,:));
                
                aNorm = squeeze(corrAvgChans(cond1,cond2,:)); % comparison correlation
                
                b = squeeze(corrAvgChans(cond1,cond1,:)); % correlation within task 1
                c = squeeze(corrAvgChans(cond2,cond2,:)); % correlation within task 2
                d = mean([b;c],'omitnan'); % Get the average correlation of both tasks being compared
                aNorm = aNorm./d; % divide over average
            end
            allAs(j,:) = a;
            allAsNorm(j,:) = aNorm;
            allCounts(j,:) = histc(a,binrange);
        end
        
        meanCounts = mean(mean(allAs,'omitnan'));
        
        
        allCounts = allCounts/sum(sum(allCounts))*100;
        
        sub(ctx) = subplot(1,3,ctx);
        h = bar(binrange,allCounts','stack','EdgeColor','none');
        colormap jet
        
        % Add mean
        hold on
        plot(meanCounts,58,'o','MarkerFaceColor',color,'MarkerEdgeColor','none')
        
        % Plot a line at the normalized correlation
        if f_norm == 0
            meanNorm = mean(mean(allAsNorm,'omitnan'));
            plot([meanNorm meanNorm],[0 60],'--','color',[160 160 160]/255,'LineWidth',0.5);
        end
        
        
        h = gca;
        h.TickDir = 'out';
        if f_norm == 1
            xlim([-1 2])
            xticks(-2:0.5:2)
            xlabel('Normalized correlation across tasks')
        else
            xlim([-0.5 1])
            xticks(-1:0.25:1)
            xlabel('Correlation across tasks')
        end
        
        if ctx == 1
            ylabel('Day comparisons (%)')
        end
        
        
    end
else
    
    figure
    combos = nchoosek(1:nCond,2); % unique combinations
    nCombos = length(combos);
    if f_norm == 1
        binrange = -2:0.4:2;
    else
        binrange = -1:0.2:1;
    end
    for ctx = 1:3 % 1 for Premotor, 2 for Motor, 3 for Sensory
        
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
        
        % Get correlations for all days of cortex channels into one matrix
        corrAvgChans = squeeze(mean(corr_AllDays(chans,:,:,:),1)); %average across channels for each day
        
        allCounts = zeros(nCombos,length(binrange));
        corrHist = [];
        for j = 1:nCombos
            cond1 = combos(j,1);
            cond2 = combos(j,2);
            
            if f_norm == 1
                a = squeeze(corrAvgChans(cond1,cond2,:)); % comparison correlation
                
                b = squeeze(corrAvgChans(cond1,cond1,:)); % correlation within task 1
                c = squeeze(corrAvgChans(cond2,cond2,:)); % correlation within task 2
                d = mean([b;c],'omitnan'); % Get the average correlation of both tasks being compared
                a = a./d; % divide over average
            else
                a = squeeze(corrAvgChans(cond1,cond2,:));
                
                aNorm = squeeze(corrAvgChans(cond1,cond2,:)); % comparison correlation
                
                b = squeeze(corrAvgChans(cond1,cond1,:)); % correlation within task 1
                c = squeeze(corrAvgChans(cond2,cond2,:)); % correlation within task 2
                d = mean([b;c],'omitnan'); % Get the average correlation of both tasks being compared
                aNorm = aNorm./d; % divide over average
            end
            allAs(j,:) = a;
            allAsNorm(j,:) = aNorm;
            allCounts(j,:) = histc(a,binrange);
        end
        
        meanCounts = mean(mean(allAs,'omitnan'));
        
        
        allCounts = allCounts/sum(sum(allCounts))*100;
        
        sub(ctx) = subplot(1,3,ctx);
        h = bar(binrange,allCounts','stack','EdgeColor','none');
        colormap jet
        
        % Add mean
        hold on
        plot(meanCounts,58,'o','MarkerFaceColor',color,'MarkerEdgeColor','none')
        
        % Plot a line at the normalized correlation
        if f_norm == 0
            meanNorm = mean(mean(allAsNorm,'omitnan'));
            plot([meanNorm meanNorm],[0 60],'--','color',[160 160 160]/255,'LineWidth',0.5);
        end
        
        
        h = gca;
        h.TickDir = 'out';
        if f_norm == 1
            xlim([-1 2])
            xticks(-2:0.5:2)
            xlabel('Normalized correlation across tasks')
        else
            xlim([-0.5 1])
            xticks(-1:0.25:1)
            xlabel('Correlation across tasks')
        end
        
        if ctx == 1
            ylabel('Day comparisons (%)')
        end
        
        
    end
end

linkaxes(sub,'y')
set(gcf,'Position',[219   213   851   209],'Renderer','painters')

%% Plot histogram for correlations for each cortex
figure
combos = nchoosek(1:nCond,2); % unique combinations
nCombos = length(combos);
if f_norm == 1
    binrange = -2:0.4:2;
else
    binrange = -1:0.2:1;
end
for ctx = 1:3 % 1 for Premotor, 2 for Motor, 3 for Sensory
    
    clear corr_AlDays
    
    % Get correlations for all days of cortex channels into one matrix
    corr_AllDays = [];
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
        
        corrAvgChans(:,:,Day) = squeeze(mean(Results(Day).corr_All(chans,:,:)));
        
    end
    
    allCounts = zeros(nCombos,length(binrange));
    corrHist = [];
    for j = 1:nCombos
        cond1 = combos(j,1);
        cond2 = combos(j,2);
        
        if f_norm == 1
            a = squeeze(corrAvgChans(cond1,cond2,:)); % comparison correlation
            
            b = squeeze(corrAvgChans(cond1,cond1,:)); % correlation within task 1
            c = squeeze(corrAvgChans(cond2,cond2,:)); % correlation within task 2
            d = mean([b;c],'omitnan'); % Get the average correlation of both tasks being compared
            a = a./d; % divide over average
        else
            a = squeeze(corrAvgChans(cond1,cond2,:));
            
            aNorm = squeeze(corrAvgChans(cond1,cond2,:)); % comparison correlation
            
            b = squeeze(corrAvgChans(cond1,cond1,:)); % correlation within task 1
            c = squeeze(corrAvgChans(cond2,cond2,:)); % correlation within task 2
            d = mean([b;c],'omitnan'); % Get the average correlation of both tasks being compared
            aNorm = aNorm./d; % divide over average
        end
        allAs(j,:) = a;
        allAsNorm(j,:) = aNorm;
        allCounts(j,:) = histc(a,binrange);
    end
    
    sub(ctx) = subplot(1,3,ctx);
    histogram(reshape(allAs,1,[]),binrange,'Normalization','probability',...
        'FaceColor',color,'EdgeColor','none')
    
    meanCounts = mean(reshape(allAs,1,[]));
    stdCounts = std(reshape(allAs,1,[]));
    
    % Add mean and std
    hold on
    plot(meanCounts,.55,'o','MarkerFaceColor',color,'MarkerEdgeColor','none')
    plot([meanCounts-stdCounts, meanCounts+stdCounts],[.55 .55],'color',color,'LineWidth',0.5)
    
    % Plot a line at the normalized correlation
    if f_norm == 0
        meanNorm = mean(mean(allAsNorm,'omitnan'));
        plot([meanNorm meanNorm],[0 .6],'--','color',[160 160 160]/255,'LineWidth',0.5);
    end
    
    
    h = gca;
    h.TickDir = 'out';
    if f_norm == 1
        xlim([-1 2])
        xticks(-2:0.5:2)
        xlabel('Normalized correlation across tasks')
    else
        xlim([-0.5 1])
        xticks(-1:0.25:1)
        xlabel('Correlation across tasks')
    end
    
    if ctx == 1
        ylabel('Day comparisons (%)')
    end
    
    
end
linkaxes(sub,'y')
set(gcf,'Position',[219   213   851   209],'Renderer','painters')

%% Plot BARPLOT for correlations for each cortex

combos = nchoosek(1:nCond,2); % unique combinations
nCombos = length(combos);
if f_norm == 1
    binrange = -2:0.4:2;
else
    binrange = -1:0.2:1;
end
for ctx = 1:3 % 1 for Premotor, 2 for Motor, 3 for Sensory
    
    clear corr_AlDays
    
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
        
        corrAvgChans(:,:,Day) = squeeze(mean(Results(Day).corr_All(chans,:,:)));
    end
    
    % Get the average correlation for each day
    for j = 1:nCombos
        cond1 = combos(j,1);
        cond2 = combos(j,2);
        
        if f_norm == 1
            a = squeeze(corrAvgChans(cond1,cond2,:)); % comparison correlation
            
            b = squeeze(corrAvgChans(cond1,cond1,:)); % correlation within task 1
            c = squeeze(corrAvgChans(cond2,cond2,:)); % correlation within task 2
            d = mean([b;c],'omitnan'); % Get the average correlation of both tasks being compared
            a = a./d; % divide over average
        else
            a = squeeze(corrAvgChans(cond1,cond2,:));
            
            aNorm = squeeze(corrAvgChans(cond1,cond2,:)); % comparison correlation
            
            b = squeeze(corrAvgChans(cond1,cond1,:)); % correlation within task 1
            c = squeeze(corrAvgChans(cond2,cond2,:)); % correlation within task 2
            d = mean([b;c],'omitnan'); % Get the average correlation of both tasks being compared
            aNorm = aNorm./d; % divide over average
        end
        allAs(j,:) = a;
    end
    
    corr_AllDays(:,ctx) = mean(allAs,1)';
    
end

% Mybar that shit
figure
mybar(corr_AllDays,'SD',[0 1],'unpaired',{'PM','M1','S1'},'Correlation (R)',color_ctx);
h = gca;
h.TickDir = 'out';
set(gcf,'Position',[360   278   560   420])

[~,pVal_PMM1] = quickStats(corr_AllDays(:,1),corr_AllDays(:,2),'paired')
[~,pVal_PMS1] = quickStats(corr_AllDays(:,1),corr_AllDays(:,3),'paired')
[~,pVal_M1S1] = quickStats(corr_AllDays(:,2),corr_AllDays(:,3),'paired')

avg = mean(corr_AllDays,1)
sd = std(corr_AllDays,1)
%% Plot histogram for all muscles but keeping conditions separate
figure
combos = nchoosek(1:nCond,2); % unique combinations
nCombos = length(combos);
if f_norm == 1
    binrange = -2:0.4:2;
else
    binrange = -1:0.2:1;
end
for ctx = 1 % 1 for Premotor, 2 for Motor, 3 for Sensory
    
    chans = muscles;
    color = color_EMG(1,:);
    
    % Get correlations for all days of cortex channels into one matrix
    corrAvgChans = squeeze(mean(corr_AllDaysEmg(chans,:,:,:),1)); %average across channels for each day

    allCounts = zeros(nCombos,length(binrange));
    corrHist = [];
    for j = 1:nCombos
        cond1 = combos(j,1);
        cond2 = combos(j,2);
        
        if f_norm == 1
            a = squeeze(corrAvgChans(cond1,cond2,:)); % comparison correlation
            
            b = squeeze(corrAvgChans(cond1,cond1,:)); % correlation within task 1
            c = squeeze(corrAvgChans(cond2,cond2,:)); % correlation within task 2
            d = mean([b;c],'omitnan'); % Get the average correlation of both tasks being compared
            a = a./d; % divide over average
        else
            a = squeeze(corrAvgChans(cond1,cond2,:));
            
            aNorm = squeeze(corrAvgChans(cond1,cond2,:)); % comparison correlation
            
            b = squeeze(corrAvgChans(cond1,cond1,:)); % correlation within task 1
            c = squeeze(corrAvgChans(cond2,cond2,:)); % correlation within task 2
            d = mean([b;c],'omitnan'); % Get the average correlation of both tasks being compared
            aNorm = aNorm./d; % divide over average
        end
        allAs(j,:) = a;
        allAsNorm(j,:) = aNorm;
        allCounts(j,:) = histc(a,binrange);
    end
    
    meanCounts = mean(mean(allAs,'omitnan'));
    
    
    allCounts = allCounts/sum(sum(allCounts))*100;
    
    sub(ctx) = subplot(1,3,ctx);
    h = bar(binrange,allCounts','stack','EdgeColor','none');
    colormap jet
    
    % Add mean
    hold on
    plot(meanCounts,58,'o','MarkerFaceColor',color,'MarkerEdgeColor','none')
    
    % Plot a line at the normalized correlation
    if f_norm == 0
        meanNorm = mean(mean(allAsNorm,'omitnan'));
        plot([meanNorm meanNorm],[0 60],'--','color',[160 160 160]/255,'LineWidth',0.5);
    end
    
    
    h = gca;
    h.TickDir = 'out';
    if f_norm == 1
        xlim([-1 2])
        xticks(-2:0.5:2)
        xlabel('Normalized correlation across tasks')
    else
        xlim([-0.5 1])
        xticks(-1:0.25:1)
        xlabel('Correlation in EMG across tasks')
    end
    
    if ctx == 1
        ylabel('Day comparisons (%)')
    end
    
    
end

linkaxes(sub,'y')
set(gcf,'Position',[219   213   851   209],'Renderer','painters')

%% Plot EMG
EMGmap_Elektra
f_plotEMG = 1;
if f_plotEMG == 1
    load('Natalya_KinNames.mat')
    figure
    for i = 1:size(EMG_Avg,1)
        subplot(4,5,i)
        for cond = 1:nCond
            a = squeeze(EMG_Avg(i,cond,:));
            plot(a,'color',color_Cyber(cond,:))
            hold on
        end
        title(trialInfo.EMGnames{i,1})
    end
    set(gcf,'Position', [79          35        1133         646],'Renderer','Painters')
end
%% Plot histogram for all kinematics but keeping conditions separate
figure
combos = nchoosek(1:nCond,2); % unique combinations
nCombos = length(combos);
if f_norm == 1
    binrange = -2:0.4:2;
else
    binrange = -1:0.2:1;
end
for ctx = 1 % 1 for Premotor, 2 for Motor, 3 for Sensory
    
    chans = goodKin;
    color = color_EMG(2,:);
    
    % Get correlations for all days of cortex channels into one matrix
    corrAvgChans = squeeze(mean(corr_AllDaysKin(chans,:,:,:),1)); %average across channels for each day
    
    allCounts = zeros(nCombos,length(binrange));
    corrHist = [];
    for j = 1:nCombos
        cond1 = combos(j,1);
        cond2 = combos(j,2);
        
        if f_norm == 1
            a = squeeze(corrAvgChans(cond1,cond2,:)); % comparison correlation
            
            b = squeeze(corrAvgChans(cond1,cond1,:)); % correlation within task 1
            c = squeeze(corrAvgChans(cond2,cond2,:)); % correlation within task 2
            d = mean([b;c],'omitnan'); % Get the average correlation of both tasks being compared
            a = a./d; % divide over average
        else
            a = squeeze(corrAvgChans(cond1,cond2,:));
            
            aNorm = squeeze(corrAvgChans(cond1,cond2,:)); % comparison correlation
            
            b = squeeze(corrAvgChans(cond1,cond1,:)); % correlation within task 1
            c = squeeze(corrAvgChans(cond2,cond2,:)); % correlation within task 2
            d = mean([b;c],'omitnan'); % Get the average correlation of both tasks being compared
            aNorm = aNorm./d; % divide over average
        end
        allAs(j,:) = a;
        allAsNorm(j,:) = aNorm;
        allCounts(j,:) = histc(a,binrange);
    end
    
    meanCounts = mean(mean(allAs,'omitnan'));
    
    
    allCounts = allCounts/sum(sum(allCounts))*100;
    
    sub(ctx) = subplot(1,3,ctx);
    h = bar(binrange,allCounts','stack','EdgeColor','none');
    colormap jet
    
    % Add mean
    hold on
    plot(meanCounts,58,'o','MarkerFaceColor',color,'MarkerEdgeColor','none')
    
    % Plot a line at the normalized correlation
    if f_norm == 0
        meanNorm = mean(mean(allAsNorm,'omitnan'));
        plot([meanNorm meanNorm],[0 60],'--','color',[160 160 160]/255,'LineWidth',0.5);
    end
    
    
    h = gca;
    h.TickDir = 'out';
    if f_norm == 1
        xlim([-1 2])
        xticks(-2:0.5:2)
        xlabel('Normalized correlation across tasks')
    else
        xlim([-0.5 1])
        xticks(-1:0.25:1)
        xlabel('Correlation in EMG across tasks')
    end
    
    if ctx == 1
        ylabel('Day comparisons (%)')
    end
    
    
end
compNames        = cell(1,size(combos,1));
for p = 1:size(combos,1)
    compNames{p} = [condNames{combos(p,1),:} ' vs ' ...
        condNames{combos(p,2),:} ];
end
legend(compNames,'Location','NorthWest','FontSize',8), legend boxoff
linkaxes(sub,'y')
set(gcf,'Position',[219   213   851   209],'Renderer','painters')

%% Plot Kinematics
f_plotKin = 1;
if f_plotKin == 1
    load('Natalya_KinNames.mat')
    figure
    for i = 1:size(Kin_Avg,1)
        subplot(4,5,i)
        for cond = 1:nCond
            a = squeeze(Kin_Avg(i,cond,:));
            plot(a,'color',color_Cyber(cond,:))
            hold on
        end
        title(kinNames{i,1})
    end
    set(gcf,'Position', [79          35        1133         646],'Renderer','Painters')
end

%% Plot confussion matrix for each cortex and comparison

% for ctx = 1:3 % 1 for Premotor, 2 for Motor, 3 for Sensory
%
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
%     % Get correlations for all days of cortex channels into one matrix
%     correlations = [];
%     for i = 1:nDays
%         correlations = [correlations; squeeze(corr_AllDays(chans,:,:,i))];
%     end
%
%     meanCorr = squeeze(mean(correlations,1,'omitnan'));
%
%     % normalize by the variability of that condition
%     if f_norm == 1
%         for col = 1:nCond
%             meanCorr(:,col) = meanCorr(:,col)/meanCorr(col,col);
%         end
%     end
%
%     if f_plotHist
%         % Plot histogram
%         figure
%         k = 1;
%         for cond1 = 1:nCond
%             for cond2 = 1:nCond
%                 corrThis = correlations(:,cond1,cond2);
%                 nb = calcnbins(corrThis, 'fd');
%                 sub(k) = subplot(nCond,nCond,k);
%                 h = histogram(corrThis,nb);
%                 h.EdgeColor = 'none';
%                 h.FaceColor = color;
%
%                 % Add mean
%                 hold on
%                 plot([meanCorr(cond1,cond2), meanCorr(cond1,cond2)],[0 max(h.Values)],'LineWidth',2,'color',[0 0 0])
%                 text(meanCorr(cond1,cond2),max(h.Values)+max(h.Values)*0.12,...
%                     mat2str(round(meanCorr(cond1,cond2),2)),'HorizontalAlignment','center')
%
%                 % Labels, Limits, and Legends
%                 xlim([-1 1])
%                 k = k+1;
%                 if cond2 == 1
%                     ylabel(condNames{cond1,1})
%                 end
%                 if cond1 == nCond
%                     xlabel(condNames{cond2,1})
%                 end
%             end
%         end
%         set(sub,'box','off');
%         set(gcf,'Position',[ 1  1 925 704],'Renderer','Painters')
%         clear sub corrThis correlations nb h
%     end
%
%     % Plot matrix for average
%     figure
%     imagesc(tril(meanCorr,-1))
%     colormap(colorConf)
%     colorbar
%     %     caxis([min(min(meanCorr)) 1])
%     caxis([0.2 1])
%     % Create strings
%     textStrings = num2str(round(meanCorr(:),1));
%     textStrings = strtrim(cellstr(textStrings));
%
%     % Create x and y coordinates for the strings and plot them
%     [x,y] = meshgrid(1:nCond);
%     hStrings = text(x(:),y(:),textStrings(:), ...
%         'HorizontalAlignment','center');
%
%     % Choose white or black for the text color of the strings so
%     % they can be easily seen over the background color
%     textColors = zeros(length(hStrings),4);
%     set(hStrings,{'Color'},num2cell(textColors,2));
%
%     % Setting the axis labels
%     set(gca,'XTick',1:nCond,...
%         'XTickLabel',condNames,...
%         'YTick',1:nCond,...
%         'YTickLabel',condNames,...
%         'TickLength',[0 0]);
%
%     set(gcf,'Position',[927 435 354 269],'Renderer','Painters')
%
%     %     % Plot heriarchal tree
%     %     figure
%     %     for cond1 = 1:nCond
%     %         for cond2 = 1:nCond
%     %             if cond1 ~= cond2
%     %                 plot([cond1 cond2],[meanCorr(cond1,cond2) meanCorr(cond1,cond2)],'color',color)
%     %                 hold on
%     %             end
%     %         end
%     %     end
%     %     set(gca,'YDir','reverse')
%     %     ylim([0 1])
%     %     xlim([0 8])
%     %     set(gca,'XTick',1:nCond,...
%     %     'XTickLabel',condNames)
%     %     set(gcf,'Position',[931     1   349   420],'Renderer','Painters')
%
%
% end

%% Plot summary histogram for each cortex
% figure
% combos = nchoosek(1:7,2); % unique combinations
% nCombos = length(combos);
% binrange = -1:0.2:1;
% for ctx = 1:3 % 1 for Premotor, 2 for Motor, 3 for Sensory
%
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
%     % Get correlations for all days of cortex channels into one matrix
%     corrCtx = zeros(length(chans),nCombos,nDays);
%     for j = 1:nCombos
%         cond1 = combos(j,1);
%         cond2 = combos(j,2);
%         if f_norm == 1
%             a = squeeze(corr_AllDays(chans,cond1,cond1,:)); % correlation within task 1
%             b = squeeze(corr_AllDays(chans,cond2,cond2,:)); % correlation within task 2
%             c = [a,b];
%             avgCorr = mean(mean(c));
%             corrCtx(:,j,:) = corr_AllDays(chans,cond1,cond2,:)/avgCorr;
%         else
%             corrCtx(:,j,:) = corr_AllDays(chans,cond1,cond2,:);
%         end
%     end
%
%     a = reshape(corrCtx,1,[]);
%     nb = calcnbins(a, 'fd');
%     sub(ctx) = subplot(1,3,ctx);
%     h = histogram(a,nb,'Normalization','probability');
%     h.EdgeColor = 'none';
%     h.FaceColor = color;
%     ax = gca;
%     ax.TickDir = 'out';
%     if f_norm == 1
%         xlim([-3 3])
%         xticks(-3:1:3)
%     else
%         xlim([-1 1])
%         xticks(-1:0.5:1)
%     end
%
%     % Add mean
%     hold on
%     plot(mean(a,'omitnan'),max(h.Values)+0.01,'o','MarkerFaceColor',color,'MarkerEdgeColor','none')
%
%     if ctx == 1
%         ylabel('Neural units (%)')
%     end
%     xlabel('Correlation across tasks')
%
% end
%
% linkaxes(sub,'y')
% set(gcf,'Position',[219   213   851   209],'Renderer','painters')

%% Plot histogram for the EMGs


for f_emg = 1:2 % 1 for flexors, 2 for extensors
    
    if f_emg == 1
        chans = muscles(1:4);
        color = color_EMG(2,:);
        colorConf = color_Flexw;
    elseif f_emg == 2
        chans = muscles(5:7);
        color = color_EMG(4,:);
        colorConf = color_Extw;
    end
    
    % Get correlations for all days of cortex channels into one matrix
    correlations = [];
    for i = 1:nDays
        % this is for flex (IL and EDL) and ext (MG and FHL) separatly 
        correlations = [correlations; squeeze(corr_AllDaysEmg(chans,:,:,i))];
        % this is the mean of all the correlation 
        %correlations = [correlations; squeeze(corr_AllDaysEmg(:,:,:,i))];
    end
    
    meanCorr = squeeze(mean(correlations,1,'omitnan'));
    
    % normalize by the variability of that condition
    if f_norm == 1
        for col = 1:nCond
            meanCorr(:,col) = meanCorr(:,col)/meanCorr(col,col);
        end
    end
    
    if f_plotHist
        % Plot histogram
        figure
        k = 1;
        for cond1 = 1:nCond
            for cond2 = 1:nCond
                corrThis = correlations(:,cond1,cond2);
                nb = calcnbins(corrThis, 'middle');
                sub(k) = subplot(nCond,nCond,k);
                h = histogram(corrThis,nb);
                h.EdgeColor = 'none';
                h.FaceColor = color;
                
                % Add mean
                hold on
                plot([meanCorr(cond1,cond2), meanCorr(cond1,cond2)],[0 max(h.Values)],'LineWidth',2,'color',[0 0 0])
                text(meanCorr(cond1,cond2),max(h.Values)+max(h.Values)*0.12,...
                    mat2str(round(meanCorr(cond1,cond2),2)),'HorizontalAlignment','center')
                
                % Labels, Limits, and Legends
                xlim([-1 1])
                k = k+1;
                if cond2 == 1
                    ylabel(condNames{cond1,1})
                end
                if cond1 == nCond
                    xlabel(condNames{cond2,1})
                end
            end
        end
        set(sub,'box','off');
        set(gcf,'Position',[ 1  1 925 704],'Renderer','Painters')
        clear sub corrThis correlations nb h
    end
    
    % Plot matrix for average
    figure
    %imagesc(tril(meanCorr,-1))
        imagesc(meanCorr)

    colormap(colorConf)
    colorbar
    %     caxis([min(min(meanCorr)) 1])
    caxis([0.2 1])
    % Create strings
    textStrings = num2str(round(meanCorr(:),1));
    textStrings = strtrim(cellstr(textStrings));
    
    % Create x and y coordinates for the strings and plot them
    [x,y] = meshgrid(1:nCond);
    hStrings = text(x(:),y(:),textStrings(:), ...
        'HorizontalAlignment','center');
    
    % Choose white or black for the text color of the strings so
    % they can be easily seen over the background color
    textColors = zeros(length(hStrings),4);
    set(hStrings,{'Color'},num2cell(textColors,2));
    
    % Setting the axis labels
    set(gca,'XTick',1:nCond,...
        'XTickLabel',condNames,...
        'YTick',1:nCond,...
        'YTickLabel',condNames,...
        'TickLength',[0 0]);
    
    set(gcf,'Position',[927 435 354 269],'Renderer','Painters')
    
    
end

%% Plot histogram for the KINematics

goodKin2=[1,2,3,10,11,12,13,14,18];
for f_emg = 1:2 % 1 for angles, 2 for speeds
    
    if f_emg == 1
        %chans = goodKin2(1:4);
        chans = goodKin(1:3);
        color = color_EMG(2,:);
        mycolor = [23,172,232]./255; % shitty blue for angles
        %mycolor = [150,0,0]./255; % bordeau for angles
        colorConf = zeros(48,3);
        for xx=1:3
            colorConf(:,xx)=linspace(mycolor(1,xx),1 ,48);
        end
        colorConf=flip(colorConf);
        %colorConf = color_Flexw;
    elseif f_emg == 2
        %chans = goodKin2(5:end);
        chans = goodKin(4:end);
        color = color_EMG(4,:);
        mycolor = [105,49,158]./255; % another shitty color for speeds
        %mycolor = [50,91,7]./255; % olive for speeds
        colorConf = zeros(48,3);
        for xx=1:3
            colorConf(:,xx)=linspace(mycolor(1,xx),1 ,48);
        end
        colorConf=flip(colorConf);
        %colorConf = color_Extw;
    end
    
    % Get correlations for all days of cortex channels into one matrix
    correlations = [];
    for i = 1:nDays
        % this is for flex (IL and EDL) and ext (MG and FHL) separatly 
        correlations = [correlations; squeeze(corr_AllDaysKin(chans,:,:,i))];
        % this is the mean of all the correlation 
        %correlations = [correlations; squeeze(corr_AllDaysEmg(:,:,:,i))];
    end
    
    meanCorr = squeeze(mean(correlations,1,'omitnan'));
    
    % normalize by the variability of that condition
    if f_norm == 1
        for col = 1:nCond
            meanCorr(:,col) = meanCorr(:,col)/meanCorr(col,col);
        end
    end
    
    if f_plotHist
        % Plot histogram
        figure
        k = 1;
        for cond1 = 1:nCond
            for cond2 = 1:nCond
                corrThis = correlations(:,cond1,cond2);
                nb = calcnbins(corrThis, 'middle');
                sub(k) = subplot(nCond,nCond,k);
                h = histogram(corrThis,nb);
                h.EdgeColor = 'none';
                h.FaceColor = color;
                
                % Add mean
                hold on
                plot([meanCorr(cond1,cond2), meanCorr(cond1,cond2)],[0 max(h.Values)],'LineWidth',2,'color',[0 0 0])
                text(meanCorr(cond1,cond2),max(h.Values)+max(h.Values)*0.12,...
                    mat2str(round(meanCorr(cond1,cond2),2)),'HorizontalAlignment','center')
                
                % Labels, Limits, and Legends
                xlim([-1 1])
                k = k+1;
                if cond2 == 1
                    ylabel(condNames{cond1,1})
                end
                if cond1 == nCond
                    xlabel(condNames{cond2,1})
                end
            end
        end
        set(sub,'box','off');
        set(gcf,'Position',[ 1  1 925 704],'Renderer','Painters')
        clear sub corrThis correlations nb h
    end
    
    % Plot matrix for average
    figure
    %imagesc(tril(meanCorr,-1))
        imagesc(meanCorr)

    colormap(colorConf)
    colorbar
    %     caxis([min(min(meanCorr)) 1])
    caxis([0.2 1])
    % Create strings
    textStrings = num2str(round(meanCorr(:),1));
    textStrings = strtrim(cellstr(textStrings));
    
    % Create x and y coordinates for the strings and plot them
    [x,y] = meshgrid(1:nCond);
    hStrings = text(x(:),y(:),textStrings(:), ...
        'HorizontalAlignment','center');
    
    % Choose white or black for the text color of the strings so
    % they can be easily seen over the background color
    textColors = zeros(length(hStrings),4);
    set(hStrings,{'Color'},num2cell(textColors,2));
    
    % Setting the axis labels
    set(gca,'XTick',1:nCond,...
        'XTickLabel',condNames,...
        'YTick',1:nCond,...
        'YTickLabel',condNames,...
        'TickLength',[0 0]);
    
    set(gcf,'Position',[927 435 354 269],'Renderer','Painters')
    
    
end



%% Plot one channel at a time for inspection
%
% for chan = 49:nChan
%     figure
%     for cond = 1:nCond
%         plot(squeeze(FR_Avg(chan,cond,:)),'color',color_Cyber(cond,:))
%         hold on
%     end
%     title(['Channel ' mat2str(chan)])
%     pause
%     close
% end

%% Correlation in FR population dynamics
for ctx = 1:3 % 1 for Premotor, 2 for Motor, 3 for Sensory
    
    if ctx == 1
        chansCtx = 1:nPM;
        color = color_PM(1,:);
        colorConf = color_PMw;
    elseif ctx == 2
        chansCtx = nPM+1:nPM+nM1;
        color = color_M1(1,:);
        colorConf = color_M1w;
    elseif ctx == 3
        chansCtx = nPM+nM1+1:nPM+nM1+nS1;
        color = color_S1(1,:);
        colorConf = color_S1w;
    end
    
    for f_PC = 1:2 % 1 for flexors, 2 for extensors
        
        if f_PC == 1
            chans = chansCtx(1) + [1,1];
            %             color = color_Grayw(end,:);
            %             colorConf = color_Grayw;
        elseif f_PC == 2
            chans = chansCtx(1)+ [2,2];
            %             color = color_Grayw(end,:);
            %             colorConf = color_Grayw;
        end
        
        % Get correlations for all days of cortex channels into one matrix
        correlations = [];
        for i = 1:nDays
            correlations = [correlations; squeeze(corr_AllDaysPopuFR(chans,:,:,i))];
        end
        
        meanCorr = squeeze(mean(correlations,1,'omitnan'));
        
        % normalize by the variability of that condition
        if f_norm == 1
            for col = 1:nCond
                meanCorr(:,col) = meanCorr(:,col)/meanCorr(col,col);
            end
        end
        
        if f_plotHist
            
            % Plot histogram
            figure
            k = 1;
            for cond1 = 1:nCond
                for cond2 = 1:nCond
                    corrThis = correlations(:,cond1,cond2);
                    nb = calcnbins(corrThis, 'middle');
                    sub(k) = subplot(nCond,nCond,k);
                    h = histogram(corrThis,nb);
                    h.EdgeColor = 'none';
                    h.FaceColor = color;
                    
                    % Add mean
                    hold on
                    plot([meanCorr(cond1,cond2), meanCorr(cond1,cond2)],[0 max(h.Values)],'LineWidth',2,'color',[0 0 0])
                    text(meanCorr(cond1,cond2),max(h.Values)+max(h.Values)*0.12,...
                        mat2str(round(meanCorr(cond1,cond2),2)),'HorizontalAlignment','center')
                    
                    % Labels, Limits, and Legends
                    xlim([-1 1])
                    k = k+1;
                    if cond2 == 1
                        ylabel(condNames{cond1,1})
                    end
                    if cond1 == nCond
                        xlabel(condNames{cond2,1})
                    end
                end
            end
            set(sub,'box','off');
            set(gcf,'Position',[ 1  1 925 704],'Renderer','Painters')
            clear sub corrThis correlations nb h
        end
        
        % Plot matrix for average
        figure('Name',['PC ' mat2str(f_PC)])
        imagesc(tril(meanCorr,-1))
        colormap(colorConf)
        colorbar
        %     caxis([min(min(meanCorr)) 1])
        caxis([0.2 1])
        % Create strings
        textStrings = num2str(round(meanCorr(:),1));
        textStrings = strtrim(cellstr(textStrings));
        
        % Create x and y coordinates for the strings and plot them
        [x,y] = meshgrid(1:nCond);
        hStrings = text(x(:),y(:),textStrings(:), ...
            'HorizontalAlignment','center');
        
        % Choose white or black for the text color of the strings so
        % they can be easily seen over the background color
        textColors = zeros(length(hStrings),4);
        set(hStrings,{'Color'},num2cell(textColors,2));
        
        % Setting the axis labels
        set(gca,'XTick',1:nCond,...
            'XTickLabel',condNames,...
            'YTick',1:nCond,...
            'YTickLabel',condNames,...
            'TickLength',[0 0]);
        
        set(gcf,'Position',[927 435 354 269],'Renderer','Painters')
        
        
    end
end
%% Correlation in EMG population dynamics

for f_PC = 1:2 % 1 for flexors, 2 for extensors
    
    if f_PC == 1
        chans = [1,1];
        color = color_Grayw(end,:);
        colorConf = color_Grayw;
    elseif f_PC == 2
        chans = [2,2];
        color = color_Grayw(end,:);
        colorConf = color_Grayw;
    end
    
    % Get correlations for all days of cortex channels into one matrix
    correlations = [];
    for i = 1:nDays
        correlations = [correlations; squeeze(corr_AllDaysPopuEMG(chans,:,:,i))];
    end
    
    meanCorr = squeeze(mean(correlations,1,'omitnan'));
    
    % normalize by the variability of that condition
    if f_norm == 1
        for col = 1:nCond
            meanCorr(:,col) = meanCorr(:,col)/meanCorr(col,col);
        end
    end
    
    if f_plotHist
        % Plot histogram
        figure
        k = 1;
        for cond1 = 1:nCond
            for cond2 = 1:nCond
                corrThis = correlations(:,cond1,cond2);
                nb = calcnbins(corrThis, 'middle');
                sub(k) = subplot(nCond,nCond,k);
                h = histogram(corrThis,nb);
                h.EdgeColor = 'none';
                h.FaceColor = color;
                
                % Add mean
                hold on
                plot([meanCorr(cond1,cond2), meanCorr(cond1,cond2)],[0 max(h.Values)],'LineWidth',2,'color',[0 0 0])
                text(meanCorr(cond1,cond2),max(h.Values)+max(h.Values)*0.12,...
                    mat2str(round(meanCorr(cond1,cond2),2)),'HorizontalAlignment','center')
                
                % Labels, Limits, and Legends
                xlim([-1 1])
                k = k+1;
                if cond2 == 1
                    ylabel(condNames{cond1,1})
                end
                if cond1 == nCond
                    xlabel(condNames{cond2,1})
                end
            end
        end
        set(sub,'box','off');
        set(gcf,'Position',[ 1  1 925 704],'Renderer','Painters')
        clear sub corrThis correlations nb h
    end
    % Plot matrix for average
    figure('Name',['PC ' mat2str(f_PC)])
    imagesc(tril(meanCorr,-1))
    colormap(colorConf)
    colorbar
    %     caxis([min(min(meanCorr)) 1])
    caxis([0.2 1])
    % Create strings
    textStrings = num2str(round(meanCorr(:),1));
    textStrings = strtrim(cellstr(textStrings));
    
    % Create x and y coordinates for the strings and plot them
    [x,y] = meshgrid(1:nCond);
    hStrings = text(x(:),y(:),textStrings(:), ...
        'HorizontalAlignment','center');
    
    % Choose white or black for the text color of the strings so
    % they can be easily seen over the background color
    textColors = zeros(length(hStrings),4);
    set(hStrings,{'Color'},num2cell(textColors,2));
    
    % Setting the axis labels
    set(gca,'XTick',1:nCond,...
        'XTickLabel',condNames,...
        'YTick',1:nCond,...
        'YTickLabel',condNames,...
        'TickLength',[0 0]);
    
    set(gcf,'Position',[927 435 354 269],'Renderer','Painters')
    
    
end

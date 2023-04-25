% a_PGP_Variance

close all
perc = 85; % This percentage of the most robust channels
f_norm = 1; % Flag to normalize the correlation by the within-task correlation
f_plotDays = 0; % Flag to make interesting plots for individual days
muscles = [2,7,5,8]; % IL, EDL, MG, FHL

%%%%%%%%%%% Multiunit
% load('sessionsElektra.mat')
% dataDir = ['/Users/ismaelseanez/Google Drive/Working Folder/G-Lab/MONKEY/'...
%     'BSI LOCOMOTION Project/Research Material/Analysis/PROCESSED DATA/Elektra/SingleUnit/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% Single Unit
% load('sessionsElektra.mat')
% % dataDir = '/Volumes/SB_ANALYSIS/_ON-GOING_ANALYSIS_DATA/spikesorting/processed/Elektra/Figure2/SingeUnitDs/';

% load('sessionsNatalya.mat')
% dataDir = ['/Users/ismaelseanez/Google Drive/Working Folder/G-Lab/MONKEY/'...
%     'BSI LOCOMOTION Project/Research Material/Analysis/PROCESSED DATA/Natalya/SingleUnit/'];
% 
% load('sessionsNatalya.mat')
Monkey = 'Natalya';
 load(['sessions' Monkey '.mat']);
% dataDir = ['/Users/ismaelseanez/Google Drive/Working Folder/G-Lab/MONKEY/'...
%     'BSI LOCOMOTION Project/Research Material/Analysis/PROCESSED DATA/' Monkey '/SingleUnitClean/'];
% 
%dataDir = ['/Users/Borgo/Google Drive/BSI LOCOMOTION Project/Research Material/Analysis/PROCESSED DATA/Natalya/SingleUnitRealNewNew/'];
%dataDir = ['/Users/Borgo/Google Drive/BSI LOCOMOTION Project/Research Material/Analysis/PROCESSED DATA/' Monkey '/SingleUnitClean/'];
dataDir = ['/Users/Borgo/Documents/MATLAB/fri_mks/BSI LOCOMOTION Project_end_of_july_2021/Research Material/Analysis/PROCESSED DATA/Natalya/SingleUnitClean/'];
f_hope='no';
if regexp(f_hope,'yes')
    load('hope.mat') %load IdxElectOut variable
end

%sessions(5,:) = [];
% sessions(8,:) = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nDays = size(sessions,1);
Results = struct();
colormap_Cyber;
color_ctx = [color_PM(1,:);color_M1(1,:);color_S1(1,:)];
for Day = 1:nDays
    clear DataAllGaits
%     load('/Users/Borgo/Documents/MATLAB/fri_mks/SpikeSorting/Plexon/Figure2/functions/ColorMapDarkHot.mat')
    
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
            
            %Neural_cyclic(:,:,j) = D(j).Neural_cyclic;
            %EMG_cyclic(:,:,j) = D(j).EMG_cyclic;
            DataAllGaits(k).duration = length(D(j).times)/100;
            DataAllGaits(k).FOpercent = D(j).RFO/size(D(j).A,1)*100;
            DataAllGaits(k).Neural_cyclic = D(j).Neural_cyclic;
%             DataAllGaits(k).EMG_cyclic = D(j).EMG_cyclic;
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
        nCtx = 3;
    elseif nChan <160
        nPM = Data.nCellPM; % Number of neurons on each cortex
        if regexp(f_hope,'yes')
            nM1 = Data.nCellM1 - length(cellIDtokickOut);
            nChan = nChan - length(cellIDtokickOut);
        else
            nM1 = Data.nCellM1;
        end
        nS1 = Data.nCellS1;
        nCtx = 3;
    end
    
    % Average setup
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
    FR_All = firingRatesAll; % bsxfun(@minus,firingRatesAll,mean(firingRatesAll,3,'omitnan'));
    FR_Avg = firingRatesAvg; % bsxfun(@minus,firingRatesAvg,mean(firingRatesAvg,3,'omitnan'));
%     EMG_All = bsxfun(@minus,EMGs,mean(EMGs,3,'omitnan'));
%     EMG_Avg = bsxfun(@minus,EMGsAvg,mean(EMGsAvg,3,'omitnan'));

        %% Remove cells from bad electrodes in Natalya
    if regexp(f_hope,'yes')
        cellIDtokickOut2=cellIDtokickOut+nPM;
        FR_Avg(cellIDtokickOut2,:,:)=[];     
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Preferred Gait Phase Analysis Starts Here
    chosenCond = 1; % Choose the standard condition
    
    % Get preferred gait phase for each channel and each gait cycle
    pgpAll = nan(nChan,nCond);
    for cond = 1:nCond
        for chan = 1:nChan
            [~,pgpAll(chan,cond)] = max(FR_Avg(chan,cond,:));
        end
    end
    pgpRad = pgpAll.*2*pi/100; % Convert gait phase to radians
    pgpDeg = rad2deg(pgpRad); % Covert to degrees
      
    % Calculate variance in preferred gait phase across conditions
%     pgpVar = circ_var(pgpRad');
    pgpVar = circ_std(pgpRad');
    
    Results(Day).pgpVar = rad2deg(pgpVar);
    Results(Day).nChans = [nPM,nM1,nS1];
    
end % The end of Days

%% Do analysis and plots for all days together
% Compute total number of neurons
nNeurons = [];
for i = 1:numel(Results)
    nNeurons = [nNeurons; Results(i).nChans];
end

%% Plot barplot for PGPvar of each cortex
figure
pgp_AllDays = zeros(nDays,3);
for ctx = 1:3 % 1 for Premotor, 2 for Motor, 3 for Sensory

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

        temp = Results(Day).pgpVar(chans); % Get results for this day
        pgp_AllDays(Day,ctx) = mean(reshape(temp,[],1));
        
    end
       
end

% Mybar that shit
mybar(pgp_AllDays,'SD',[],'unpaired',{'PM','M1','S1'},'Variance in PGP',color_ctx);
h = gca;
h.TickDir = 'out';
set(gcf,'Position',[360   278   560   420])

[~,pVal_PMM1] = quickStats(pgp_AllDays(:,1),pgp_AllDays(:,2),'paired')
[~,pVal_PMS1] = quickStats(pgp_AllDays(:,1),pgp_AllDays(:,3),'paired')
[~,pVal_M1S1] = quickStats(pgp_AllDays(:,2),pgp_AllDays(:,3),'paired')

avg = mean(pgp_AllDays,1)
sd = std(pgp_AllDays,1)

% a_preferredGaitPhase
% This code will run the analysis for one session on the preferred gait
% phase and change in preferred gait phase.
% Run first the m1_TrainingMatrices file for this monkey which will create
% the average firingRatesAll etc.

clearvars -except diffThreshSD
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

Monkey = 'Elektra';
load(['sessions' Monkey '.mat']);
% dataDir = ['/Users/ismaelseanez/Google Drive/Working Folder/G-Lab/MONKEY/'...
%     'BSI LOCOMOTION Project/Research Material/Analysis/PROCESSED DATA/' Monkey '/SingleUnitClean/'];
 %dataDir = ['/Users/Borgo/Google Drive/BSI LOCOMOTION Project/Research Material/Analysis/PROCESSED DATA/' Monkey '/SingleUnitClean/'];
 dataDir = ['/Users/Borgo/Documents/MATLAB/fri_mks/BSI LOCOMOTION Project_end_of_july_2021/Research Material/Analysis/PROCESSED DATA/' Monkey '/SingleUnitClean/'];

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
        nM1 = Data.nCellM1;
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
    FR_All = firingRatesAll; 
    %FR_ALL = bsxfun(@minus,firingRatesAll,mean(firingRatesAll,3,'omitnan'));
    FR_Avg = firingRatesAvg; 
    %FR_Avg = bsxfun(@minus,firingRatesAvg,mean(firingRatesAvg,3,'omitnan'));
%     EMG_All = bsxfun(@minus,EMGs,mean(EMGs,3,'omitnan'));
%     EMG_Avg = bsxfun(@minus,EMGsAvg,mean(EMGsAvg,3,'omitnan'));
    
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
    
    %% Get the percentage of cells tuned to each gait phase
    
    f_equalWinA = 0; % Flag to make them equal-length windows or not
    
    % Count how many neuorns are tunned to each gait phase
    if f_equalWinA == 1
        regFS = [87.5 12.5]; % Foot Strike
        regStance = [12.5 37.5]; % Early Stance
        regFO = [37.5 62.5]; %  Late stance
        regSwing = [62.5 87.5]; % Swing
    else
        regFS = [100-7.5 7.5]; % 7.5+7.5 = 15%
        regStance = [7.5 50.5];% 50.5-7.5 = 43%
        regFO = [50.5 65.5]; % 65.5-50.5 = 15%
        regSwing = [65.5 92.5]; % 92.5 - 65.5 = 27%
    end
    
    gaitPhasePref = zeros(4,nCtx,nCond);
    
    for cond = 1:nCond
        n1 = find(pgpAll(:,cond)>regFS(1));
        n2 = find(pgpAll(:,cond)<regFS(2));
        regFSchans = cat(1,n1,n2);
        regFS_PM = length(find(regFSchans<=nPM));
        regFS_M1 = length(find(regFSchans>nPM & regFSchans<=nPM+nM1));
        regFS_S1 = length(find(regFSchans>nPM+nM1));
        
        regStanceChans = find(pgpAll(:,cond)>regStance(1) & pgpAll(:,cond)<=regStance(2));
        regStan_PM = length(find(regStanceChans<=nPM));
        regStan_M1 = length(find(regStanceChans>nPM & regStanceChans<=nPM+nM1));
        regStan_S1 = length(find(regStanceChans>nPM+nM1));
        
        regFOChans = find(pgpAll(:,cond)>regFO(1) & pgpAll(:,cond)<regFO(2));
        regFO_PM = length(find(regFOChans<=nPM));
        regFO_M1 = length(find(regFOChans>nPM & regFOChans<=nPM+nM1));
        regFO_S1 = length(find(regFOChans>nPM+nM1));
        
        regSwingChans = find(pgpAll(:,cond)>regSwing(1) & pgpAll(:,cond)<regSwing(2));
        regSwing_PM = length(find(regSwingChans<=nPM));
        regSwing_M1 = length(find(regSwingChans>nPM & regSwingChans<=nPM+nM1));
        regSwing_S1 = length(find(regSwingChans>nPM+nM1));
        
        
        gaitPhasePref(:,:,cond) = [regFS_PM,regFS_M1,regFS_S1;...
            regStan_PM,regStan_M1,regStan_S1;...
            regFO_PM,regFO_M1,regFO_S1;...
            regSwing_PM,regSwing_M1,regSwing_S1];
        
    end
       
    % Normalize by number of units on each cortex or also by window length
    if f_equalWinA == 1
        % Normalize by number of cells on each cortex
        gaitPhasePref = gaitPhasePref./[nPM nM1 nS1]*100;
    else
        winLength = [regFS(2)+(100-regFS(1)), regStance(2)-regStance(1), regFO(2)-regFO(1), regSwing(2)-regSwing(1)];
        gaitPhasePref = gaitPhasePref./winLength';
        
        % As a percentage of total
        a = sum(gaitPhasePref,1);
        b = gaitPhasePref./a*100;
        gaitPhasePref = b;
    end
        
    %% Identify units that change tunning to an "antagonist" phase
    
    f_equalWin = 1; % Flag to make them equal-length windows or not
    
    % Count how many neuorns are tunned to each gait phase
    if f_equalWin == 1
        regFS = [87.5 12.5]; % Foot Strike
        regStance = [12.5 37.5]; % Early Stance
        regFO = [37.5 62.5]; %  Late stance
        regSwing = [62.5 87.5]; % Swing
    else
        regFS = [100-7.5 7.5]; % 7.5+7.5 = 15%
        regStance = [7.5 50.5];% 50.5-7.5 = 43%
        regFO = [50.5 65.5]; % 65.5-50.5 = 15%
        regSwing = [65.5 92.5]; % 92.5 - 65.5 = 27%
    end
    
    neuronWindow = zeros(nChan,nCond);
    
    for cond = 1:nCond
        n1 = find(pgpAll(:,cond)>regFS(1));
        n2 = find(pgpAll(:,cond)<regFS(2));
        
        regFSChans = cat(1,n1,n2); % All channels in this window
        regStanceChans = find(pgpAll(:,cond)>regStance(1) & pgpAll(:,cond)<=regStance(2));
        regFOChans = find(pgpAll(:,cond)>regFO(1) & pgpAll(:,cond)<regFO(2));
        regSwingChans = find(pgpAll(:,cond)>regSwing(1) & pgpAll(:,cond)<regSwing(2));
        
        % Indicate each neuron's tunning phase
        neuronWindow(regFSChans,cond) = 1;
        neuronWindow(regStanceChans,cond) = 2;
        neuronWindow(regFOChans,cond) = 3;
        neuronWindow(regSwingChans,cond) = 4;
        
    end
    
    % Count how many times a neuron changes from 1 to 3, from 2 to 4 (or vice
    % versa), or no "antagonist" change
    countChange = zeros(nChan,3); % 1-(1 to 3); 2-(2-4); 3-(else)
    
    combos = nchoosek(1:nCond,2); % Possible combinations
    
    for i = 1:length(combos)
        
        % Choose tasks for comparison
        cond1 = combos(i,1);
        cond2 = combos(i,2);
        
        for chan = 1:nChan
            if (neuronWindow(chan,cond1) == 1 && neuronWindow(chan,cond2) == 3)...
                    || (neuronWindow(chan,cond1) == 3 && neuronWindow(chan,cond2) == 1)
                countChange(chan,1) = countChange(chan,1)+1;
            elseif (neuronWindow(chan,cond1) == 2 && neuronWindow(chan,cond2) == 4)...
                    || (neuronWindow(chan,cond1) == 4 && neuronWindow(chan,cond2) == 2)
                countChange(chan,2) = countChange(chan,2)+1;
            else
                countChange(chan,3) = countChange(chan,3)+1;
            end
        end
    end
    
    % Convert to percentage
    countChange = countChange./length(combos)*100;
    
    % Count percentage of units that fall in each category
    antagoPM = sum(countChange(1:nPM,1:2))/nPM;
    antagoM1 = sum(countChange(nPM+1:nPM+nM1,1:2))/nM1;
    antagoS1 = sum(countChange(nPM+nM1+1:nPM+nM1+nS1,1:2))/nS1;
    
    antagoChange = [antagoPM;antagoM1;antagoS1];
    
  
    %% change in preferred gait phase for all channels
    pgp_diffAllChans = nan(3,nCond,max([nPM, nM1, nS1]));
    for cond = 1:nCond
        pgp_diffAllChans(1,cond,1:nPM) = abs(angdiff(pgpRad(1:nPM,chosenCond),pgpRad(1:nPM,cond)));
        pgp_diffAllChans(2,cond,1:nM1) = abs(angdiff(pgpRad(nPM+1:nPM+nM1,chosenCond),pgpRad(nPM+1:nPM+nM1,cond)));
        pgp_diffAllChans(3,cond,1:nS1) = abs(angdiff(pgpRad(nPM+nM1+1:nPM+nM1+nS1,chosenCond),pgpRad(nPM+nM1+1:nPM+nM1+nS1,cond)));
    end
    
    %% Make all day plots
    if f_plotDays == 1
        
        %% Plot preferred gait phase distribution
        figure
        for i = 1:3 % 3 cortices
            subplot(3,1,i);
            b = bar(squeeze(gaitPhasePref(:,i,:))','FaceColor','k','EdgeColor','none');
            for j = 2:4
                b(j).FaceAlpha = 1-j*0.2;
                if i == 1
                    
                    ylabel('Premotor')
                elseif i == 2
                    
                    ylabel('Motor')
                elseif i == 3
                    
                    ylabel('Sensory')
                    set(gca,'XTick',1:nCond,'XTickLabel',condNames)
                    
                end
            end
        end
        
        subplot(3,1,1)
        if f_equalWinA == 1
            title('Number of channels tuned to Foot strike, Early stance, Late stance, Swing')
        else
            title('Number of channels tuned to FS, Stance, FO, Swing')
        end
        
        %% Plot the distribution of differences from mean for each cortex and
        % highlight bottom percentile
        
        figure('Name',['Distribution of SD in gait phase differences for ' condNames{chosenCond,1}])
        subplot(1,3,1),histogram(pgp_dMean_sd(1:nPM,chosenCond),10,'FaceColor',color_PM(1,:));
        title('Premotor')
        hold on
        a = gca;
        plot([pgp_bottomPerc(1,chosenCond) pgp_bottomPerc(1,chosenCond)],[0 a.YLim(2)],'k')
        
        subplot(1,3,2),histogram(pgp_dMean_sd(nPM+1:nPM+nM1,chosenCond),10,'FaceColor',color_M1(1,:));
        title('Motor')
        hold on
        a = gca;
        plot([pgp_bottomPerc(2,chosenCond) pgp_bottomPerc(2,chosenCond)],[0 a.YLim(2)],'k')
        
        subplot(1,3,3),histogram(pgp_dMean_sd(nPM+nM1+1:nPM+nM1+nS1,chosenCond),10,'FaceColor',color_S1(1,:));
        title('Sensory')
        hold on
        a = gca;
        plot([pgp_bottomPerc(3,chosenCond) pgp_bottomPerc(3,chosenCond)],[0 a.YLim(2)],'k')
        
        set(gcf,'Position',[238   440   800   198],'Renderer','Painters')
        
        %% Plot the good channels for that condition in Polar coordinates
        cond = chosenCond;
        figure
        for chan = 1:nPM
            if find(chansPM == chan)
                color = [1 0 0];
            else
                color = [0 0 0];
            end
            a = squeeze(pgpRad(chan,cond,1:countCond(cond)-1));
            b = repmat(a,1,2);
            c = ones(size(b));
            c(:,1) = 0;
            subplot(nCond,nCond,chan),polarplot(b',c','Color',color_PM(1,:),'LineWidth',0.07)
            
            rticks({})
            rticklabels({})
            thetaticks({})
            
            % Add average preferred gait phase
            hold on
            polarplot([pgpRad(chan,cond) pgpRad(chan,cond)],[0 1],'color',color,'LineWidth',2)
        end
        set(gcf,'Position',[1    12   408   691],'Renderer','painters')
        
        figure
        for chan = nPM+1:nPM+nM1+1
            if find(chansM1 == chan)
                color = [1 0 0];
            else
                color = [0 0 0];
            end
            a = squeeze(pgpRad(chan,cond,1:countCond(cond)-1));
            b = repmat(a,1,2);
            c = ones(size(b));
            c(:,1) = 0;
            subplot(nCond,nCond,chan-nPM),polarplot(b',c','Color',color_M1(1,:),'LineWidth',0.07)
            
            rticks({})
            rticklabels({})
            thetaticks({})
            
            % Add average preferred gait phase
            hold on
            polarplot([pgpRad(chan,cond) pgpRad(chan,cond)],[0 1],'color',color,'LineWidth',2)
        end
        set(gcf,'Position',[411    14   408   691],'Renderer','painters')
        
        figure
        for chan = nPM+nM1+1:nPM+nM1+nS1
            if find(chansS1 == chan)
                color = [1 0 0];
            else
                color = [0 0 0];
            end
            a = squeeze(pgpRad(chan,cond,1:countCond(cond)-1));
            b = repmat(a,1,2);
            c = ones(size(b));
            c(:,1) = 0;
            subplot(8,8,chan-nPM-nM1),polarplot(b',c','Color',color_S1(1,:),'LineWidth',0.07)
            
            rticks({})
            rticklabels({})
            thetaticks({})
            
            % Add average preferred gait phase
            hold on
            polarplot([pgpRad(chan,cond) pgpRad(chan,cond)],[0 1],'color',color,'LineWidth',2)
        end
        set(gcf,'Position',[820    14   408   691],'Renderer','painters')
        
        %% Plot the difference in preferred
        figure('Name',['Change in preferred gait phase for ' mat2str(perc) '% most robust channels'])
        for ctx = 1:nCtx
            if ctx == 1
                chans = 1:length(chansPM);
                name = 'Premotor';
            elseif ctx == 2
                chans = 1:length(chansM1);
                name = 'M1';
            elseif ctx == 3
                chans = 1:length(chansS1);
                name = 'S1';
            end
            subplot(3,1,ctx)
            
            a = rad2deg(squeeze(pgp_diffGoodChans(ctx,:,chans))');
            [pvalue] = mybar(a,'SD',[0 100],'unpaired',condNames','\Delta in PGP',color_Cyber(1:nCond,:));
            title(name)
        end
            
        set(gcf,'Position',[154   113   856   592],'Renderer','Painters')
        
        
        % Plot the PGP of one neuron for the different conditions
        bestChan = 1; % What number of the good channels to plot
        figure
        for cond = 1:nCond
            a = pgpRad(chansPM(bestChan),cond);
            subplot(3,nCond,cond),polarplot([a a],[0 1],'color',color_PM(1,:),'LineWidth',2)
            
            rticks({})
            rticklabels({})
            thetaticks({})
            
            a = pgpRad(chansM1(bestChan),cond);
            subplot(3,nCond,cond+nCond),polarplot([a a],[0 1],'color',color_M1(1,:),'LineWidth',2)
            
            rticks({})
            rticklabels({})
            thetaticks({})
            
            a = pgpRad(chansS1(bestChan),cond);
            subplot(3,nCond,cond+nCond*2),polarplot([a a],[0 1],'color',color_S1(1,:),'LineWidth',2)
            
            rticks({})
            rticklabels({})
            thetaticks({})
        end
    end % End of plots
    
    % Save all important results
    Results(Day).gaitPhasePref = gaitPhasePref;
    Results(Day).antagoChange = antagoChange;
    Results(Day).pgp_diffAllChans = pgp_diffAllChans;
    Results(Day).nChans = [nPM,nM1,nS1];
    
end % The end of Days

%% Do analysis and plots for all days together
% Compute total number of neurons
nNeurons = [];
for i = 1:numel(Results)
    nNeurons = [nNeurons; Results(i).nChans];
end

% Gait Phase Preference
all_gaitPhasePref = zeros(4,nCtx,nCond,nDays);
for i = 1:numel(Results)
    all_gaitPhasePref(:,:,:,i) = Results(i).gaitPhasePref;
end
avg_gaitPhasePref = mean(all_gaitPhasePref,4);
std_gaitPhasePref = std(all_gaitPhasePref,0,4);

figure('Name','Preferred Gait Phase for each condition')

for i = 1:nCtx % 3 cortices
    subplot(nCtx,1,i);
    model_series = squeeze(avg_gaitPhasePref(:,i,:))';
    model_error = squeeze(std_gaitPhasePref(:,i,:))';
    b = bar(model_series,'grouped','FaceColor','k','EdgeColor','none');
    for j = 2:4
        b(j).FaceAlpha = 1-j*0.2;
        if i == 1

            ylabel('Premotor')
        elseif i == 2

            ylabel('Motor')
        elseif i == 3

            ylabel('Sensory')
            set(gca,'XTick',1:nCond,'XTickLabel',condNames)

        end
    end
    
    % Add errorbars
    hold on
    % Find the number of groups and the number of bars in each group
    ngroups = size(model_series, 1);
    nbars = size(model_series, 2);
    % Calculate the width for each bar group
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    % Set the position of each error bar in the centre of the main bar
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    for i = 1:nbars
        % Calculate center of each bar
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
    end
    hold off
    ylim([0 90])
    set(gca,'TickDir','out')
end

subplot(3,1,1)
if f_equalWin == 1
    title('Number of channels tuned to Foot strike, Early stance, Late stance, Swing')
else
    title('Number of channels tuned to FS, Stance, FO, Swing')
end
 set(gcf,'Renderer','painters')       
%% Antagonist change
all_antagoChange = zeros(nCtx,2,nDays);
for i = 1:numel(Results)
    all_antagoChange(:,:,i) = Results(i).antagoChange;
end

figure('Name','Percentage of units that change to antagonist phase')
[p_TaskInd] = mybar(squeeze(all_antagoChange(1,1:2,:))','SD',[],'unpaired',{'FS-FO','Stance-Swing'},'Percent Cells',color_Cyber(1:2,:));
set(gca,'TickDir','out')
set(gcf,'Position',[509   200   227   196],'Renderer','painters')


%% Change in preferred gait phase for all channels, KEEPING CONDITIONS SEPARATE
nonStandardCond = 1:nCond;
nonStandardCond(chosenCond) = [];
    
figure
for ctx = 1:nCtx
    all_pgpAllChans = [];
    for i = 1:numel(Results)
        temp = squeeze(Results(i).pgp_diffAllChans(ctx,:,:)); % Get results for this day
        removeCol = find(squeeze(isnan(temp(1,:))));  % columns to remove
        temp(:,removeCol) = [];
        all_pgpAllChans = [all_pgpAllChans temp];
    end
    
    subplot(1,nCtx,ctx)
    % Prepare data for histogram
    ratioAll = rad2deg(all_pgpAllChans)';
    binrange = 0:10:180;
    
     %figure(101)
     %cmap=color_Cyber(1:nCond,:);
     %mybar(ratioAll,'SD',[],'unpaired',condNames',[],cmap)
     size(ratioAll)
    
    allCounts = [];
    for cond1 = 1:nCond
        [allCounts(cond1,:)] = histc(ratioAll(:,cond1),binrange);
    end
   
    % get data for stack bar graph
    allCounts(chosenCond,:) = [];
    allCounts = allCounts/sum(sum(allCounts))*100;
    
    %plot stack bargraph
    b = bar(binrange,allCounts','stack','EdgeColor','none');
    hold on  
    for i=1:length(b)
        b(i).FaceColor=color_Cyber(nonStandardCond(i),:);
    end
    
    ylim([0 45])
    % Add text with the total number of neurons for this histogram
    text(120,10,['n = ' mat2str(sum(nNeurons(:,ctx)))])
    set(gca,'TickDir','out')
end
set(gcf,'Position',[258   319   894   212],'Renderer','painters')

%% Change in preferred gait phase for all channels, IGNORING CONDITIONS
nonStandardCond = 1:nCond;
nonStandardCond(chosenCond) = [];
    
figure
for ctx = 1:nCtx
    all_pgpAllChans = [];
    for i = 1:numel(Results)
        temp = squeeze(Results(i).pgp_diffAllChans(ctx,:,:)); % Get results for this day
        removeCol = find(squeeze(isnan(temp(1,:))));  % columns to remove
        temp(:,removeCol) = [];
        all_pgpAllChans = [all_pgpAllChans temp];
    end
    
    subplot(1,nCtx,ctx)
    % Prepare data for histogram
    ratioAll = rad2deg(all_pgpAllChans)';
    ratioAll(:,chosenCond) = [];
    binrange = 0:20:180;
    
    histogram(reshape(ratioAll,1,[]),binrange,'Normalization','probability',...
        'FaceColor',color_ctx(ctx,:),'EdgeColor','none')
 
    meanCounts = mean(reshape(ratioAll,1,[]));
    stdCounts = std(reshape(ratioAll,1,[]));

    % Add mean and std
    hold on
    plot(meanCounts,.68,'o','MarkerFaceColor',color_ctx(ctx,:),'MarkerEdgeColor','none')
    plot([meanCounts-stdCounts, meanCounts+stdCounts],[.68 .68],'color',color_ctx(ctx,:),'LineWidth',0.5)

    % Add text with the total number of neurons for this histogram
    text(120,.30,['n = ' mat2str(sum(nNeurons(:,ctx)))])
    set(gca,'TickDir','out')
    ylim([0 0.7])
    xlim([0 200])
    
end
set(gcf,'Position',[258   319   894   212],'Renderer','painters')

%% Plot barplot for dPGP of each cortex
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

        temp = squeeze(Results(Day).pgp_diffAllChans(ctx,:,:)); % Get results for this day
        removeCol = find(squeeze(isnan(temp(1,:))));  % columns to remove
        temp(:,removeCol) = [];
        temp(chosenCond,:) = [];
        pgp_AllDays(Day,ctx) = median(reshape(temp,[],1));
        
    end
       
end

% Mybar that shit
pgp_AllDays = rad2deg(pgp_AllDays);
mybar(pgp_AllDays,'SD',[0 90],'unpaired',{'PM','M1','S1'},'Median dPGP',color_ctx);
h = gca;
h.TickDir = 'out';
set(gcf,'Position',[360   278   560   420])

[~,pVal_PMM1] = quickStats(pgp_AllDays(:,1),pgp_AllDays(:,2),'paired')
[~,pVal_PMS1] = quickStats(pgp_AllDays(:,1),pgp_AllDays(:,3),'paired')
[~,pVal_M1S1] = quickStats(pgp_AllDays(:,2),pgp_AllDays(:,3),'paired')

avg = mean(pgp_AllDays,1)
sd = std(pgp_AllDays,1)

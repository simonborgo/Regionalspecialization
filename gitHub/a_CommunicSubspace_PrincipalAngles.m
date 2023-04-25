% a_Subspace_BetweenTasks
% Is the output potent manifold between a cortex and EMGs related to
% its task-dependent or task-independet components?
% First, we compute the output potent manifold between each cortex and EMGs
% Then, we compute the principal angles between the output potent manifold of
% a cortex and it's task-dependent components and task-independent
% manifolds
% Finally, we compute the ratio of the principal angles task-independent/task-dependent
% A ratio smaller than 1 - means the output-potent manifold is more similar to
% the task-independent manifold
% A ratio greater than 1 - means the output-potent manifold is more similar to
% the task-dependent manifold

clear all
close all
clc

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

mDimOut = 3; % This is the minimum number of dimensions in output-potent space across Elektra

warning('off','stats:pca:ColRankDefX')

%%%%%%%%%%% Multiunit
load('sessionsElektra')
%load('sessionsNatalya')
% sessions = ['20190415';'20190425';'20190426';'20190430';'20190506';...
%     '20190514';'20190522';'20190531';'20190604'];

% sessions = ['20190415';'20190425'];
%dataDir = ['/Users/ismaelseanez/Google Drive/Working Folder/G-Lab/MONKEY/'...
 %   'BSI LOCOMOTION Project/Research Material/Analysis/PROCESSED DATA/Elektra/MultiUnit/'];
%dataDir = ['C:\Users\physio\Google Drive\BSI LOCOMOTION Project\Research Material\Analysis\PROCESSED DATA\Elektra\MultiUnit\'];
dataDir = ['/Users/Borgo/Documents/MATLAB/fri_mks/BSI LOCOMOTION Project_end_of_july_2021/Research Material/Analysis/PROCESSED DATA/Elektra/MultiunitClean'];
%dataDir = ['/Users/Borgo/Documents/MATLAB/fri_mks/BSI LOCOMOTION Project_end_of_july_2021/Research Material/Analysis/PROCESSED DATA/Natalya/MultiunitClean'];

trainingDir = [dataDir filesep 'TrainingMatricesSimon' filesep]; %TrainingMatrices2 for Nat TrainingMatricesSimon for Elektra
commSpaceDir = [dataDir filesep 'CommunicationSpace' filesep];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% Single Unit
% sessions = ['20190425';'20190425'];
% dataDir = ['/Users/ismaelseanez/Google Drive/Working Folder/G-Lab/MONKEY/'...
%     'BSI LOCOMOTION Project/Research Material/Analysis/PROCESSED DATA/Elektra/SingleUnit/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nDays = size(sessions,1);
Results = struct();
all_prinAngles = zeros(nDays,3,7,7,mDim);
all_prinAnglesSurr = zeros(nDays,3,7,7,mDim);
EMG_prinAngles = zeros(nDays,7,7,mDimMusc);
EMG_prinAnglesSurr = zeros(nDays,7,7,mDimMusc);

%% Load important data and run functions
load('ColorMapDarkHot.mat')
colormap_Cyber;
color_ctx = [color_PM(1,:);color_M1(1,:);color_S1(1,:)];
Yfinal = [];
close all
for myFiles = 1:3
    % Choose the one
    if myFiles == 1
        ctxNames = {'M1'; 'S1'};
        ctxTitles = {'M1toPM';'S1toPM'};
        myCSdata = '_B_pred_fromM1toPMandS1toPM.mat';
        
    elseif myFiles == 2
        ctxNames = {'PM'; 'S1'};
        ctxTitles = {'PMtoM1';'S1toM1'};
        myCSdata = '_B_pred_fromPMtoM1andS1toM1.mat';
        
        
    elseif myFiles ==3
        ctxNames = {'PM'; 'M1'};
        ctxTitles = {'PMtoS1';'M1toS1'};
        myCSdata = '_B_pred_fromPMtoS1andM1toS1.mat';
    end
    
    for Day = 1:nDays
        
        % Load output potent matrices
        load([commSpaceDir 'Elektra_' sessions(Day,:) myCSdata])
        %load([commSpaceDir 'Natalya_' sessions(Day,:) myCSdata])
        nCond = size(B_pred,2);
        
        %% For each cortex
        for ctx = 1:2
            
            % Load training matrix for this cortex
           load([trainingDir 'Training_Mat_Elektra_' sessions(Day,:) '_' ctxNames{ctx,:} '.mat'])
            %load([trainingDir 'Training_Mat_Natalya_' sessions(Day,:) '_' ctxNames{ctx,:} '.mat'])
            
            % Compute principal angle between
            iTI = Matrices.idx_taskInd;
            iTD = Matrices.idx_taskDep;
            
            for cond1 = 1:nCond
                angTI = rad2deg(principal_angles(B_pred{ctx,cond1}(:,1:mDimOut),Matrices.W_All(:,iTI(1:mDimOut))));
                %angTI = princAngle(B_pred{ctx,cond1}(:,1:mDimOut),Matrices.W_All(:,iTI(1:mDimOut)));
                angTD = rad2deg(principal_angles(B_pred{ctx,cond1}(:,1:mDimOut),Matrices.W_All(:,iTD(1:mDimOut))));
                %angTD = princAngle(B_pred{ctx,cond1}(:,1:mDimOut),Matrices.W_All(:,iTD(1:mDimOut)));
                
                OutPot_TI(Day,ctx,cond1,:) = angTI;
                OutPot_TD(Day,ctx,cond1,:) = angTD;
            end
        end
    end % The end of days                 (heyoooo!)
    
    % Calulate averages across days
    avg_OutPot_TI = squeeze(mean(OutPot_TI,1));
    avg_OutPot_TD = squeeze(mean(OutPot_TD,1));
%     %%
%     % Calulate averages across task
%     avg_OutPot_TI = squeeze(mean(OutPot_TI,3));
%     avg_OutPot_TD = squeeze(mean(OutPot_TD,3));
%     dim = 1;
%     comp = nan(7,4); %7days
%     comp(1:7,1) = squeeze(avg_OutPot_TI(:,1,dim)); % CS-TI angle ctx1 vs target
%     comp(1:7,2) = squeeze(avg_OutPot_TI(:,2,dim)); % CS-TI angle ctx2 vs target
%     comp(1:7,3) = squeeze(avg_OutPot_TD(:,1,dim)); % CS-TD angle ctx1 vs target
%     comp(1:7,4) = squeeze(avg_OutPot_TD(:,2,dim)); % CS-TD angle ctx2 vs target
%     
%     % subplot(3,1,myFiles)
%     %mybar(comp,'SD',[],'unpaired',{['CS-TI; ' ctxTitles{1}],['CS-TI; ' ctxTitles{2}], ['CS-TD; ' ctxTitles{1}], ['CS-TD; ' ctxTitles{2}]},'angle(°)',[]);
%     
%     %     cmap{1} = {[.5 .5 .5],[.5 .5 .5],[.5 .5 .5],[.5 .5 .5]};
%     %     cmap{2} = {[.5 .5 .5],[.5 .5 .5],[.5 .5 .5],[.5 .5 .5]};
%     %     cmap{3} = {[.5 .5 .5],[.5 .5 .5],[.5 .5 .5],[.5 .5 .5]};
%     %
%     %
%     %     myscatter(comp,{['CS-TI; ' ctxTitles{1}],['CS-TI; ' ctxTitles{2}], ['CS-TD; ' ctxTitles{1}], ['CS-TD; ' ctxTitles{2}]},cmap{myFiles})
%     %     ylim([50 90])
%     
%     % Plot task-independent and task-dependent angles
%     %figure
%     TIavg = [];
%     TDavg = [];
%     for Day = 1:nDays
%         for ctx = 1:2
%             %subplot(1,2,ctx)
%             for cond1 = 1:nCond
%                 TI = squeeze(OutPot_TI(Day,ctx,cond1,:));
%                 TD = squeeze(OutPot_TD(Day,ctx,cond1,:));
%                 %plot(TI,'color',color_TI)
%                 %hold on
%                 %plot(TD,'color',color_TD)
%                 TIavg = [TIavg TI];
%                 TDavg = [TDavg TD];
%             end
%             
%             title(ctxTitles{ctx,:})
%             if ctx == 2
%                 legend({'Locomotion','Task'},'Location','SouthEast'), legend boxoff
%             elseif ctx == 1
%                 ylabel('Principal angle (deg)')
%                 
%             end
%             xlabel('Neural mode')
%             
%         end
%         
%         
%     end
%     pos = [1,3,5];
%     subplot(3,2,pos(myFiles))
%     plot(mean(TIavg(:,1:35),2),'color',color_TI,'Linewidth',4)
%     hold on
%     plot(mean(TDavg(:,1:35),2),'color',color_TD,'Linewidth',4)
%     ylim([50 90])
%     subplot(3,2,pos(myFiles)+1)
%     plot(mean(TIavg(:,36:70),2),'color',color_TI,'Linewidth',4)
%     hold on
%     plot(mean(TDavg(:,36:70),2),'color',color_TD,'Linewidth',4)
%     ylim([50 90])
%     
%     set(gcf,'Position',[39   143  1197  269],'Renderer','painters')
    %%
    %figure(myFiles)
    for ctx = 1:2
        %sub(ctx) = subplot(1,2,ctx);
        ratioAll  = [];
        clear ratios
        for Day = 1:nDays
            for cond1 = 1:nCond
                for dim = 1:3
                    TI = squeeze(OutPot_TI(Day,ctx,cond1,:));
                    TD = squeeze(OutPot_TD(Day,ctx,cond1,:));
                    %ratios(cond1,dim) = (TI(dim)./TD(dim)-1)*100; % Normalize and convert to percentage
                    ratios(cond1,dim) = (TD(dim)-TI(dim)); % distance
                    ratiosD(Day,cond1,dim,ctx,myFiles) = (TD(dim)-TI(dim)); 
                end
                
            end
            
            ratioAll = [ratioAll ratios];
        end
        ratioAll = ratioAll./max(max(ratioAll));
        
        Y = reshape(ratioAll,[105 1]); %105 if 3 dim, 35 if 1 dim and 70 if 2 dim
%         [y,x]=ksdensity(Y);
%         plot(x,y)
        
        Yfinal = [Yfinal, Y];
        
%         hold on
%         
%         
%         title(ctxTitles{ctx,1})
%         colormap(color_Cyber(1:nCond,:))
%         %xlim([-100 100])
%         % ylim([0 70])
%         xlabel('Similarity ratio (%)')
%         if ctx == 1
%             ylabel('Normalized Count (%)')
%         end
    end
    
    
    
    
end
%%
cnt = 1;
cmap(:,1) = color_M1(1,:)';
cmap(:,2) = color_S1(1,:)';
cmap(:,3) = color_PM(1,:)';
cmap(:,4) = color_S1(1,:)';
cmap(:,5) = color_PM(1,:)';
cmap(:,6) = color_M1(1,:)';
for ctx = 1:3 
    subplot(1,4,ctx)
    for j = 1:2
    [y,x] = ksdensity(Yfinal(:,cnt));
    [~,pvalue] = ttest(Yfinal(:,cnt))
    %plot(x,y)
    a = area(x,y,'FaceColor',cmap(:,cnt),...
        'EdgeColor',cmap(:,cnt));
    a.FaceAlpha = 0.1;
    hold on 
    plot([0 0],[0 1.5],'-','Color','k')
    cnt = cnt+1;
    end
    xlim([-1.5 1.5])
    ylim([0 2])
    %camroll(270)
end
set(gcf,'Position',[39   143  1197  269],'Renderer','painters')

% figure
% mybar_gray_ND(Yfinal,'SEM',[],'unpaired',{'m1-pm','s1-pm','pm-m1','s1-m1','pm-s1','m1-s1'},'fuck',[])
% %
% figure
Yfinal2(1:210,1) = [Yfinal(:,1); Yfinal(:,2)];
Yfinal2(1:210,2) = [Yfinal(:,3); Yfinal(:,4)];
Yfinal2(1:210,3) = [Yfinal(:,5); Yfinal(:,6)];
% mybar_gray_ND(Yfinal2,'SEM',[],'unpaired',{'to PM','to M1','to S1'},'fuck',[])
% %
% % mega dirty: don't like it because not same unit than histo so
% size(ratiosD) %ratiosD = Day,cond1,dim,ctx,myFiles
% a = squeeze(mean(ratiosD,2));
% a = squeeze(mean(a,2));
% a = squeeze(mean(a,2));
% figure
% mybar(a,'SEM',[],'unpaired',{'to PM','to M1','to S1'},'fuck',[])
cnt = 1;
cnt2 = 30;
clear a
for i = 1:7
    a(i,:) = mean(Yfinal2(cnt:cnt2,:),1);
    cnt = cnt + 30;
    cnt2 = cnt2 + 30;
end
subplot(1,4,4)
mybar(a,'SEM',[-0.1 0.5],'unpaired',{'to PM','to M1','to S1'},'alignment index',color_ctx)

pause


%% quick stat
%mybar([a(:,1) a(:,3)],'SEM',[],'paired',{'to PM','to S1'},'alignment index',[])




%% Plot effect size
effectSize = OutPot_TD-OutPot_TI;
meanEffect2 = squeeze(mean(effectSize,3));
meanEffect = squeeze(mean(meanEffect2,3));

[p,h] = signrank(meanEffect(:,1))
[p,h] = signrank(meanEffect(:,2))
figure
axis('square')
pVal = mybar([squeeze(meanEffect(:,1)), squeeze(meanEffect(:,2))],'SD',[],'unpaired',{'PM','M1'},...
    'Diff. in principal angle (TD-TI)',...
        [color_PM(1,:); color_M1(1,:); color_S1(1,:)])
set(gca,'TickDir','out')
set(gcf,'Position',[456   250   241   261],'Renderer','painters')

ylim([-10,30])
axis('square')

%% Plot task-independent and task-dependent angles
figure
TIavg = [];
TDavg = [];
for Day = 1:nDays
    for ctx = 1:2
        subplot(1,2,ctx)
        for cond1 = 1:nCond
            TI = squeeze(OutPot_TI(Day,ctx,cond1,:));
            TD = squeeze(OutPot_TD(Day,ctx,cond1,:));
%             plot(TI,'color',color_TI) 
%             hold on
%             plot(TD,'color',color_TD)
            TIavg = [TIavg TI];
            TDavg = [TDavg TD];
        end
        
        title(ctxTitles{ctx,:})
        if ctx == 2
            legend({'Locomotion','Task'},'Location','SouthEast'), legend boxoff
        elseif ctx == 1
            ylabel('Principal angle (deg)')
            
        end
        xlabel('Neural mode')

    end
    
    
end
subplot(1,2,1)
        plot(mean(TIavg(:,1:35),2),'color',color_TI,'Linewidth',4)
        hold on
        plot(mean(TDavg(:,1:35),2),'color',color_TD,'Linewidth',4)
        subplot(1,2,2)
        plot(mean(TIavg(:,36:70),2),'color',color_TI,'Linewidth',4)
         plot(mean(TDavg(:,36:70),2),'color',color_TD,'Linewidth',4)

set(gcf,'Position',[39   143  1197  269],'Renderer','painters')

%% Plot all ratios together all days
figure
for Day = 1:nDays
    for ctx = 1:2
        subplot(2,2,ctx)
        for cond1 = 1:nCond
            TI = squeeze(OutPot_TI(Day,ctx,cond1,:));
            TD = squeeze(OutPot_TD(Day,ctx,cond1,:));
            plot(TI./TD,'color',color_Cyber(cond1,:))
            hold on
        end
        
        title(ctxTitles{ctx,:})
        ylim([0 2])
        plot([1 mDimOut],[1 1],'--k')
        if ctx == 1
            ylabel('TI/TD ratio')
        end
        xlabel('Neural mode')
    end
    
    
end

%% Plot average ratio for each condition
for ctx = 1:2
    subplot(2,2,ctx+2)
    for cond1 = 1:nCond
        TI = squeeze(avg_OutPot_TI(ctx,cond1,:));
        TD = squeeze(avg_OutPot_TD(ctx,cond1,:));
        plot(TI./TD,'color',color_Cyber(cond1,:))
        hold on
    end
    title(ctxTitles{ctx,:})
    ylim([0 2])
    plot([1 mDimOut],[1 1],'--k')
    
    if ctx == 1
            ylabel('TI/TD average ratio')
    end
    xlabel('Neural mode')
    
end


set(gcf,'Position',[79   151   1173    516],'Renderer','painters')
%%

%h = myhist(x,'sim','cnt',[.5 .5 .5]');
%plot(conv(h.BinEdges, [0.5 0.5], 'valid'), h.BinCounts)

%%
x = rand(1,100)
h = histogram(x);
hold on
plot(conv(h.BinEdges, [[0.5 0.5]], 'valid'), h.BinCounts)
%%
Y=randn(1000,1)
hist(Y)
[y,x]=hist(Y); hold on
plot(x,y,'r')
%%
[y,x]=ksdensity(Y)
plot(x,y)
%% Plot histograms: my way
figure
for ctx = 1:2
    sub(ctx) = subplot(1,2,ctx);
    ratioAll  = [];
    clear ratios
    for Day = 1:nDays
        for cond1 = 1:nCond
            for dim = 1:3
                TI = squeeze(OutPot_TI(Day,ctx,cond1,:));
                TD = squeeze(OutPot_TD(Day,ctx,cond1,:));
                %ratios(cond1,:) = (TI./TD-1)*100; % Normalize and convert to percentage
                ratios(cond1,dim) = (TD(dim)-TI(dim)); % distance
                %hold on
            end
            %ratios = zscore(ratios,[],2);
        end
        
        ratioAll = [ratioAll ratios];
    end
    %ratioAll = zscore(ratioAll,[],2);
    ratioAll = ratioAll./max(max(ratioAll));
    %ratioAll = zscore(ratioAll);
    % Do a histogram separating by condition
    %binrange = -100:10:100;
    %     binrange = -1:0.1:1;
    %     allCounts = [];
    %     for cond1 = 1:nCond
    %         [allCounts(cond1,:)] = histc(ratioAll(cond1,:),binrange);
    %     end
    %
    %     allCounts = allCounts/sum(sum(allCounts))*100;
    %     b1 = bar(binrange,allCounts','stack','EdgeColor','none','FaceColor','Flat');
    %     for i = 1:5
    %       b1(i).CData = repmat(color_Cyber(i,:),size(b1(i).CData,1),1);
    %     end
    % hold on
    
    Y = reshape(ratioAll,[105 1]);
    [y,x]=ksdensity(Y);
    plot(x,y)
    
    
    title(ctxTitles{ctx,1})
    colormap(color_Cyber(1:nCond,:))
    %xlim([-100 100])
    % ylim([0 70])
    xlabel('Similarity ratio (%)')
    if ctx == 1
        ylabel('Normalized Count (%)')
    end
    
     %% Add stats!
   a = reshape(allCounts,1,[]);
   [~,pVal] = ttest(a); % Compare population to 0
    if pVal >= 0.05
        pValtxt = '';
    elseif (pVal < 0.05 && pVal >= 0.005)
        pValtxt = '*';
    elseif (pVal < 0.005 && pVal >= 0.0005)
        pValtxt = '**';
    elseif (pVal < 0.0005)
        pValtxt = '***';
    end

    [yPos,xPos] = max(sum(allCounts,1));
    text(binrange(xPos),yPos+5,pValtxt,'HorizontalAlignment','center');
end

linkaxes(sub,'y');
set(gcf,'Position',[35         156        1225         305],'Renderer','painters')



%% Plot histograms: Nicolo's way
figure
for ctx = 1:2
    sub(ctx) = subplot(1,2,ctx);
    ratioAll  = [];
    for Day = 1:nDays
        for cond1 = 1:nCond
            TI = squeeze(OutPot_TI(Day,ctx,cond1,:));
            TD = squeeze(OutPot_TD(Day,ctx,cond1,:));
            ratios(cond1,:) = (TI./TD-1)*100; % Normalize and convert to percentage
            hold on
        end
        
        ratioAll = [ratioAll ratios];  
    end
    
    % Do a histogram separating by condition
    binrange = -100:10:100;
    allCounts = [];
    for cond1 = 1:nCond
        [allCounts(cond1,:)] = histc(ratioAll(cond1,:),binrange);
    end

%         allCounts = allCounts/sum(sum(allCounts))*100;
    allCounts = allCounts/sum(sum(allCounts))*100;
    b1 = bar(binrange,allCounts','stack','EdgeColor','none','FaceColor','Flat');
    for i = 1:5
      b1(i).CData = repmat(color_Cyber(i,:),size(b1(i).CData,1),1);
    end
    hold on  


    title(ctxTitles{ctx,1})
    colormap(color_Cyber(1:nCond,:))
    xlim([-100 100])
    ylim([0 70])
    xlabel('Similarity ratio (%)')
    if ctx == 1
        ylabel('Normalized Count (%)')
    end
    
     % Add stats!
   a = reshape(allCounts,1,[]);
   [~,pVal] = ttest(a); % Compare population to 0
    if pVal >= 0.05
        pValtxt = '';
    elseif (pVal < 0.05 && pVal >= 0.005)
        pValtxt = '*';
    elseif (pVal < 0.005 && pVal >= 0.0005)
        pValtxt = '**';
    elseif (pVal < 0.0005)
        pValtxt = '***';
    end

    [yPos,xPos] = max(sum(allCounts,1));
    text(binrange(xPos),yPos+5,pValtxt,'HorizontalAlignment','center');
end

linkaxes(sub,'y');
set(gcf,'Position',[35         156        1225         305],'Renderer','painters')

%% Plot histograms for Task-Dependent and task-Independent
figure
for ctx = 1:2
    a = reshape(OutPot_TI(:,ctx,:,1),1,[]);
    b = reshape(OutPot_TD(:,ctx,:,1),1,[]);
    nba = calcnbins(a, 'middle');
    nbb = calcnbins(b, 'middle');
    sub(ctx) = subplot(1,2,ctx);
    hb = histogram(b,nbb,'Normalization','probability');
    hb.FaceColor = color_TD;
    hb.EdgeColor = 'none';
    hb.FaceAlpha = 0.6;
    hold on
    ha = histogram(a,nba,'Normalization','probability');
    ha.FaceColor = color_TI;
    ha.EdgeColor = 'none';
    ha.FaceAlpha = 0.6;
    title(ctxTitles{ctx,1})
    xlim([0 100])
    set(gca,'TickDir','out');
    xlabel('Principal angle (deg)')
    
    % Add stats
    maxY = max([ha.Values hb.Values]);
    
%     plot(mean(a),maxY*0.10+maxY,'o','MarkerSize',8,...
%         'MarkerFaceColor',color_TI,'MarkerEdgeColor','none')
%     plot([mean(a)-std(a) mean(a)+std(a)],...
%          [maxY*0.10+maxY maxY*0.10+maxY],'-','color',color_TI,'LineWidth',1)
%     plot(mean(b),maxY*0.15+maxY,'o','MarkerSize',8,...
%         'MarkerFaceColor',color_TD,'MarkerEdgeColor','none')
%     plot([mean(b)-std(b) mean(b)+std(b)],...
%          [maxY*0.15+maxY maxY*0.15+maxY],'-','color',color_TD,'LineWidth',1)
%      
%     [h,pVal] = quickStats(a, b,'unpaired');
%     if h == 1
%        sigstar([mean(a) mean(b)],pVal)
%     end
%     
    if ctx == 1
        ylabel('Normalized count')
    elseif ctx == 2
        legend([ha hb],{'task-indep','task-dep'},'Location','north','box','off')
    end
  
end
linkaxes(sub,'y')
set(gcf,'Position',[1         163        1255         322],'Renderer','painters')
%% Plot histograms for Task-Dependent and task-Independent
figure(1001)
for ctx = 1:2
subplot(1,2,ctx)
a = reshape(OutPot_TI(:,ctx,:,1),1,[]);
    b = reshape(OutPot_TD(:,ctx,:,1),1,[]);
    
    mybar_gray_ND([a;b]','SD',[0 90],'unpaired',{'task-indep','task-dep'},'angles (°)',[color_TI;color_TD]);
    

    title(ctxTitles{ctx,1})
    set(gca,'TickDir','out');

   
end
set(gcf,'Position',[1         163        1255         322],'Renderer','painters')

% load handel
% sound(y,Fs)
%%% This script plots crosscorrelation coefficients and lags for the
%%% specified condition and frequency band. It also plots a histogram of the
%%% same values.
%%%
%%% Osorio & Assaneo, 2025

% initialize BST
cd('F:\Matlab\brainstorm3');
brainstorm
%%
clear, clc, close,

condition2analyze = 'Speech';   % 'speech' or 'music'
band2analyze      = 'HFB';      % SFB (1-8 Hz) or HFB (70-120 Hz) 
plot_oncortex     = 1;          % whether to plot effect on cortical surface
plot_histogram    = 0;          % whether to plot histogram of plotted values
effect2plot       = 'lag';      % rho (corrcoefficients) or lag (xcorr lags)

% set paths
project_dir = 'PATH';
data_dir = [project_dir,filesep,'Data'];

% load data to plot
load([data_dir,filesep,'CROSdata_',band2analyze,'_whitenoise.mat']);

% these are the parameters after optimization per condition and freq band
if strcmpi(condition2analyze,'speech') && strcmpi(band2analyze,'SFB')
    mindist = 0.006; minpoints = 16;
elseif strcmpi(condition2analyze,'speech') && strcmpi(band2analyze,'HFB')
    mindist = 0.006; minpoints = 12;
elseif strcmpi(condition2analyze,'music') && strcmpi(band2analyze,'SFB')
    mindist = 0.012; minpoints = 12;
elseif strcmpi(condition2analyze,'music') && strcmpi(band2analyze,'HFB')
    mindist = 0.014; minpoints = 12;
end

% locate data in separate variables and plot
if strcmpi(effect2plot,'rho')
    effect4speech = dataMat(:,:,1);
    effect4music  = dataMat(:,:,2);
elseif strcmpi(effect2plot,'lag')
    effect4speech = LagMat(:,:,1);
    effect4music  = LagMat(:,:,2);
end

% get MNI cortical surface using brainstorm. We will plot data here.
SurfaceFile   = 'F:\MATLAB\brainstorm_db\iEEG\anat\@default_subject\tess_cortex_pial_low.mat';

% create arrays containing the data of interest
if strcmpi(condition2analyze,'speech')
    counter = 1;
    for sub_i=1:length(sub2plot)    
        clear ThisSubStruct
        ThisSubStruct = cell2struct(AllChannelLabels{sub_i},names4fields,2);        
        for idx=1:length(dataMat)
            if ~isnan(effect4speech(idx,sub_i)) %% && isnan(rhos4music(idx,sub_i)) %%
                testMat(counter,:)   = [ThisSubStruct(idx).Loc(1),ThisSubStruct(idx).Loc(2),ThisSubStruct(idx).Loc(3)];
                ValRange(counter,:)  = effect4speech(idx,sub_i);
                subIDelec(counter,:) = [sub_i,idx];
                counter = counter + 1;
            end
        end
    end    
elseif strcmpi(condition2analyze,'music')    
    counter = 1;    
    for sub_i=1:length(sub2plot)       
        clear ThisSubStruct
        ThisSubStruct = cell2struct(AllChannelLabels{sub_i},names4fields,2);        
        for idx=1:length(dataMat)
            if  ~isnan(effect4music(idx,sub_i)) %% isnan(rhos4speech(idx,sub_i)) &&
                testMat(counter,:)   = [ThisSubStruct(idx).Loc(1),ThisSubStruct(idx).Loc(2),ThisSubStruct(idx).Loc(3)];
                ValRange(counter,:)  = effect4music(idx,sub_i);
                subIDelec(counter,:) = [sub_i,idx];
                counter = counter + 1;
            end
        end
    end 
end

% now get only the desired electrodes that belong to the specified cluster
testClusters = dbscan(testMat,mindist,minpoints);
tmptable     = sortrows(tabulate(testClusters),-3);
deleteThis   = 0;

% identify clusters where number of elecs < minpoints
if any(tmptable(:,2) < minpoints)
    deleteThis = tmptable(tmptable(:,2) < minpoints,1);
end

% group smaller surviving clusters into one single cluster
if size(tmptable,1)-1 > 3
    testClusters(ismember(testClusters,tmptable(4:end,1))) = 3;
end

% get rid of non-clustered electrodes
ValRange(ismember(testClusters,[-1; deleteThis]))     = [];
testMat(ismember(testClusters,[-1; deleteThis]),:,:)  = [];
subIDelec(ismember(testClusters,[-1; deleteThis]),:)  = [];
testClusters(ismember(testClusters,[-1; deleteThis])) = []; 

% get mean rho vals at the subject level
if strcmpi(condition2analyze,'music')
    for i=1:30
        meanRho(i)  = mean(dataMat(subIDelec(subIDelec(:,1) == i,2),i,2));
        meanPval(i) = mean(PValMat(subIDelec(subIDelec(:,1) == i,2),i,2));
        meanLag(i)  = mean(LagMat(subIDelec(subIDelec(:,1) == i,2),i,2));
    end
else 
    for i=1:30
        meanRho(i)  = mean(dataMat(subIDelec(subIDelec(:,1) == i,2),i,1));
        meanPval(i) = mean(PValMat(subIDelec(subIDelec(:,1) == i,2),i,1));
        meanLag(i)  = mean(LagMat(subIDelec(subIDelec(:,1) == i,2),i,1));        
    end
end

% save data
save([data_dir,filesep,'stats_clustered_',band2analyze,'_',condition2analyze,'.mat'],'meanRho','meanPval','meanLag');

% convert data from samples to ms
if strcmpi(effect2plot,'Lag')
    ValRange = ValRange/200;
end

% plot statistical effect on cortex
if plot_oncortex == 1
    [hFig, iDS, iFig] = view_surface(SurfaceFile);
    hFig.Color = [1 1 1];
    hold on;
    
    sh  = scatter3(testMat(:,1),testMat(:,2),testMat(:,3),60,ValRange,'filled');
    sh.MarkerFaceAlpha = .8;
    sh.SizeData = 80;
    ch = colorbar;
    ch.FontSize = 15;
    if strcmpi(effect2plot,'rho') && strcmpi(condition2analyze,'Speech')
        set(gca,'clim',[0.08 0.12])  
        colormap hot
        hFig.Colormap =  hFig.Colormap(6:length(hFig.Colormap)-6,:);
    elseif strcmpi(effect2plot,'rho') && strcmpi(condition2analyze,'Music')
        set(gca,'clim',[0.08 0.10])
        colormap hot
        hFig.Colormap =  hFig.Colormap(6:length(hFig.Colormap)-6,:);
    else
        set(gca,'clim',[-1.5 1.5])
        hFig.Colormap(22:64-21,:) = [[linspace(hFig.Colormap(21,1),1,11)', linspace(hFig.Colormap(21,2),1,11)', linspace(hFig.Colormap(21,3),1,11)']; ...
                                     [linspace(1,hFig.Colormap(64-21,1),11)', linspace(1,hFig.Colormap(64-21,2),11)', linspace(1,hFig.Colormap(64-21,3),11)']];
    end
end

% histogram for ploted values
if plot_histogram == 1
    % colors per condition (1,:) music, (2,:) speech, (3,:) music and speech
    colors 	= [0.7176 0.2745 1.0000; ...
               0.9412 0.5804 0.3373; ...
               0.4667 0.6745 0.1882];  
           
    figure,clf
    hh = histogram(ValRange,8); hold on;
    plot([1 1]*mean(ValRange),get(gca,'ylim'),'k','linew',1)
    plot([1 1]*median(ValRange),get(gca,'ylim'),'k--','linew',1)
    legend({'data','mean','median'}, 'Box','off')
    if strcmpi(effect2plot,'rho')
        ylabel('Electrode count'); %ylim([0 30]);
        xlabel('Correlation coefficient'); xlim([0 0.3]);
    else
        ylabel('Electrode count'); %ylim([0 50]);
        xlabel('Lag (ms)'); xlim([-2 2]);
    end
    if strcmpi(condition2analyze,'music')
        hh.FaceColor = colors(1,:);
    else
        hh.FaceColor = colors(2,:);
    end       
    hh.EdgeColor = [1 1 1];
    set(gca,'FontSize',12,'FontName','Arial');
    box off
end

% get some informative data
if strcmpi(condition2analyze,'music')
    pvals = PValMat(:,:,2);
    if strcmpi(effect2plot,'rho')
        caxis([0.08 0.10])
        ch.Ticks = 0.08:0.005:0.10;
        ch.Label.String = 'Rho';
    else
        ch.Label.String = 'Lag (s)';
        caxis([-1.5 1.5]);
        ch.Ticks = -1.5:.5:1.5;
    end
else
    pvals = PValMat(:,:,1);
    if strcmpi(effect2plot,'rho')
        caxis([0.08 0.12])
        ch.Ticks = 0.08:0.01:0.12;
    else
        ch.Label.String = 'Lag (s)';
        caxis([-1.5 1.5]);
        ch.Ticks = -1.5:.5:1.5;
    end
end

pvals = pvals(~isnan(pvals));

% print some descriptive statistics
if strcmpi(effect2plot,'rho')
    disp(['n = ' num2str(sum(~isnan(meanRho))), ...
          ', mean = ' num2str(mean(ValRange),2), ...
          ', median = ' num2str(median(ValRange),2), ...
          ', elecs = ' num2str(length(ValRange)), ...
          ', p = ' num2str(mean(pvals),2), ...
          ', SD = ' num2str(std(ValRange),2)]);
else
    disp(['n = ' num2str(sum(~isnan(meanRho))), ...
          ', mean = ' num2str(mean(ValRange),2), ...
          ', median = ' num2str(median(ValRange),2), ...
          ', elecs = ' num2str(length(ValRange)), ...
          ', p = ' num2str(mean(pvals),2), ...
          ', SD = ' num2str(std(ValRange),2)]);
end

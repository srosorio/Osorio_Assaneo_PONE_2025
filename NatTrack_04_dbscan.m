%%% This script plots statistically significant electrodes in common space
%%% (MNI) per condition, and conducts a data-driven cluster analysis to identify
%%% cortical regions driving statistical effects. This code uses the
%%% Brainstorm toolbox. 
%%%
%%% Osorio & Assaneo, 2025

% initialize brainstorm (for data visualization)
cd('PATH\brainstorm3');
brainstorm
%%
clear, clc, close

% -------------------------------------------------------------------------
%                             SET PARAMETERS
% -------------------------------------------------------------------------
plot_neteffect    = 1;  % plot effect prior to statistics
dbscan_param      = 1;  % plot dbscan optimization process
condition2analyze = 'Music';
band2analyze      = 'SFB';
% -------------------------------------------------------------------------

% set paths
project_dir = 'PATH';
data_dir = [project_dir,filesep,'Data'];

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

% ------------------- nothing to change from here on ----------------------

% colors per condition (1,:) music, (2,:) speech, (3,:) music and speech
colors  = [0.7176 0.2745 1.0000; ...
           0.9412 0.5804 0.3373];              
    
% load data to plot
load([data_dir,filesep,'CROSdata_',band2analyze,'_whitenoise.mat']);

% locate data in separate variables and plot
rhos4speech = dataMat(:,:,1);
rhos4music  = dataMat(:,:,2);

% get MNI cortical surface using brainstorm. We will plot data here.
SurfaceFile   = 'PATH\anat\@default_subject\tess_cortex_pial_low.mat';

n_subs  = length(sub2plot);
n_elecs = length(dataMat);

% first, we get the info we need and put it in a single array 
if strcmpi(condition2analyze,'speech')
    counter = 1;
    for sub_i=1:n_subs    
        % get the subject structure with electrode info
        clear ThisSubStruct
        ThisSubStruct = cell2struct(AllChannelLabels{sub_i},names4fields,2);    
        for idx=1:n_elecs
            % electrodes selective to speech 
            if ~isnan(rhos4speech(idx,sub_i))
                datamat(counter,:)   = [ThisSubStruct(idx).Loc(1),ThisSubStruct(idx).Loc(2),ThisSubStruct(idx).Loc(3)];
                ValRange(counter,:)  = rhos4speech(idx,sub_i);
                subIDelec(counter,:) = [sub_i,idx];
                counter = counter + 1;
            end
        end
    end    
elseif strcmpi(condition2analyze,'music')    
    counter = 1;   
    for sub_i=1:n_subs       
        clear ThisSubStruct
        ThisSubStruct = cell2struct(AllChannelLabels{sub_i},names4fields,2);
        for idx=1:n_elecs
            % electrodes selective to music
            if  ~isnan(rhos4music(idx,sub_i))
                datamat(counter,:)   = [ThisSubStruct(idx).Loc(1),ThisSubStruct(idx).Loc(2),ThisSubStruct(idx).Loc(3)];
                ValRange(counter,:)  = rhos4music(idx,sub_i);
                subIDelec(counter,:) = [sub_i,idx];
                counter = counter + 1;
            end
        end
    end 
%     save([data_dir,filesep,'dbscan_results_' condition2analyze '_' band2analyze '.mat'],'datamat','ValRange','subIDelec')
end

% plot all electrodes (Net effect, without cluster analysis)
if plot_neteffect == 1
    % bst will throw an error every time the cortical surface is ploted for
    % the first time. This try catch statement overrides this error. Close
    % the loading bar if it doesn't disappear on its own. 
    try
        [hFig, iDS, iFig] = view_surface(SurfaceFile);
    catch
        [hFig, iDS, iFig] = view_surface(SurfaceFile);
    end
    hFig.Color = [1 1 1];
    hold on;
    
    % plot electrodes
    sh  = scatter3(datamat(:,1),datamat(:,2),datamat(:,3),60,'filled');
    sh.MarkerFaceAlpha = .8;
    sh.SizeData = 100;
    % set colors according to condition
    if strcmpi(condition2analyze,'music')
        sh.CData = colors(1,:);
    elseif strcmpi(condition2analyze,'speech')
        sh.CData = colors(2,:);
    elseif strcmpi(condition2analyze,'both')
        sh.CData = colors(3,:);
    end
    sh.MarkerFaceAlpha = .5;
    sh.MarkerEdgeColor  = [0 0 0];
    sh.MarkerEdgeAlpha = .5;
    sh.SizeData = 50;
end

% this is to determine the best parameters for dbscan algorithm. We are
% looking for the highest possible value in the y axis and the lowest
% possible value in the x axis
if dbscan_param == 1
    clear testClusters prctelecs numclusters numoutliers
    dsts   = 0:0.002:0.04;
    points = 12:2:20;

    for idx=1:length(dsts)
        clear clusterdata
        for jdx=1:length(points)
            tmp = dbscan(datamat,dsts(idx),points(jdx));
            clusterdata = tabulate(tmp);
            if numel(clusterdata(:,1)) == 1 && (clusterdata(1,1) == -1)
                prctelecs(jdx,idx)   = 0;
                prctoutliers(jdx,idx) = clusterdata(1,3);
            elseif numel(clusterdata(:,1)) > 1 && (clusterdata(1,1) == -1)
                prctelecs(jdx,idx)   = sum(clusterdata(2:end,3));
                numclusters(jdx,idx) = sum(clusterdata(:,1) ~= -1);
                prctoutliers(jdx,idx) = clusterdata(1,3);
            elseif numel(clusterdata(:,1)) == 1 && (clusterdata(1,1) == 1)
                prctelecs(jdx,idx)   = sum(clusterdata(:,3));
                numclusters(jdx,idx) = sum(clusterdata(:,1) ~= -1);
                prctoutliers(jdx,idx) = 0;
            end
        end
    end
    
    opt_idx = prctelecs ./ (dsts*1000) .* points' .* numclusters;  
    max_val = max(max(opt_idx));
    opt_idx = opt_idx / max_val; 
    
    % plot data
    figure(3), clf
    plot(dsts,opt_idx,'LineWidth',1.5)
    ylabel('OI_{norm}');
    xlabel('epsilon');
    legend(num2str(points'),'Location','SouthEast'); legend boxoff
    ylim([0 1.02]); set(gca,'FontSize',20);
    title([condition2analyze ' - ' band2analyze],'FontWeight','normal')
    box off
end
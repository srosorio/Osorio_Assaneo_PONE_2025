%%% This code creates a csv file with all the data necessary for LMER
%%% analyses of the single models  (Music,  Speech). The output file is to 
%%% be used as the input dataset for the R code used to conduct our 
%%% statistical analyses.
%%% 
%%% Osorio & Assaneo, 2025

clear, clc

% set paths
project_dir = 'PATH';
data_dir    = [project_dir,filesep,'Data'];
cd(data_dir)

fnames = dir('Electrode_AnatLabels*');
print_summary = 1;  % whether to print a table with the total numbers per region

if print_summary == 1
    for idx=1:length(fnames)
        load(fnames(idx).name)
        AnatLabels = cellstr(AnatLabels);
        info = split(fnames(idx).name,'_');
        cond = info{3};
        band = split(info{4},'.mat');
        band = band{1};
        disp([cond, '-' ,band]);
        tbl = sortrows(tabulate(cellstr(AnatLabels)),3,'descend');
        if idx == 1
           tmptbl = tbl;
        else 
            tmptbl = [tmptbl; tbl];
        end
        disp(tbl);
    end
end

% create table with the number of electrodes per anatomical region (total,
% across conditions)
regions = unique(tmptbl(:,1));
for idx=1:length(regions)
    this_region = strcmpi(regions{idx},cellstr(tmptbl(:,1)));
    num_elec    = sum(cell2mat(tmptbl(this_region,2)));
    datacell{idx,1} = regions{idx};
    datacell{idx,2} = num_elec;
end

% writetable(cell2table(sortrows(datacell,2,'descend'),'VariableNames',{'Region','Electrodes'}),[data_dir,filesep,'ElectrodeCountPerAnatRegions.xlsx'])

condition2analyze = {'Speech','Music'};   % 'speech' or 'music'
band2analyze      = {'SFB','HFB'};      % SFB (1-8 Hz) or HFB (70-120 Hz) 
effect2plot       = {'rho','lag'};      % rho (corrcoefficients) or lag (xcorr lags)


for cond_i=1:length(condition2analyze)
    for band_i=1:length(band2analyze)
        % load data
        load([data_dir,filesep,'\CROSdata_',band2analyze{band_i},'_whitenoise.mat']);
        
        % dbscan parameters
        if strcmpi(condition2analyze{cond_i},'speech') && strcmpi(band2analyze{band_i},'SFB')
            mindist = 0.006; minpoints = 16;
        elseif strcmpi(condition2analyze{cond_i},'speech') && strcmpi(band2analyze{band_i},'HFB')
            mindist = 0.006; minpoints = 12;
        elseif strcmpi(condition2analyze{cond_i},'music') && strcmpi(band2analyze{band_i},'SFB')
            mindist = 0.012; minpoints = 12;
        elseif strcmpi(condition2analyze{cond_i},'music') && strcmpi(band2analyze{band_i},'HFB')
            mindist = 0.014; minpoints = 12;
        end
        
        % locate data in separate variables and plot
        rho4speech = dataMat(:,:,1);
        rho4music  = dataMat(:,:,2);
        lag4speech = LagMat(:,:,1);
        lag4music  = LagMat(:,:,2);
        
        % create arrays containing the data of interest
        if strcmpi(condition2analyze{cond_i},'speech')
            counter = 1;
            for sub_i=1:length(sub2plot)
                clear ThisSubStruct
                ThisSubStruct = cell2struct(AllChannelLabels{sub_i},names4fields,2);
                for idx=1:length(dataMat)
                    if ~isnan(rho4speech(idx,sub_i)) %% && isnan(rhos4music(idx,sub_i)) %%
                        testMat(counter,:)   = [ThisSubStruct(idx).Loc(1),ThisSubStruct(idx).Loc(2),ThisSubStruct(idx).Loc(3)];
                        RhoRange(counter,:)  = rho4speech(idx,sub_i);
                        LagRange(counter,:)  = lag4speech(idx,sub_i);
                        subIDelec(counter,:) = [sub_i,idx];
                        counter = counter + 1;
                    end
                end
            end
        elseif strcmpi(condition2analyze{cond_i},'music')
            counter = 1;
            for sub_i=1:length(sub2plot)
                clear ThisSubStruct
                ThisSubStruct = cell2struct(AllChannelLabels{sub_i},names4fields,2);
                for idx=1:length(dataMat)
                    if  ~isnan(rho4music(idx,sub_i)) %% isnan(rhos4speech(idx,sub_i)) &&
                        testMat(counter,:)   = [ThisSubStruct(idx).Loc(1),ThisSubStruct(idx).Loc(2),ThisSubStruct(idx).Loc(3)];
                        RhoRange(counter,:)  = rho4music(idx,sub_i);
                        LagRange(counter,:)  = lag4music(idx,sub_i);
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
        RhoRange(ismember(testClusters,[-1; deleteThis]))     = [];
        LagRange(ismember(testClusters,[-1; deleteThis]))     = [];
        testMat(ismember(testClusters,[-1; deleteThis]),:,:)  = [];
        subIDelec(ismember(testClusters,[-1; deleteThis]),:)  = [];

        % get anatomical labels to create anova table
        load([data_dir,filesep,'Electrode_AnatLabels_' condition2analyze{cond_i} '_' band2analyze{band_i} '.mat'])
        AnatLabels = cellstr(AnatLabels);
        disp([condition2analyze{cond_i}, '-' ,band2analyze{band_i}]);
        tbl = sortrows(tabulate(cellstr(AnatLabels)),3,'descend');
        disp(tbl);
        
        % get only the data we need
        if strcmpi(condition2analyze{cond_i},'Speech')
            these_elecs = contains(AnatLabels,'mtg') | ... 
                contains(AnatLabels,'stg') | ...
                contains(AnatLabels,'somatomotor') | ...
                contains(AnatLabels,'somatosensory') | ...
                contains(AnatLabels,'ifg') | ...
                contains(AnatLabels,'supramarginal');
        elseif strcmpi(condition2analyze{cond_i},'Music')
            these_elecs = contains(AnatLabels,'stg') | ...
                contains(AnatLabels,'mtg') | ...
                contains(AnatLabels,'pfc') | ...
                contains(AnatLabels,'supramarginal');
        end
        
        % create and save table for this effect
        if band_i==1 %&& cond_i == 1
            tbl4lmer = table(subIDelec(these_elecs,1), ...
                    repmat(convertCharsToStrings(band2analyze{band_i}),sum(these_elecs),1), ...
                    repmat(convertCharsToStrings(condition2analyze{cond_i}),sum(these_elecs),1), ...
                    AnatLabels(these_elecs)', ...
                    RhoRange(these_elecs), LagRange(these_elecs), 'VariableNames', {'subject','frequency','condition','region','r','lag'});
        else
             newtbl4lmer = table(subIDelec(these_elecs,1), ...
                    repmat(convertCharsToStrings(band2analyze{band_i}),sum(these_elecs),1), ...
                    repmat(convertCharsToStrings(condition2analyze{cond_i}),sum(these_elecs),1), ...
                    AnatLabels(these_elecs)', ...
                    RhoRange(these_elecs), LagRange(these_elecs), 'VariableNames', {'subject','frequency','condition','region','r','lag'});           
             tbl4lmer = [tbl4lmer; newtbl4lmer];
             tbl4lmer.lag = tbl4lmer.lag / 200;
             writetable(tbl4lmer,[data_dir,filesep,'AnovaNatTrackTable_wn_ANAT_' condition2analyze{cond_i} '.csv'])
        end
    end
end

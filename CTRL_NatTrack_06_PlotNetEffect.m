%%% This code corresponds to the analyses requested during peer review. It
%%% plots non-thresholded r and lag values for all electrodes. These
%%% results are presented in the suplementary materials.
%%%
%%% Osorio & Assaneo, 2025

% -------------------------------------------------------------------------
% initialize BST
cd('PATH\brainstorm3');
brainstorm
clear, clc, close
%%
% set paths
project_dir = 'PATH';
data_dir    = [project_dir,filesep,'Data'];

band2analyze      = 'HFB';      % SFB (1-8 Hz) or HFB (70-120 Hz)
condition2analyze = 'Music';    % 'speech' or 'music'
effect2plot       = 'rho';

% load data
load([data_dir,filesep,'xcorr_',band2analyze,'.mat'])
load([data_dir,filesep,'xcorr_',band2analyze,'_whitenoise_PERM.mat'])

% Data cleaning: this is to merge High Density Grids (which were preprocessed separately) 
% to ecog data of their corresponding participants
for sub_i=1:length(sub2plot)
    if contains(sub2plot{sub_i},'HD')
        r_music{sub_i-1}       = [r_music{sub_i-1}; r_music{sub_i}]; 
        r_music_perm{sub_i-1}  = [r_music_perm{sub_i-1}; r_music_perm{sub_i}];
        lag_music{sub_i-1}     = [lag_music{sub_i-1}; lag_music{sub_i}];
        r_speech{sub_i-1}      = [r_speech{sub_i-1}; r_speech{sub_i}];
        r_speech_perm{sub_i-1} = [r_speech_perm{sub_i-1}; r_speech_perm{sub_i}];
        lag_speech{sub_i-1}    = [lag_speech{sub_i-1}; lag_speech{sub_i}];
        AllChannelLabels{sub_i-1} = [AllChannelLabels{sub_i-1};  AllChannelLabels{sub_i}];
    end
end

% more cleaning
r_music(contains(sub2plot,'HD'))          = [];   r_music_perm(contains(sub2plot,'HD'))   = [];
lag_music(contains(sub2plot,'HD'))        = [];   r_speech(contains(sub2plot,'HD'))       = [];
r_speech_perm(contains(sub2plot,'HD'))    = [];   lag_speech(contains(sub2plot,'HD'))     = [];
AllChannelLabels(contains(sub2plot,'HD')) = [];   sub2plot(contains(sub2plot,'HD'))       = [];

% initialize arrays 
[pvals_speech, pvals_music, lags_speech, lags_music, fdr_speech, fdr_music] = deal(cell(1,length(sub2plot)));

% get a couple of variables we need for our loops
n_subs    = size(sub2plot,2);
n_trials  = size(r_speech{1},2);
n_conds   = 2;

% loop through subjects and electrodes to get corresponding p values
counter = 1;
for sub_i=1:n_subs
    n_elecs  = size(r_speech{sub_i},1);
    for elec_i=1:n_elecs
        CoordMat(counter,:) = AllChannelLabels{sub_i}{elec_i,5}';
        if strcmpi(condition2analyze,'speech') && strcmpi(effect2plot,'rho')
            ValRange(counter,:) = mean(r_speech{sub_i}(elec_i,:),2);
        elseif strcmpi(condition2analyze,'speech') && strcmpi(effect2plot,'lag')
            ValRange(counter,:) = mean(lag_speech{sub_i}(elec_i,:),2);
        elseif strcmpi(condition2analyze,'music') && strcmpi(effect2plot,'rho')
            ValRange(counter,:) = mean(r_music{sub_i}(elec_i,:),2);
        elseif strcmpi(condition2analyze,'music') && strcmpi(effect2plot,'lag')
            ValRange(counter,:) = mean(lag_music{sub_i}(elec_i,:),2);
        end
        counter = counter + 1;
    end
end

if strcmpi(effect2plot,'Lag')
    ValRange = ValRange/200;
end

% plot all electrodes on the common cortical template
SurfaceFile   = 'BrainstormPATH\anat\@default_subject\tess_cortex_pial_low.mat';

[hFig, iDS, iFig] = view_surface(SurfaceFile);
hFig.Color = [1 1 1];
hold on;

sh  = scatter3(CoordMat(:,1),CoordMat(:,2),CoordMat(:,3),60,ValRange,'filled');
sh.MarkerFaceAlpha = .8;
sh.SizeData = 80;
ch = colorbar;
ch.FontSize = 15;

if strcmpi(effect2plot,'rho') && strcmpi(condition2analyze,'Speech')
    set(gca,'clim',[0 0.18])
    colormap hot
    hFig.Colormap =  hFig.Colormap(6:length(hFig.Colormap)-6,:);
elseif strcmpi(effect2plot,'rho') && strcmpi(condition2analyze,'Music')
    set(gca,'clim',[0 0.16])
    colormap hot
    hFig.Colormap =  hFig.Colormap(6:length(hFig.Colormap)-6,:);
else
    set(gca,'clim',[-2 2])
    hFig.Colormap(22:64-21,:) = [[linspace(hFig.Colormap(21,1),1,11)', linspace(hFig.Colormap(21,2),1,11)', ...
        linspace(hFig.Colormap(21,3),1,11)']; [linspace(1,hFig.Colormap(64-21,1),11)', ...
        linspace(1,hFig.Colormap(64-21,2),11)', linspace(1,hFig.Colormap(64-21,3),11)']];
end

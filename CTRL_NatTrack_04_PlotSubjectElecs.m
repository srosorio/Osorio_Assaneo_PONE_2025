%%% This code plots the electrode location for each one of the subjects
%%% included in the study
%%%
%%% Osorio & Assaneo, 2025

% initialize brainstorm (for data visualization)
cd('F:\Matlab\brainstorm3');
brainstorm
%%
clear, clc, close
% -------------------------------------------------------------------------
%                             SET PARAMETERS
% -------------------------------------------------------------------------
band2analyze      = 'SFB';  % any of the two, it doesn't matter.
% -------------------------------------------------------------------------

% set paths
iEEG_dir = 'F:\Matlab\IEEG';
data_dir = [iEEG_dir,filesep,'Data'];          
    
% load data to plot
load(['F:\Matlab\IEEG\Data\CROSdata_',band2analyze,'_whitenoise.mat']);
    
% locate data in separate variables and plot
rhos4speech = dataMat(:,:,1);
rhos4music  = dataMat(:,:,2);

% get MNI cortical surface using brainstorm. We will plot data here.
SurfaceFile   = 'F:\MATLAB\brainstorm_db\iEEG\anat\@default_subject\tess_cortex_pial_low.mat'; %['C:\Users\andre\OneDrive\Documentos\MATLAB\brainstorm_db\iEEG\anat\' sub2plot '\tess_cortex_central_low.mat'];

n_subs  = length(sub2plot);
n_elecs = length(dataMat);
    
for sub_i=1:n_subs
    
    try
        [hFig, iDS, iFig] = view_surface(SurfaceFile);
    catch
        [hFig, iDS, iFig] = view_surface(SurfaceFile);
    end
    hFig.Color = [1 1 1];
    hold on;

    clear ThisSubStruct
    ThisSubStruct = cell2struct(AllChannelLabels{sub_i},names4fields,2);
    n_elecs = length({ThisSubStruct.Name});
    for elec_i=1:n_elecs
        datamat(elec_i,:)   = [ThisSubStruct(elec_i).Loc(1),ThisSubStruct(elec_i).Loc(2),ThisSubStruct(elec_i).Loc(3)];
    end
    elecspersub(sub_i) = length(datamat);
    
    % plot electrodes
    sh  = scatter3(datamat(:,1),datamat(:,2),datamat(:,3),60,'filled');
    sh.MarkerFaceAlpha = .5;
    sh.SizeData = 50;
    sh.CData    = [0    0.4470    0.7410];
    sh.MarkerEdgeColor  = [0 0 0];
    sh.MarkerEdgeAlpha = .5;
    title(sub2plot{sub_i}, 'FontSize', 18, 'FontWeight', 'Normal')
    clear datamat
    saveas(hFig,[iEEG_dir,filesep,'FiguresCROSCORR\FigureS1',filesep,sub2plot{sub_i},'.tiff'])
    close all
end


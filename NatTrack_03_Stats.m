%%% This code estimates the critical value corresponding to a predetermined 
%%% alpha and optionally performs FDR correction. If specified, it also plots 
%%% the total number of electrodes that survive statistics and saves subject-level 
%%% values.
%%%
%%% Osorio & Assaneo, 2025

clear, clc,

% set paths
project_dir = 'PATH';
data_dir    = [project_dir,filesep,'Data'];

band2analyze  = 'HFB';      % SFB (1-8 Hz) or HFB (70-120 Hz)
plot_bars     = 1;          % bar plot (1 = yes, 0 = no)
alpha         = 0.001;

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
[pvals_speech, pvals_music, lags_speech, lags_music] = deal(cell(1,length(sub2plot)));

% get a couple of variables we need for our loops
n_subs    = size(sub2plot,2);
n_trials  = size(r_speech{1},2);
n_conds   = 2;

% loop through subjects and electrodes to get corresponding p values
for sub_i=1:n_subs
    n_elecs  = size(r_speech{sub_i},1);
    n_perms  = length(r_speech_perm{sub_i});
    % obtain the mean p value across trials for each electrode within each subject
    for trial_i=1:n_trials           
        for elec_i=1:n_elecs
            % speech
            null_dist   = mean(squeeze(r_speech_perm{sub_i}(elec_i,:,:)),1);
            observed_r  = mean(r_speech{sub_i}(elec_i,:,:),2);
            pvals_speech{sub_i}(elec_i,:) = sum(null_dist >= observed_r) / n_perms;
            % music
            null_dist   = mean(squeeze(r_music_perm{sub_i}(elec_i,:,:)),1);
            observed_r  = mean(r_music{sub_i}(elec_i,:,:),2);
            pvals_music{sub_i}(elec_i,:)  = sum(null_dist >= observed_r) / n_perms;
        end
    end
 
    % finally, obtain the mean lag (ms) values across trials per electrode
    lags_speech{sub_i} = mean(mean(lag_speech{sub_i},3),2);
    lags_music{sub_i}  = mean(mean(lag_music{sub_i},3),2);
end


% number of electrodes per subject
elecspersubject = NaN(1,n_subs);
for sub_i=1:n_subs
    elecspersubject(sub_i) = numel(pvals_speech{sub_i});  
end

TotalElecs = sum(elecspersubject); % total electrodes
maxnelecs  = max(elecspersubject); % max number of electrodes across all subjects

% create a numeric array for all p values and another one for all lags 
PValMat    = NaN(maxnelecs,n_subs,n_conds); 
LagMat     = NaN(maxnelecs,n_subs,n_conds);

for sub_i=1:n_subs
    clear tmpSpeech tmpMusic
    tmpSpeech = pvals_speech{sub_i};
    tmpMusic  = pvals_music{sub_i};

    PValMat(1:length(tmpSpeech),sub_i,1) = tmpSpeech;
    LagMat(1:length(tmpSpeech),sub_i,1)  = lags_speech{sub_i};
    PValMat(1:length(tmpMusic),sub_i,2)  = tmpMusic;    
    LagMat(1:length(tmpSpeech),sub_i,2)  = lags_music{sub_i};    
end
 
% now filter rhos to keep only those whose p is significant
r_data = cell(n_subs,n_conds);

for sub_i=1:n_subs
    % subject-specific data
    r_data_speech = NaN(maxnelecs,n_trials);
    r_data_music  = NaN(maxnelecs,n_trials);   
    for trial_i=1:n_trials
        n_elecs = size(r_speech{sub_i},1);
        % whatever is not significant becomes NaN
        for pval_i=1:n_elecs
            if pvals_speech{sub_i}(pval_i) < alpha
                r_data_speech(pval_i,trial_i) = mean(r_speech{sub_i}(pval_i,trial_i,:),3);
            end
            if pvals_music{sub_i}(pval_i) < alpha
                r_data_music(pval_i,trial_i) = mean(r_music{sub_i}(pval_i,trial_i,:),3);
            end
        end
    end    
    % get the mean rho coefficient across trials
    r_data{sub_i,1} = nanmean(r_data_speech,2);
    r_data{sub_i,2} = nanmean(r_data_music,2);
end

% For simplicity, turn all the data above (pvalues, rho values and lag values)
% into arrays where dim 1 = electrodes, dim 2 = subjects and dim 3 = conditions).
% For conditions, (:,:,1) = Speech and (:,:,2) = music
dataMat   = NaN(maxnelecs,n_subs,n_conds);

for sub_i=1:size(r_data,1)
    clear tmpSpeech tmpMusic
    tmpSpeech = r_data{sub_i,1};
    tmpMusic  = r_data{sub_i,2};

    dataMat(1:length(tmpSpeech),sub_i,1) = tmpSpeech;
    LagMat(isnan(r_data{sub_i,1}),sub_i,1) = NaN;
    PValMat(isnan(r_data{sub_i,1}),sub_i,1) = NaN;
    dataMat(1:length(tmpMusic),sub_i,2)  = tmpMusic;    
    LagMat(isnan(r_data{sub_i,2}),sub_i,2) = NaN;
    PValMat(isnan(r_data{sub_i,2}),sub_i,2) = NaN;
end

% How many subjects show a statistically significant effect per condition?
subeffect   = sum(any(dataMat(:,:,1)));
disp(['SPEECH: Statistically significant data in ' num2str(subeffect) ' out of ' num2str(length(sub2plot)) ' subjects'])
subeffect   = sum(any(dataMat(:,:,2)));
disp(['MUSIC: Statistically significant data in ' num2str(subeffect) ' out of ' num2str(length(sub2plot)) ' subjects'])

% save data
save([data_dir,filesep,'CROSdata_' band2analyze '_whitenoise.mat'],  ...
    'dataMat','LagMat','PValMat','sub2plot','TotalElecs','AllChannelLabels','names4fields');

% bar plots for total number of electrodes per condition
if plot_bars == 1

    ElecsperSub_Speech = sum(~isnan(dataMat(:,:,1)));
    ElecsperSub_Music  = sum(~isnan(dataMat(:,:,2)));
    ElecsperSub_Both   = sum(~isnan(dataMat(:,:,1)) & ~isnan(dataMat(:,:,2)));
    
    % save data
    save([data_dir,filesep,'ElecsPerSub_',band2analyze],'ElecsperSub_Speech','ElecsperSub_Music','ElecsperSub_Both')
    
    % create bar plot
    colors = [0.7176 0.2745 1.0000; ...
              1.0000 0.4118 0.1608; ...
              0.4667 0.6745 0.1882];    
    ElecsPerSub = [ElecsperSub_Music',ElecsperSub_Speech',ElecsperSub_Both'];
    
    % let's create a table with basic statistics
    RhosMusic  = dataMat(:,:,2);
    RhosMusic  = RhosMusic(~isnan(RhosMusic));
    RhosSpeech = dataMat(:,:,1);
    RhosSpeech = RhosSpeech(~isnan(RhosSpeech));
    pMusic     = PValMat(:,:,2);
    pMusic     = pMusic(~isnan(dataMat(:,:,2)));
    pSpeech    = PValMat(:,:,1);
    pSpeech    = pSpeech(~isnan(dataMat(:,:,1)));
    RhosBoth   = dataMat(~isnan(dataMat));
    pBoth      = PValMat(~isnan(dataMat));
    
    disp(['N music = ' num2str(length(RhosMusic))]);
    disp(['N speech = ' num2str(length(RhosSpeech))]);

    StatsTable = table([mean(RhosMusic); std(RhosMusic); mean(pMusic)], ...
                       [mean(RhosSpeech);std(RhosSpeech);mean(pSpeech)], ...
                       [mean(RhosBoth); std(RhosBoth); mean(pBoth)], ...
                       'VariableNames',{'Music','Speech','Both'}, ...
                       'RowNames',{'mean rho','SD rho','mean p'});
    display(StatsTable);
    
    figure(1), clf
    bh1 = bar(sum(ElecsPerSub),'FaceColor','flat');
    bh1.CData = colors;
    bh1.EdgeColor = [1 1 1];
    xticklabels({'Music','Speech','Both'});
    xtickangle(45);
    ylabel('Electrode count');
    set(gca,'FontSize',20,'FontName','Arial');
    title(band2analyze,'FontWeight','normal');
    box off
    
    % save number of electrodes per subject as a table
    ElecsPersSub = array2table(ElecsPerSub,'RowNames',sub2plot,'VariableNames',{'Music','Speech','Both'});
    writetable(ElecsPersSub,[data_dir,filesep,'ElecsPerSub_',band2analyze,'.xlsx'])
end
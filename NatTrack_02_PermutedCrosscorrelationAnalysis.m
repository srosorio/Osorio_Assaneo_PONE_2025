%%% This code estimates a null distribution for the crosscorrelation function 
%%% between electrophysiological signals and naturalistic acoustic (speech
%%% and music signals). 
%%%
%%% This code is computationally demanding and was therefor run a computing
%%% cluster, separately for each subject
%%%
%%%% Osorio & Assaneo, 2025

clear, clc,
warning('off','all') % turn down warnings
rng('shuffle')       % initialize this matlab session with a random seed 
s = rng;             % keep track of the seeds used for each subject

% set path to relevant directories
project_dir = 'PATH';
cd(project_dir)
data_dir = [project_dir,filesep,'Matlab_data'];

% frequency band to analyze (run each band separately)
band2analyze = 'SFB';

% list of subjects
subjects     = {'sub-02','sub-03','sub-05','sub-06','sub-10','sub-12','sub-16','sub-18', ...
                'sub-19','sub-20','sub-22','sub-24','sub-25','sub-26','sub-27','sub-34', ...
                'sub-36','sub-36HD','sub-39','sub-40','sub-45','sub-45HD','sub-46','sub-48', ...
                'sub-51','sub-54','sub-55','sub-58','sub-59','sub-60','sub-61','sub-63'};
            
% subject for which permutations will be run
sub2permute  = 'sub-02';
sub2perm_idx = strcmpi(sub2permute,subjects);

% load neural data
if strcmpi(band2analyze,'SFB')
    load([data_dir,filesep,'fieldtrip_structures_SFB']);
elseif strcmpi(band2analyze,'HFB')
    load([data_dir,filesep,'fieldtrip_structures_HFB']);
end

% load acoustic envelopes
load([data_dir,filesep,'envelopes_music.mat']);
load([data_dir,filesep,'envelopes_speech.mat']);

% all this must be the same as in NatTrack_02_ObservedCrosscorrelationAnalysis.m
fs            = 250;    % sampling rate
n_trials      = length(AllDataStructuresFT{1,1}.trial);
n_conditions  = size(AllDataStructuresFT,2);
maxlag        = 400;    % maximum pos and neg lag for crosscorrelation
nperms        = 1000;   % 1000 for whitenoise

%initialize cell arays for the data we need
[r_music_perm,lag_music_perm,r_speech_perm,lag_speech_perm] = deal(cell(1,1));

% get only the subject of interest
this_subject = {AllDataStructuresFT{sub2perm_idx,1},AllDataStructuresFT{sub2perm_idx,2}};
n_electrodes = size(this_subject{1,1}.trial{1},1); % number of electrodes (subject-dependent)

% Print message to keep track of subjects and make sure everything is okay
disp(['NatTrack >>> Estimating null distribution for crosscorrelation coefficients for subject ' sub2permute]);

% trim signal length if not already done
for cond_i=1:n_conditions
    for trial_i=1:n_trials
        this_subject{1,cond_i}.trial{1,trial_i} = this_subject{1,cond_i}.trial{1,trial_i}(:,1:length(envelope_speech));
        this_subject{1,cond_i}.time{trial_i}    = this_subject{1,cond_i}.time{trial_i}(:,1:length(envelope_speech));
    end
end

for perm_i=1:nperms
    
    % get the time vector
    time_vector  = this_subject{1,cond_i}.time{trial_i};
    
    % create 6 white noise segments of 30 seconds per condition
    wn_musicsignal  = resample(wgn(6,length(envelope_speech),1)',time_vector,16000)';
    wn_speechsignal = resample(wgn(6,length(envelope_speech),1)',time_vector,16000)';
    
    % now we apply the cochlear filter to the wn signals as we did the real signals
    for sig_i=1:size(wn_speechsignal,1)
        
        % v corresponds to the signals in the 128 frequencies between 180 and 7142 Hz
        [music_v,~]  = wav2aud2(wn_musicsignal(sig_i,:),[5 8 -2 0]);
        [speech_v,~] = wav2aud2(wn_speechsignal(sig_i,:),[5 8 -2 0]);
        
        % detrend the signals
        music_v  = detrend(music_v);
        speech_v = detrend(speech_v);
        
        % get the envelope by averaging across frequencies and replace the original envelopes with the dummy envelopes
        envelope_music(sig_i,:)  = mean(music_v');
        envelope_speech(sig_i,:) = mean(speech_v');
    end
    
    % do this for each trial and electrode
    for trial_i=1:n_trials
        for elec_i=1:n_electrodes
            % speech
            acoustic_signal  = envelope_speech(trial_i,:);
            brain_signal     = this_subject{1,1}.trial{trial_i}(elec_i,:);
            [tempr,templags] = xcorr(zscore(brain_signal), ...
                zscore(acoustic_signal),maxlag,'normalized');
            r_speech_perm{1}(elec_i,trial_i,perm_i) = max(tempr);
            lag_speech_perm{1}(elec_i,trial_i,perm_i)= templags(tempr == max(tempr));
            % music
            acoustic_signal  = envelope_music(trial_i,:);
            brain_signal     = this_subject{1,2}.trial{trial_i}(elec_i,:);
            [tempr,templags] = xcorr(zscore(brain_signal), ...
                zscore(acoustic_signal),maxlag,'normalized');
            r_music_perm{1}(elec_i,trial_i,perm_i) = max(tempr);
            lag_music_perm{1}(elec_i,trial_i,perm_i)= templags(tempr == max(tempr));
        end
    end
end

% save data
disp(['NatTrack >>> saving xcorr_' band2analyze '_whitenoise_PERM.mat']);
save([data_dir,filesep,subjects{sub2perm_idx} '_xcorr_' band2analyze '_whitenoise_PERM.mat'], ...
    'r_speech_perm','lag_speech_perm','r_music_perm', ...
    'lag_music_perm','band2analyze','subjects','AllChannelLabels','names4fields','s');

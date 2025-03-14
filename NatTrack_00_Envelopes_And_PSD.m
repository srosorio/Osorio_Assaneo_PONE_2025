%%% This code estimates the cochlear envelopes and plots the PSD of 7 the 
%%% music and speech signals
%%%
%%% Osorio & Assaneo, 2025

clear, clc
cd('PATH\Scripts');

% add NSL toolbox to path
addpath('PATH\NSL');

% create a directory to save data
project_dir = 'PATH';
data_dir = [project_dir,filesep,'Data'];
Fig_Dir  = [project_dir,filesep,'FiguresCROSCORR'];

% load audio file
[y,fs] = audioread('PATH\audio.wav');
dt     = 1/fs;                   
tfull  = 0:dt:(length(y)*dt)-dt;  % time vector

% Cut full audio signal into speech and music segments
segstart  = 1;
seglength = 30;
totalsegments = floor((length(tfull)/48000) / 30);

% initiallize variables
speech_segments = cell(1,totalsegments);
music_segments  = cell(1,totalsegments);

for idx=1:totalsegments
    if mod(idx,2)==0
        % speech
        speech_segments{idx} = y(dsearchn(tfull',segstart):dsearchn(tfull',seglength*idx),:);
        segstart = segstart + seglength;
    elseif mod(idx,2)~=0
        % music 
        music_segments{idx}  = y(dsearchn(tfull',segstart):dsearchn(tfull',seglength*idx),:);
        segstart = segstart + seglength;
    end
end

% get rid of empty cells
music_segments(cellfun(@isempty,music_segments))   = [];
speech_segments(cellfun(@isempty,speech_segments)) = [];

% create a new time vector for segmented data
tseg = 0:dt:length(music_segments{1,1})*dt-dt; 

%% Get cochlear envelopes
for idx=1:length(music_segments)-1
    % resample signals
    musicsignal   = resample(music_segments{idx},tseg,16000);  %%%% resamples and then pads to 12 secs in total
    speechsignal  = resample(speech_segments{idx},tseg,16000);
    
    % v corresponds to the signals in the 128 frequencies between 180 and 7142 Hz
    [music_v_all(:,:,idx),freqs]      = wav2aud2(musicsignal,[5 8 -2 0]);
    [speech_v_all(:,:,idx),~]     = wav2aud2(speechsignal,[5 8 -2 0]);
    
    music_v  = detrend(music_v_all(:,:,idx));
    speech_v = detrend(speech_v_all(:,:,idx));
    
    % get the envelope by averaging across frequencies
    envelope_music(idx,:)  = mean(music_v');
    envelope_speech(idx,:) = mean(speech_v');
end

% save envelopes
save([data_dir,filesep,'envelopes_music.mat'],'envelope_music');
save([data_dir,filesep,'envelopes_speech.mat'],'envelope_speech');

%% Get PSDs
window   = 1200;
noverlap = 600;

% get powers spectrum density for music
[psd_music, frex_music] = pwelch(envelope_music',window,noverlap,[],200);
% normalize spectrum
psd_music = psd_music ./ sum(psd_music);

% plot
figure(1), clf
h1 = plot(frex_music,psd_music(:,1:6),'linew',1);
xlim([0 8]); 
ylim([0 .3]); yticks(0:0.05:0.3);
xlabel('Frequency (Hz)'); ylabel('Normalized power');
title('Music','FontWeight','normal');
set(gca,'FontSize',12);
box off, hold on

for idx=1:numel(h1)-1
    h1(idx).Color = [h1(idx).Color 0.3];
end

h2 = plot(frex_music,mean(psd_music(:,1:6),2),'linew',2); 
legend({'Segment 1', 'Segment 2', 'Segment 3', 'Segment 4', 'Segment 5', 'Segment 6','Average'},'FontSize',14,'Box','off');
saveas(gca,[Fig_Dir,filesep,'MusicStimuliPSD.tiff']);

% now get the PSD for speech and normalize spectrum
[psd_speech, frex_speech] = pwelch(envelope_speech',window,noverlap,[],200);
psd_speech = psd_speech ./ sum(psd_speech);

% plot
figure(2), clf
h1 = plot(frex_speech,psd_speech,'linew',1);
xlim([0 8]); 
ylim([0 .1]); yticks(0:0.02:0.1);
xlabel('Frequency (Hz)'); ylabel('Normalized power');
title('Speech','FontWeight','normal');
set(gca,'FontSize',12);
box off, hold on

for idx=1:numel(h1)
    h1(idx).Color = [h1(idx).Color 0.3];
end

plot(frex_speech,mean(psd_speech,2),'linew',2,'Color', h2.Color);
legend({'Segment 1', 'Segment 2', 'Segment 3', 'Segment 4', 'Segment 5', 'Segment 6','Average'},'FontSize',14,'Box','off');
saveas(gca,[Fig_Dir,filesep,'SpeechStimuliPSD.tiff']);
%% Plot raw acoustic signals
close all
colors  = [0.7176 0.2745 1.0000; ...
           0.9412 0.5804 0.3373];   
            
% plot music and speech signals
for idx=1:length(music_segments)-1
    figure,clf
    subplot(3,3,1:3)
    ph   = plot(tseg,music_segments{idx}(:,1),'color',colors(1,:));
    set(gca,'ylim',[-.8 .8],'xlim',[0.9 28.90],'box','off','XTick',[],'FontSize',12)
    ylabel('Norm. amplitude')
    subplot(3,3,4:6)
    imagesc(linspace(1,30,length(music_v)),freqs,music_v_all(:,:,idx)')
    set(gca,'YDir','normal')
    set(gca,'clim',[0 1.5], 'xlim', [0.9 28.90],'box','off','XTick',[],'FontSize',12)
    ylabel('Frequency (Hz)')
    saveas(gca,[Fig_Dir,filesep,'Acoustics_music_segment' num2str(idx) '.tiff'])
    subplot(3,3,7:9)
    envh = plot(linspace(0,30,length(envelope_music)),envelope_music(idx,:),'k');
    set(gca,'ylim',[-.8 .8],'xlim',[0.9 28.90],'box','off','FontSize',12)
    ylabel('Norm. amplitude'); xlabel('Time (s)')
    saveas(gca,[Fig_Dir,filesep,'Acoustics_music_segment' num2str(idx), '.tiff'])
end

for idx=1:length(speech_segments)
    figure,clf
    subplot(3,3,1:3)
    ph   = plot(tseg,speech_segments{idx}(:,1),'color',colors(2,:));
    set(gca,'ylim',[-.8 .8],'xlim',[0.9 28.90],'box','off','XTick',[],'FontSize',12)
    ylabel('Norm. amplitude')
    subplot(3,3,4:6)
    imagesc(linspace(1,30,length(speech_v)),freqs,speech_v_all(:,:,idx)')
    set(gca,'YDir','normal')
    set(gca,'clim',[0 1.5], 'xlim',[0.9 28.90],'XTick',[],'FontSize',12)
    ylabel('Frequency (Hz)')
    subplot(3,3,7:9)
    envh = plot(linspace(0,30,length(envelope_speech)),envelope_speech(idx,:),'k');
    set(gca,'ylim',[-.8 .8],'xlim',[0.9 28.90],'box','off','FontSize',12)
    ylabel('Norm. amplitude'); xlabel('Time (s)')  
    saveas(gca,[Fig_Dir,filesep,'Acoustics_speech_segment' num2str(idx), '.tiff'])
end

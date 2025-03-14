%%% This code plots the statistical effects for mean r and lag, as per the
%%% results of the LMER analyses (joint model)
%%%
%%% Osorio & Assaneo, 2025

% set paths
clear,clc
project_dir = 'PATH';
data_dir    = [project_dir,filesep,'Data'];

tbl   = readtable([data_dir,filesep,'AnovaNatTrackTable_wn_ANAT.csv']);
means = [mean(tbl.r(contains(tbl.condition,'Music'))),mean(tbl.r(contains(tbl.condition,'Speech')))];
err   = [std(tbl.r(contains(tbl.condition,'Music')))/sqrt(length(tbl.r(contains(tbl.condition,'Music')))),std(tbl.r(contains(tbl.condition,'Speech')))/sqrt(length(tbl.r(contains(tbl.condition,'Speech'))))];

% plot effect of condition on r coefficients
subplot(1,3,2)
bh = bar(means, 'FaceColor','Flat','BarWidth', .5); 
hold on;
bh.EdgeColor = [1 1 1];
bh.CData(1,:) = [0.7176 0.2745 1.0000];
bh.CData(2,:) = [1.0000 0.4118 0.1608];
bh.FaceAlpha = .8;
ylim([0 .16]); xlim([0 3]); yticks(0:0.04:0.16)
xtickangle(45)
errorbar(1,means(1),err(1),'k')
errorbar(2,means(2),err(2),'k')
xticklabels({'Music','Speech'})
ylabel('Mean cross-correlation (r)')
set(gca,'FontSize',12)
box off
axis square
%

means = [[mean(tbl.lag(contains(tbl.frequency,'SFB') & contains(tbl.condition,'Music'))) ; ...
         mean(tbl.lag(contains(tbl.frequency,'HFB') & contains(tbl.condition,'Music')))], ...
         [mean(tbl.lag(contains(tbl.frequency,'SFB') & contains(tbl.condition,'Speech'))); ...
         mean(tbl.lag(contains(tbl.frequency,'HFB') & contains(tbl.condition,'Speech')))]];
     
err = [[std(tbl.lag(contains(tbl.frequency,'SFB') & contains(tbl.condition,'Music'))) / sqrt(length(tbl.lag(contains(tbl.frequency,'SFB') & contains(tbl.condition,'Music')))); ...
         std(tbl.lag(contains(tbl.frequency,'HFB') & contains(tbl.condition,'Music'))) / sqrt(length(tbl.lag(contains(tbl.frequency,'HFB') & contains(tbl.condition,'Music'))))], ...
         [std(tbl.lag(contains(tbl.frequency,'SFB') & contains(tbl.condition,'Speech'))) / sqrt(length(tbl.lag(contains(tbl.frequency,'SFB') & contains(tbl.condition,'Speech')))); ...
         std(tbl.lag(contains(tbl.frequency,'HFB') & contains(tbl.condition,'Speech'))) / sqrt(length(tbl.lag(contains(tbl.frequency,'HFB') & contains(tbl.condition,'Speech'))))]];

% plot effect of condition on mean lag values
subplot(1,3,3)
bh = bar(means, 'FaceColor','Flat','BarWidth', .5); 
hold on;
bh(1).EdgeColor = [1 1 1];
bh(2).EdgeColor = [1 1 1];
bh(1).CData = [0.7176 0.2745 1.0000; 1.0000 0.4118 0.1608];
bh(2).CData = [0.8157 0.6000 0.9490; 0.9804 0.5804 0.4118];
bh(1).FaceAlpha = .8;
bh(2).FaceAlpha = .8;
ylim([-.2 1]); 
xtickangle(45)
plot(get(gca,'xlim'),[0 0],'k');

for k = 1:numel(bh)                                                      % Earlier MATLAB Versions
    ctr(k,:) = bsxfun(@plus, bh(k).XData, [bh(k).XOffset]');
    ydt(k,:) = bh(k).YData;
end

for idx=1:length(ctr)
    errorbar(ctr(1,idx), ydt(1,idx), err(1,idx), 'k', 'MarkerSize',0.1)
    errorbar(ctr(2,idx), ydt(2,idx), err(2,idx), 'k', 'MarkerSize',0.1)
end
xticklabels({'Music','Speech'})
ylabel('Mean lag (s)')
set(gca,'FontSize',12)
box off
axis square
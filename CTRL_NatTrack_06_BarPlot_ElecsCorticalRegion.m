%%% This code creates bar plots showing the number of significant
%%% electrodes per anatomical regions
%%%
%%% Osorio & Assaneo, 2025

clear, clc
project_dir = 'PATH';
data_dir    = [project_dir,filesep,'Data'];

%%% ------- bar plot for music --------
% SFB
load([data_dir,filesep,'Electrode_AnatLabels_Music_SFB.mat'])
Music_SFB = sortrows(tabulate(cellstr(AnatLabels)),1,'ascend');
Music_SFB(2,:) = [];
% HFB
load([data_dir,filesep,'Electrode_AnatLabels_Music_HFB.mat'])
Music_HFB = sortrows(tabulate(cellstr(AnatLabels)),1,'ascend');
Music_HFB(4,:) = [];

Values = [[cell2mat(Music_SFB(1,2));0;cell2mat(Music_SFB(2:3,2));0; ... 
            cell2mat(Music_SFB(4:end,2))], cell2mat(Music_HFB(:,2))];
        
subplot(1,3,2)
bh = bar(Values,'BarWidth',1,'FaceColor','flat');
bh(1).CData = repmat([0.7176 0.2745 1.0000],9,1);
bh(2).CData = repmat([0.8157 0.6000 0.9490],9,1);
bh(1).FaceAlpha = .8; bh(2).FaceAlpha = .8; 
bh(1).EdgeColor = [1 1 1]; bh(2).EdgeColor = [1 1 1];
xticklabels({'IFG','IP','MTG','PFC','Premotor','Precentral','Postcentral','STG','SG'})
xtickangle(45)
legend({'SFB','HFB'}, 'box', 'off', 'Orientation', 'horizontal', 'Location', 'north')
ylabel('Electrode count')
title('Music','FontWeight','normal')
box off
set(gca,'FontSize',12, 'ylim', [0 50]);
axis square
tmptbl = table({'IFG','IP','MTG','PFC','Premotor','Precentral','Postcentral','STG','SG'}', ...
    Values(:,1), Values(:,2), 'VariableNames',{'Region','SFB','HFB'});
% writetable(tmptbl,[data_dir,filesep,'ElecCountPerAnatRegion_Music.xlsx']);

%%% ------- bar plot for speech --------
% SFB
load([data_dir,filesep,'Electrode_AnatLabels_Speech_SFB.mat'])
Speech_SFB = sortrows(tabulate(cellstr(AnatLabels)),1,'ascend');
% HFB
load([data_dir,filesep,'Electrode_AnatLabels_Speech_HFB.mat'])
Speech_HFB = sortrows(tabulate(cellstr(AnatLabels)),1,'ascend');

Values = [cell2mat(Speech_SFB(:,2)), [cell2mat(Speech_HFB(1:3,2));0; ... 
    cell2mat(Speech_HFB(4:end,2))]];

subplot(1,3,3)
bh = bar(Values,'BarWidth',1,'FaceColor','flat');
bh(1).CData = repmat([1.0000 0.4118 0.1608],8,1);
bh(2).CData = repmat([0.9804 0.5804 0.4118],8,1);
bh(1).FaceAlpha = .8; bh(2).FaceAlpha = .8; 
bh(1).EdgeColor = [1 1 1]; bh(2).EdgeColor = [1 1 1];
xticklabels({'IFG','MTG','PFC','Premotor','Precentral','Postcentral','STG','SG'})
xtickangle(45)
legend({'SFB','HFB'}, 'box', 'off', 'Orientation', 'horizontal', 'Location', 'north')
ylabel('Electrode count')
title('Speech','FontWeight','normal')
box off
set(gca,'FontSize',12, 'ylim', [0 250], 'xlim',[0 9])
axis square

tmptbl = table({'IFG','MTG','PFC','Premotor','Precentral','Postcentral','STG','SG'}' ...
    ,Values(:,1), Values(:,2), 'VariableNames',{'Region','SFB','HFB'});
% writetable(tmptbl,[data_dir,filesep,'ElecCountPerAnatRegion_Speech.xlsx']);
%%
%%% ------- This is for both music and speech --------
clf
subplot(1,3,1);
Music_SFB  = Music_SFB(contains(Music_SFB(:,1),{'stg','mtg','supramarginal'}),:);
Music_HFB  = Music_HFB(contains(Music_HFB(1),{'stg','mtg','supramarginal'}),:);
Speech_SFB = Speech_SFB(contains(Speech_SFB(:,1),{'stg','mtg','supramarginal'}),:);
Speech_HFB = Speech_HFB(contains(Speech_HFB(:,1),{'stg','mtg','supramarginal'}),:);

Values = [cell2mat(Music_SFB(:,2)), cell2mat(Music_HFB(:,2)); cell2mat(Speech_SFB(:,2)), ...
    cell2mat(Speech_HFB(:,2)) ];

bh = bar(Values,'BarWidth',1,'FaceColor','flat');
xticklabels({'STG','MTG','SG','STG','MTG','SG'})
xtickangle(45)
set(gca,'FontSize',12, 'ylim', [0 250], 'xlim',[0 7])
bh(1).CData = [repmat([0.7176 0.2745 1.0000],3,1); repmat([1.0000 0.4118 0.1608],3,1)];
bh(2).CData = [repmat([0.8157 0.6000 0.9490],3,1); repmat([ 0.9804 0.5804 0.4118],3,1)];
bh(1).FaceAlpha = .8; bh(2).FaceAlpha = .8; 
bh(1).EdgeColor = [1 1 1]; bh(2).EdgeColor = [1 1 1];
ylabel('Electrode count');
box off
axis square

tmptbl = table({'STG','MTG','SG','STG','MTG','SG'}', ...
    Values(:,1), Values(:,2), 'VariableNames',{'Region','SFB','HFB'});
% writetable(tmptbl,[data_dir,filesep,'ElecCountPerAnatRegion_Both.xlsx']);

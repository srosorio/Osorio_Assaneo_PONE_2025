%%% This code plots the number of significant electrodes after statistical 
%%% thresholding and DBSCAN analyses.
%%%
%%% Osorio & Assaneo, 2025

condition2analyze = {'music','speech','both'};
band2analyze      = {'HFB','SFB'};
project_dir       = 'PATH';
data_dir          = [project_dir,filesep,'Data'];

colors = [0.7176 0.2745 1.0000; ...
    1.0000 0.4118 0.1608; ...
    0.4667 0.6745 0.1882];

for band_i=1:length(band2analyze)
    for cond_i = 1:length(condition2analyze)
        load([data_dir,filesep,'dbscan_results_' condition2analyze{cond_i} '_' band2analyze{band_i} '.mat'])
        elec_count(cond_i) = length(datamat);
    end
    figure(band_i),
    bh = bar(elec_count, 'FaceColor','Flat');
    bh.EdgeColor = [1 1 1];
    bh.CData = colors;
    ylim([0 1000]);
    xticklabels({'Music','Speech','Both'});
    ylabel('Electrode count');
    set(gca,'FontSize',20,'FontName','Arial');
    xtickangle(45);
    box off
end
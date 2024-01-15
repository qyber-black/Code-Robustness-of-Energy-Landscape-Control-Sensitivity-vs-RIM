% SPDX-FileCopyrightText: Copyright (C) 2023 Sean Patrick O'Neil <seanonei@usc.edu>
% SPDX-FileCopyrightText: Copyright (C) 2024 S M Shermer <lw1660@gmail.com>
% SPDX-License-Identifier: CC-BY-SA-4.0 
%
function [] = plot_all_heatmaps_rings()

% This routine takes chain controllers from results/rings
% and the corresponding kde results from results/kde- and produces
% heatmaps of the individual controllers saved in figures/heatmaps-rings

if ~exist('../figures/heatmaps-rings/','dir')
    mkdir('../figures/heatmaps-rings/');
end

num  = 100;
in   = 1;
type = {'fidelity','dephasing','overlap'};
index = 1:10;

for N = 5:6

    if N == 5
        target = [2 3];
    else
        target = [2 3 4];
    end

    for x = 1:length(target)

	    out = target(x);
	    
        for q = 1:3
			
            % Load data 
		    opt = type{q};
    	
			loadtag1 = sprintf('../results/controllers-rings/%s_%d-ring_1-%d.mat',opt,N,out);
 			load(loadtag1);
			loadtag2 = sprintf('../results/kde-rings/log_kde_%s_%d-ring_1-%d.mat',opt,N,out);
 			load(loadtag2);

            for k=1:length(index)

                if ~isfield(density{k},'min_err')
                    Density{k}.min_err = min(min(Error{k}));
                end
                if ~isfield(density{k},'max_err')
                    Density{k}.max_err = max(max(Error{k}));
                end
                PlotHeatMaps(density{k},sys{k})

                set(gcf,'PaperSize',[8 6],'PaperPosition',[0 0 8 6]);
		        set(findall(gcf,'-property','FontSize'),'FontSize',12);

    		    figtag = sprintf('../figures/heatmaps-rings/%s_%d-ring_1-%d_ctrl-%d.fig',opt,N,out,k);
	    	    savefig(figtag);
		        saveas(gcf,sprintf('%s.png',figtag(1:end-4)),'png');
                close all
            end

		end
	end
end

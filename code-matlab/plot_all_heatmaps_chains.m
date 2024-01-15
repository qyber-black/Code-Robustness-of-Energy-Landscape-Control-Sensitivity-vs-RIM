% SPDX-FileCopyrightText: Copyright (C) 2023 Sean Patrick O'Neil <seanonei@usc.edu>
% SPDX-FileCopyrightText: Copyright (C) 2024 S M Shermer <lw1660@gmail.com>
% SPDX-License-Identifier: CC-BY-SA-4.0 

function [] = plot_all_heatmaps_chains(noise_level)

% This routine takes chain controllers from results/controllers-chains
% and the corresponding kde results from results/kde-chains and produces
% heatmaps of the individual controllers saved in figures/heatmaps-chains
%
if ~exist('../figures/heatmaps-chains/','dir')
    mkdir('../figures/heatmaps-chains/');
end

num   = 100;
in    = 1;
type  = {'nmplus','ppo','snob'};
index = 1:10; % must be less length(density)=100 but plotting every controller not recommended!
for N = 5:6

    if N == 5
        target = [3 5];
    else
        target = [4 6];
    end

    for x = 1:2

	    out = target(x);
	    
        for q = 1:3
			
            % Load data 
		    opt = type{q};

			loadtag1 = sprintf('../results/controllers-chains/noise_%s_%s_%d-chain_1-%d.mat',noise_level,opt,N,out);
 			load(loadtag1);
			loadtag2 = sprintf('../results/kde-chains/log_kde_%s_%s_%d-chain_1-%d.mat',noise_level,opt,N,out);
 			load(loadtag2);

            for k=index 

                if ~isfield(density{k},'min_err')
                    density{k}.min_err = min(min(Error{k}));
                end
                if ~isfield(density{k},'max_err')
                    density{k}.max_err = max(max(Error{k}));
                end
                PlotHeatMaps(density{k},sys{k})

                set(gcf,'PaperSize',[8 6],'PaperPosition',[0 0 8 6]);
		        set(findall(gcf,'-property','FontSize'),'FontSize',10);

    		    figtag = sprintf('../figures/heatmaps-chains/noise_%s_%s_%d-chain_1-%d_ctrl-%d.fig',noise_level,opt,N,out,k);
	    	    savefig(figtag);
		        saveas(gcf,sprintf('%s.png',figtag(1:end-4)),'png');
                close all
            end

		end
	end
end

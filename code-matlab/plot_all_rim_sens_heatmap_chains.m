%% Plot RIM heat maps for chains (rim_sens_heat_map_chain.m)
%
function [] = plot_all_rim_sens_heatmap_chains(noise_)
%
% SPDX-FileCopyrightText: Copyright (C) 2023 Sean Patrick O'Neil <seanonei@usc.edu>
% SPDX-FileCopyrightText: Copyright (C) 2023 S Shermer <lw1660@gmail.com>, Swansea University
% SPDX-License-Identifier: CC-BY-SA-4.0 

% This function produces heatmaps of the RIM versus dephasing strength for
% chains

if ~exist('../figures/heat_maps/','dir')
    mkdir('../figures/heat_maps/');
end

in = 1;
noise_level = '0.1';

type = {'nmplus','ppo','snob'};
index = 1:100;

for N = 5:6

    if N == 5
        target = [3 5];
    else
        target = [4 6];
    end

    for x = 1:2
	    out = target(x);

		for q = 1:3
    		opt = type{q};

			load(sprintf('../results/log_sens-chains/log_sens_%s_%s_%d-chain_1-%d.mat',noise_level,opt,N,out));
			load(sprintf('../results/rim-chains/rim_noise_%s_%s_%d-chain_1-%d',noise_level,opt,N,out));
			rim_data=rim_data';
			err_rim=rim_data(:,1);
			sens = abs(sens(:,1001));
			arim = rim_data - err_rim;
			index = 1:100;
			delta = 0:0.001:0.1;
			Q = [sens arim];
			Z = sortrows(Q);
			sens = Q(:,1);
			arim = Q(:,2:102);
			pcolor(sens,delta,arim'), shading flat, colorbar;
			set(gca,'XScale','linear','YScale','linear')
			xlabel('differential sensitivity at \delta=0');
			ylabel('perturbation strength \delta');
			hCB = colorbar;
			hCB.Label.String='Adjusted RIM_1(\delta)'
			set(gcf,'PaperSize',[6 4],'PaperPosition',[0 0 6 4]);
			set(findall(gcf,'-property','FontSize'),'FontSize',12);
			figtag = sprintf('../figures/heat_maps/heat_map_chain_%s_%d-%d',opt,N,out);
			saveas(gcf,figtag,'png');
	    end
    end
end

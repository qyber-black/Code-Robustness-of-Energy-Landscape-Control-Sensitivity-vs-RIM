%% Plot RIM heat maps for rings (rim_sens_heat_map_ring.m)

function [] = plot_all_rim_sens_heatmap_rings()

% SPDX-FileCopyrightText: Copyright (C) 2023 Sean Patrick O'Neil <seanonei@usc.edu>
% SPDX-FileCopyrightText: Copyright (C) 2023 S Shermer <lw1660@gmail.com>, Swansea University
% SPDX-License-Identifier: CC-BY-SA-4.0 

% This function produces heatmaps of the RIM versus dephasing strength for
% rings

if ~exist('../figures/heat_maps/','dir')
    mkdir('../figures/heat_maps/');
end

option = {'dephasing';'fidelity';'overlap'};

for x = 1:3
    opt = option{x}
	for N = 5:6
		for out = 2:floor(N/2)+1

		load(sprintf('../results/log_sens-rings/log_sens_%s_%d-ring_1-%d.mat',opt,N,out));
		load(sprintf('../results/rim-rings/rim_%s_%d-ring_1-%d',opt,N,out));
		clear sens

		for k = 1:100
    		sens(k,1) = log_sens(k,1001)*err(k);
		end
		rim_data=rim_data';
		err_rim=rim_data(:,1);
		arim = rim_data - err_rim;
		delta = 0:0.001:0.1;
		Q = [sens arim];
		Z = sortrows(Q);
		sens = Z(:,1);
		arim = Z(:,2:102);
		pcolor(sens,delta,arim'), shading flat, colorbar;
		hCB = colorbar;
		hCB.Label.String='Adjusted RIM_1(\delta)'
		set(gca,'XScale','linear','YScale','linear')
		%titletag=sprintf('Ring %d-1-%d - %s',N,out,opt);
		%title(titletag); 
		xlabel('differential sensitivity at \delta=0');
		ylabel('perturbation strength \delta');
		set(gcf,'PaperSize',[6 4],'PaperPosition',[0 0 6 4]);
		set(findall(gcf,'-property','FontSize'),'FontSize',12);
		figtag = sprintf('../figures/heat_maps/heat_map_ring_%s_%d-%d',opt,N,out);
		saveas(gcf,figtag,'png');

		end
	end
end


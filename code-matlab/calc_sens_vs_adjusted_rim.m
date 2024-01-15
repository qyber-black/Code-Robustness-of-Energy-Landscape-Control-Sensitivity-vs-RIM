function [] = calc_sens_vs_adjusted_rim(noise_level)

% SPDX-FileCopyrightText: Copyright (C) 2023 Sean Patrick O'Neil <seanonei@usc.edu>
% SPDX-FileCopyrightText: Copyright (C) 2023 S Shermer <lw1660@gmail.com>, Swansea University
% SPDX-License-Identifier: CC-BY-SA-4.0 

% This function calculates and saves the relative error between the RIM at
% the smallest perturbation strength and the differential sensitivity as
% \delta = 0 for chains and rings

%noise_level = '0.1';

num   = 100;
index = 1:num;
in = 1;
type  = {'nmplus','ppo','snob'};
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
		    opt = type{q}
    
		    % load sensitivity and RIM data 
		    loadtag1 = sprintf('../results/log_sens-chains/log_sens_%s_%s_%d-chain_1-%d.mat',noise_level,opt,N,out);
		    load(loadtag1);
		    loadtag2 = sprintf('../results/kde-chains/log_kde_%s_%s_%d-chain_1-%d.mat',noise_level,opt,N,out);
		    load(loadtag2);
    
		    for k=1:100
			    rim_kde(k,:) = density{k}.mean;
		    end
	    
		    analytic_sens = abs(sens(:,1001));
    
		    % compute RIM based on kde data 
		    a_rim = rim_kde-rim_kde(:,1);
		    rim_sens = a_rim(:,2)*10000;
		    diff = 100*(analytic_sens - rim_sens)./analytic_sens;
		    suffix = sprintf('%s_%d_%d',opt,N,out);
		    rel_error.(suffix) = diff;
		    mean_rel_error.(suffix) = mean(diff);
		end
	end
end


option = {'dephasing';'fidelity';'overlap'};

for x = 1:3
    opt = option{x}
    
	for N = 5:6
	
		for out = 2:floor(N/2)+1
    
		    % load data 
		    loadtag1 = sprintf('../results/log_sens-rings/log_sens_%s_%d-ring_1-%d.mat',opt,N,out);
		    load(loadtag1);
		    loadtag2 = sprintf('../results/kde-rings/log_kde_%s_%d-ring_1-%d.mat',opt,N,out);
		    load(loadtag2);
		        
		    for k=1:100
    			rim_kde(k,:) = density{k}.mean;
    			kde_err(k,1) = density{k}.mean(1);
    		end
    
		    analytic_sens = abs(sens(:,1001));
    
		    % compute RIM based on kde data 
    		a_rim = rim_kde-rim_kde(:,1);
    		rim_sens = a_rim(:,2)*10000;
		    diff = 100*(analytic_sens - rim_sens)./analytic_sens;
    		suffix = sprintf('%s_%d_%d',opt,N,out);
    		rel_error.(suffix) = diff;
    		mean_rel_error.(suffix) = mean(diff);
    		clear 'rim_kde'
    
		end
	end
end

savetag = sprintf('../results/sens_v_rim_sens');
save(savetag,'rel_error','mean_rel_error');


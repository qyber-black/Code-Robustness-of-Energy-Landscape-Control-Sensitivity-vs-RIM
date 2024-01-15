%% plot composite sensitivity and adjusted RIM data for chains (plot_all_chain_v2.m)

% SPDX-FileCopyrightText: Copyright (C) 2023 Sean Patrick O'Neil <seanonei@usc.edu>
% SPDX-FileCopyrightText: Copyright (C) 2023 S Shermer <lw1660@gmail.com>, Swansea University
% SPDX-License-Identifier: CC-BY-SA-4.0 

function [] = plot_all_chains_v2(noise_level)

% This function takes as input the log-sensitivity data in the
% ../results/log_sens-chains/ and ../results/kde-chains/ and RIM data in the
% ../results/rim-chains/ and produces plots of the sensitivity and adjusted RIM
%
% Input: noise_level -- requires that log_sens, kde and rim data for the
% selected noise level has been computed.  Default 0.1
%
% On input:
% N - chain size 
% out - target spin
% opt - optimization option (nmplus, ppo, or snob)

% Output saved to ../figures/sens_composite/ directory 
% sens_comparison_chain_opt_N-out - plot of sensitity, adjusted RIM, and 
% fidelity error verus controller index 
% sens_scatter_chain_opt_N-out - scatter plot of sensitivity and adjusted RIM versus error

if ~exist('noise_level','var')
   noise_level = '0.1';
end
   
if ~exist('../figures/composite-chains-v2/','dir')
    mkdir('../figures/composite-chains-v2/');
end

num   = 100;
in    = 1;
type  = {'lbfgs','nmplus','ppo','snob'};
index = 1:100;

for N = 5:6

	if N == 5
        target = [3 5];
    else
        target = [4 6];
    end

    for x = 1:2
	    out = target(x);

    	for q = 2:4

			% Load data
		    opt = type{q};
    
		    loadtag1 = sprintf('../results/log_sens-chains/log_sens_%s_%s_%d-chain_1-%d.mat',noise_level,opt,N,out);
		    load(loadtag1);
		    loadtag2 = sprintf('../results/kde-chains/log_kde_%s_%s_%d-chain_1-%d.mat',noise_level,opt,N,out);
		    load(loadtag2);
		    loadtag3 = sprintf('../results/rim-chains/rim_noise_%s_%s_%d-chain_1-%d.mat',noise_level,opt,N,out);
		    load(loadtag3);

		    temp = arrayfun(@(n) density{n}.sensitivity(1),1:num,'UniformOutput',false);
		    log_sens_kde = cell2mat(temp');
		    log_sens_analytic = mean_log_sens;

		    % Compute sensitivity and adjusted RIM
		    sens_kde=log_sens_kde.*err;
		    sens_analytic=log_sens_analytic.*err;
		    rim = controller_rim-err;

		    % Plot Sensitivity, adjusted RIM, and Error versus Controller Index 
		    figure;
		    yyaxis left;
		    plot(index,log10(sens_analytic),'b','Linewidth',0.5);
		    hold on;
		    plot(index,log10(sens_kde),'r--','Linewidth',1);
		    plot(index,log10(rim),"-k");
		    xlabel('Controller Index');
		    ylabel('log(de(T)/d\delta) and log(RIM_1)'); 
		    yyaxis right;
		    plot(index,(err));
		    ylabel('e(T)');
			grid on;
			legend('sens_{analytic}','sens_{kde}','RIM_1','e(T)','location','best')
			titletag=sprintf('%s-optimized controllers %d-%d transfer',opt,N,out);
			title(titletag);
			set(gcf,'PaperSize',[6 4],'PaperPosition',[0 0 6 4]);
			set(findall(gcf,'-property','FontSize'),'FontSize',14);
			figtag = sprintf('../figures/composite-chains-v2/comparison_v2_%s_%s_%d-chain_1-%d.fig',noise_level,opt,N,out);
			savefig(figtag);
			saveas(gcf,sprintf('%s.png',figtag(1:end-4)),'png');

			% Plot Sensitivity and Adjusted RIM versus Error 
			figure;
			loglog(err,sens_analytic,'+',err,sens_kde,'square',err,rim,'*');
			grid on;
			xlabel('log(e(T))');
			ylabel('log(de(T)/d\delta and log(RIM_1)');
			legend('sens_{analytic}','sens_{kde}','RIM_1','location','best')
			titletag = sprintf('%s-optimized controllers %d-%d transfer',opt,N,out);
			title(titletag);
			set(gcf,'PaperSize',[6 4],'PaperPosition',[0 0 6 4]);
			set(findall(gcf,'-property','FontSize'),'FontSize',14);
			figtag = sprintf('../figures/composite-chains-v2/scatter_v2_%s_%s_%d-chain_1-%d.fig',noise_level,opt,N,out);
	    	savefig(figtag);
	    	saveas(gcf,sprintf('%s.png',figtag(1:end-4)),'png');
	    	close all;
	    end
    end
end


%% Plot composite data (log-sensitivity and RIM) for chains (plot_all_chain.m)

% SPDX-FileCopyrightText: Copyright (C) 2023 Sean Patrick O'Neil <seanonei@usc.edu>
% SPDX-FileCopyrightText: Copyright (C) 2024 S M Shermer <lw1660@gmail.com>
% SPDX-License-Identifier: CC-BY-SA-4.0 

function [] = plot_all_chain(noise_level)
%
% This function takes as input the log-sensitivity data in the
% ../results/log_sens-chains/ directory and RIM data in the
% ../results/rim_chains/ directory and produces plots of the results

% On input:
% N - chain size 
% out - target spin
% opt - optimization option (nmplus, ppo, or snob)

% On output - saved to ../figures/composite-chains/ directory 
% comparison_chain_opt_N-out - plot of log-sensitity, RIM, and 
% fidelity error verus controller index  
% scatter_chain_opt_N-out - scatter plot of log-sensitivity and RIM versus error 

% Input: noise_level -- requires that log_sens, kde and rim data for the
% selected noise level has been computed.  Default 0.1

if ~exist('noise_level','var')
   noise_level = '0.1';
end

if ~exist('../figures/composite-chains/','dir')
    mkdir('../figures/composite-chains/');
end

num=100;
in = 1;
type = {'lbfgs','nmplus','ppo','snob'};
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

			% Plot Log-Sensitivity and Error versus Controller Index 
			figure;
			yyaxis left;
			plot(index,log10(log_sens_analytic),'b','Linewidth',0.5);
			hold on;
			plot(index,log10(log_sens_kde),'r--','Linewidth',1);
			xlabel('Controller Index');
			ylabel('log(s_a(S,T)) and log(s_k(S,T))'); 
			yyaxis right;
			plot(index,err,index,controller_rim);
			ylabel('e(T) and RIM_1');
			grid on;
			set(gcf,'PaperSize',[6 4],'PaperPosition',[0 0 6 4]);
			set(findall(gcf,'-property','FontSize'),'FontSize',14);
			legend('s_a(S,T)','s_k(S,T)','e(T)','RIM_1','location','southeast')
			figtag = sprintf('../figures/composite-chains/comparison_%s_%s_%d-chain_1-%d.fig',noise_level,opt,N,out);
			savefig(figtag);
			saveas(gcf,sprintf('%s.png',figtag(1:end-4)),'png');

			% Plot Log-Sensitivity versus Error 
			figure;
			loglog(err,log_sens_analytic,'+',err,log_sens_kde,'square',err,controller_rim,'*');
			grid on;
			xlabel('log(e(T))');
			ylabel('log(s_{{a,k}}(S,T)) and log(RIM_1)');
			legend('s_a(S,T)','s_k(S,T)','RIM_1','location','best')
			set(gcf,'PaperSize',[6 4],'PaperPosition',[0 0 6 4]);
			set(findall(gcf,'-property','FontSize'),'FontSize',14);
			titletag = sprintf('%s-optimized controllers %d-%d transfer',opt,N,out);
			title(titletag);
			figtag = sprintf('../figures/composite-chains/scatter_%s_%s_%d-chain_1-%d.fig',noise_level,opt,N,out);
			savefig(figtag);
			saveas(gcf,sprintf('%s.png',figtag(1:end-4)),'png');
			close all;
		    
		end
	end
end



    


    

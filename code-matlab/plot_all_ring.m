%% Plot composite data (log-sensitivity and RIM) for rings (plot_all_ring.m)

% SPDX-FileCopyrightText: Copyright (C) 2023 Sean Patrick O'Neil <seanonei@usc.edu>
% SPDX-FileCopyrightText: Copyright (C) 2024 S M Shermer <lw1660@gmail.com>
% SPDX-License-Identifier: CC-BY-SA-4.0 

function [] = plot_all_ring()
%
% This function takes as input the log-sensitivity data in the
% ../results/log_sens-rings/ and ../results/kde_rings directory and RIM data in the
% ../results/rim-rings/ directory and produces plots of the results

% On input:
% N - chain size 
% out - target spin
% opt - optimization option (nmplus, ppo, or snob)

% On output - saved to ../figures/composite/ directory 
% comparison_ring_opt_N-out - plot of log-sensitity, RIM, and fidelity error verus controller index  
% scatter_ring_opt_N-out - scatter plot of log-sensitivity and RIM versus error 

in  = 1;
num = 1000;
tol = 10^-12;
option = {'dephasing';'fidelity';'overlap'};

if ~exist('../figures/composite-rings/','dir')
    mkdir('../figures/composite-rings/');
end

% Plot figures
index = 1:100;
option = {'dephasing';'fidelity';'overlap'};
for x = 1:3
    opt = option{x};
	for N = 5:6
	    for out = 2:floor(N/2)+1
    
		    tag1 = sprintf('../results/log_sens-rings/log_sens_%s_%d-ring_1-%d',opt,N,out);
		    load(tag1);
		    num = length(log_sens)-1;
		    log_sens_2 = log_sens(:,num+1);
		    tag2 = sprintf('../results/kde-rings/log_kde_%s_%d-ring_1-%d',opt,N,out);
		    load(tag2);
		    log_sens_1 = arrayfun(@(n) density{n}.sensitivity,1:100)';
		    err_1 = arrayfun(@(n) density{n}.mean(1),1:100)';
		    tag3 = sprintf('../results/rim-rings/rim_%s_%d-ring_1-%d.mat',opt,N,out);
		    load(tag3);
    
		    % Sort by error saved in Controller{n} to match sequence of analytical
		    % calculations 
		    Z = [err_1 log_sens_1];
		    Z = sortrows(Z);
		    err_1 = Z(:,1);
		    log_sens_1 = Z(:,2);
	    
		    % Plot Log-Sensitivity, RIM, and Error versus Controller Index 
		    figure;
		    yyaxis left;
		    plot(index,log10(log_sens_2),'b','Linewidth',0.5);
		    hold on;
		    plot(index,log10(log_sens_1),'r--','Linewidth',1);
		    plot(index,log10(controller_rim),"-k");
		    xlabel('Controller Index');
		    ylabel('log(s_{\{a,k\}}(S,T)) and log(RIM_1)'); 
		    yyaxis right;
		    plot(index,(err));
		    ylabel('e(T)');
		    grid on;
		    set(gcf,'PaperSize',[6 4],'PaperPosition',[0 0 6 4]);
		    set(findall(gcf,'-property','FontSize'),'FontSize',14);
		    legend('s_a(S,T)','s_k(S,T)','RIM_1','e(T)','location','best')
		    titletag=sprintf('%s-optimized controllers %d-%d transfer',opt,N,out);
		    title(titletag);
		    figtag = sprintf('../figures/composite-rings/comparison_%s_%d-ring_1-%d.fig',opt,N,out);
		    savefig(figtag);
		    saveas(gcf,sprintf('%s.png',figtag(1:end-4)),'png');

		    % Plot Log-Sensitivity versus Error 
		    figure;
		    loglog(err,log_sens_2,'+',err,log_sens_1,'square',err,controller_rim,'*');
		    grid on;
		    xlabel('log(e(T))');
		    ylabel('log(s_{\{a,k\}}(S,T) and log(RIM_1)');
		    legend('s_a(S,T)','s_k(S,T)','RIM_1','location','best')
		    set(gcf,'PaperSize',[6 4],'PaperPosition',[0 0 6 4]);
		    set(findall(gcf,'-property','FontSize'),'FontSize',14);
		    titletag = sprintf('%s-optimized controllers %d-%d transfer',opt,N,out);
		    title(titletag);
		    figtag = sprintf('../figures/composite-rings/scatter_%s_%d-ring_1-%d.fig',opt,N,out);
		    savefig(figtag);
		    saveas(gcf,sprintf('%s.png',figtag(1:end-4)),'png');
		    close all;
	    end
	end
end


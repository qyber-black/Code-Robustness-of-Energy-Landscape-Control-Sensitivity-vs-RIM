%% Plot composite data (sensitivity and adjusted RIM) for rings (plot_all_ring.m)

% SPDX-FileCopyrightText: Copyright (C) 2023 Sean Patrick O'Neil <seanonei@usc.edu>
% SPDX-FileCopyrightText: Copyright (C) 2024 S M Shermer <lw1660@gmail.com>
% SPDX-License-Identifier: CC-BY-SA-4.0 
%
function [] = plot_all_rings()
%
% This function takes as input the log-sensitivity data in the
% ../results/log_sens-rings/ and ../results/kde-rings directory and RIM data in the
% ../results/rim-rings/ directory and produces plots of the results for
% the sensitivity and adjusted RIM

% On input:
% N - chain size 
% out - target spin
% opt - optimization option (nmplus, ppo, or snob)

% On output - saved to ../figures/sens_composite/ directory 
% sens_comparison_ring_opt_N-out - plot of sensitity, adjusted RIM, and fidelity error verus controller index  
% sens_scatter_ring_opt_N-out - scatter plot of sensitivity and adjusted RIM versus error 

in  = 1;
num = 1000;
tol = 10^-12;
option = {'dephasing';'fidelity';'overlap'};

if ~exist('../figures/composite-rings-v2/','dir')
    mkdir('../figures/composite-rings-v2/');
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
    		sens_analytic = log_sens_2.*err;
    		sens_kde = log_sens_2.*err;
    		rim = controller_rim-err;

		    % Plot Sensitivity, RIM, and Error versus Controller Index 
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
		    figtag = sprintf('../figures/composite-rings-v2/comparison_v2_%s_%d-ring_1-%d.fig',opt,N,out);
		    savefig(figtag);
		    saveas(gcf,sprintf('%s.png',figtag(1:end-4)),'png');

		    % Plot Sensitivity versus Error 
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
		    figtag = sprintf('../figures/composite-rings-v2/scatter_v2_%s_%d-ring_1-%d.fig',opt,N,out);
		    savefig(figtag);
		    saveas(gcf,sprintf('%s.png',figtag(1:end-4)),'png');
		    close all;
    
	    end
	end
end


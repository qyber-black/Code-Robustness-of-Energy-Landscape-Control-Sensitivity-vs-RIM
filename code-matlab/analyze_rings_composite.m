%% Analyze log-sensitivity and RIM trends for rings -  (analyze_rings_composite.m)

function [] = analyze_rings_composite()

% SPDX-FileCopyrightText: Copyright (C) 2023 Sean Patrick O'Neil <seanonei@usc.edu>
% SPDX-FileCopyrightText: Copyright (C) 2023 Frank C Langbein <frank@langbein.org>
% SPDX-FileCopyrightText: Copyright (C) 2023 SM Shermer <lw1660@gmail.com>
% SPDX-License-Identifier: CC-BY-SA-4.0  
%
% This script takes as input the log-sensitivity data (analytic) in the 
% ../results/log_sens-rings/ and ../results/kde-rings and RIM data from the 
% ../results/rim-rings/and calcuates the correlations between the RIM, 
% log-sensitivity, and both versus the fidelity error. The output 
% is saved as a spreadsheet in the ../results/ directory
%
% On Input
% N - size of spin-ring
% out - target spin for transfer 
% option - optimization option (dephasing, fidelity, or dephasing)
%
% On output
% spreadsheet of correlation data organized in two sheets with the
% following columns:
%   "log_sens_v_rim":
%   1 - Kendall tau for log-sensitivity measures
%   2 - Kendall tau test statistic for (1)
%   3 - p-value for test statistic in (3) 
%   4 - Pearson r for log-sensitivity measures
%   5 - Pearson r test statistic for (4)
%   6 - p-value for test statistic in (5) 
%   7 - Kendall tau for log-sensitivity (analyitc) and RIM 
%   8 - Kendall tau test statistic for (7)
%   9 - p-value for test statistic in (8) 
%   10 - Pearson r for log-sensitivity (analytic) and RIM
%   11 - Pearson r test statistic for (10)
%   12 - p-value for test statistic in (11) 
%   "log_sens_and_rim_v_error"
%   1 - Kendall tau for log-sensitivity (kde) vs. error
%   2 - Kendall tau test statistic for (1)
%   3 - p-value for test statistic in (3) 
%   4 - Kendall tau for log-sensitivity (analytic) vs. error
%   5 - Kendall tau test statistic for (4)
%   6 - p-value for test statistic in (5) 
%   7 - Pearson r for log-sensitivity (kde) vs. error (log-log)
%   8 - Pearson r test statistic for (7)
%   9 - p-value for test statistic in (8) 
%   10 - Pearson r for log-sensitivity (analytic) vs. error (log-log)
%   11 - Pearson r test statistic for (10)
%   12 - p-value for test statistic in (11) 
%   13 - Kendall tau for RIM vs. error
%   14 - Kendall tau test statistic for (13)
%   15 - p-value for test statistic in (14) 
%   16 - Pearson r for RIM vs. error (log-log)
%   17 - Pearson r test statistic for (16)
%   18 - p-value for test statistic in (17)

num = 1000;
count = 0;
option = {'dephasing';'fidelity';'overlap'};

for x = 1:3
    opt = option{x}
    for N = 5:6
        for out = 2:floor(N/2)+1
            count = count+1;
            % load log-sensitivity data
            tag1 = sprintf('../results/log_sens-rings/log_sens_%s_%d-ring_1-%d.mat',opt,N,out);
            tag2 = sprintf('../results/kde-rings/log_kde_%s_%d-ring_1-%d.mat',opt,N,out);
            tag3 = sprintf('../results/controllers-rings/%s_%d-ring_1-%d.mat',opt,N,out);
            load(tag1); % load analytical results 
            load(tag2); % load data set 2 density data 
            load(tag3); % load data files for controller
            
            err_1 = arrayfun(@(n) 1-sys{n}.fidelity,1:100)';  
            test = arrayfun(@(n) density{n}.mean(1,1),1:100)';  
            log_sens_1 = arrayfun(@(n) density{n}.sensitivity,1:100)';
            log_sens_2 = log_sens(:,1001);

    		% sort error and log-sensitivity saved in Controller and Density
    		Z = [err_1 log_sens_1 log_sens_2];
			Z = sortrows(Z);
	        err_1 = Z(:,1);
	        log_sens_1 = abs(Z(:,2));
	        log_sens_2 = Z(:,3);
	        %log_sens_2 = abs(log_sens(:,num+1));

	        % load RIM data
	        loadtag3 = sprintf('../results/rim-rings/rim_%s_%d-ring_1-%d.mat',opt,N,out);
	        load(loadtag3);
    
	        % Calculate Kendall tau statistics for log_sens agreement
	        correl_1(count,1) = corr(log_sens_1,log_sens_2,'type','kendall');
	        n = length(err);
	        sigma_k = sqrt(2*(2*n+5)/(9*n*(n-1)));
	        correl_1(count,2) = correl_1(count,1)/sigma_k;
    
	        if correl_1(count,2) > 0
	           correl_1(count,3) = 1*(1-normcdf(correl_1(count,2)));
	        else
	            correl_1(count,3) = 1*normcdf(correl_1(count,2));
	        end
    
	        % Calculate Pearson r statistics for log_sens agreement on linear scale
	        correl_1(count,4) = corr((log_sens_2),(log_sens_1),'type','pearson');
	        sigma_p = sqrt(n-2);
	        correl_1(count,5) = correl_1(count,4)*sigma_p/sqrt(1 - correl_1(count,4)^2);
    
	        if correl_1(count,5) > 0
	           correl_1(count,6) = 1*(1-tcdf(correl_1(count,5),n-2));
	        else
	          correl_1(count,6) = 1*tcdf(correl_1(count,5),n-2);
	        end
    
	        % Calculate Kendall tau statistics for log_sens and rim
	        correl_1(count,7) = corr(log_sens_2,controller_rim,'type','kendall');
	        n = length(err);
	        sigma_k = sqrt(2*(2*n+5)/(9*n*(n-1)));
	        correl_1(count,8) = correl_1(count,7)/sigma_k;
    
	        if correl_1(count,8) > 0
	            correl_1(count,9) = 1*(1-normcdf(correl_1(count,8)));
	        else
	            correl_1(count,9) = 1*normcdf(correl_1(count,8));
	        end
    
	        % Calculate Pearson r statistics for log_sens and rim on linear scale
	        correl_1(count,10) = corr((log_sens_2),(controller_rim),'type','pearson');
	        sigma_p = sqrt(n-2);
	        correl_1(count,11) = correl_1(count,10)*sigma_p/sqrt(1 - correl_1(count,10)^2);
    
	        if correl_1(count,11) > 0
	            correl_1(count,12) = 1*(1-tcdf(correl_1(count,11),n-2));
	        else
	            correl_1(count,12) = 1*tcdf(correl_1(count,11),n-2);
	        end

	        % Calculate Kendall tau statistics for log_sens_1 (mc) versus error
	        correl_2(count,1) = corr(err,log_sens_1,'type','kendall');
	        sigma_k = sqrt(2*(2*n+5)/(9*n*(n-1)));
	        correl_2(count,2) = correl_2(count,1)/sigma_k;
	    
	        if correl_2(count,2) > 0
	        	correl_2(count,3) = 1*(1-normcdf(correl_2(count,2)));
	        else
            correl_2(count,3) = 1*normcdf(correl_2(count,2));
	        end
    
	        % Calculate Kendall tau statistics for log_sens_2 (analytic) versus error
	        correl_2(count,4) = corr(err,log_sens_2,'type','kendall');
	        sigma_k = sqrt(2*(2*n+5)/(9*n*(n-1)));
	        correl_2(count,5) = correl_2(count,4)/sigma_k;
    
	        if correl_2(count,5) > 0
		        correl_2(count,6) = 1*(1-normcdf(correl_2(count,5)));
		    else
		        correl_2(count,6) = 1*normcdf(correl_2(count,5));
		    end

		    % Calculate Pearson r statistics for log_sens_1 (mc) versus err on
		    % log-log scale
		    correl_2(count,7) = corr(log10(err),log10(log_sens_1),'type','pearson');
		    sigma_p = sqrt(n-2);
		    correl_2(count,8) = correl_2(count,7)*sigma_p/sqrt(1 - correl_2(count,7)^2);
    
		    if correl_2(count,8) > 0
		        correl_2(count,9) = 1*(1-tcdf(correl_2(count,8),n-2));
		    else
		        correl_2(count,9) = 1*tcdf(correl_2(count,8),n-2);
		    end

		    % Calculate Pearson r statistics for log_sens_2 (analytic) versus err on
		    % log-log scale
		    correl_2(count,10) = corr(log10(err),log10(log_sens_2),'type','pearson');
		    sigma_p = sqrt(n-2);
		    correl_2(count,11) = correl_2(count,10)*sigma_p/sqrt(1 - correl_2(count,10)^2);
    
		    if correl_2(count,11) > 0
		        correl_2(count,12) = 1*(1-tcdf(correl_2(count,11),n-2));
		    else
		        correl_2(count,12) = 1*tcdf(correl_2(count,11),n-2);
		    end

		    % Calculate Kendall tau statistics for rim v. error
		    correl_2(count,13) = corr(err,controller_rim,'type','kendall');
		    sigma_k = sqrt(2*(2*n+5)/(9*n*(n-1)));
		    correl_2(count,14) = correl_2(count,13)/sigma_k;
    
		    if correl_2(count,14) > 0
		        correl_2(count,15) = 1*(1-normcdf(correl_2(count,14)));
		    else
		        correl_2(count,15) = 1*normcdf(correl_2(count,14));
		    end

		    % Calculate Pearson r statistics for rim v. error
		    correl_2(count,16) = corr(log10(err),log10(controller_rim),'type','pearson');
		    sigma_p = sqrt(n-2);
		    correl_2(count,17) = correl_2(count,16)*sigma_p/sqrt(1 - correl_2(count,16)^2);
    
    		if correl_2(count,17) > 0
    		    correl_2(count,18) = 1*(1-tcdf(correl_2(count,17),n-2));
    		else
    		    correl_2(count,18) = 1*tcdf(correl_2(count,17),n-2);
    		end
    		rowtag = sprintf('N=%d out=%d %s',N,out,opt);
    		rowname{count,1} = rowtag;
	
    		close
		end
	end
end

headings_1 = {'tau - log_sens_mc v. mean(log_sens_2)','Z_tau1','p-value1','Pearson - log_sens_mc v. mean(log_sens_2)','T_2','p-value2','tau - log_sens v. rim','Z_tau3','p-value3','Pearson - log_sens v. rim','T_4','p-value4'};
headings_2 = {'tau - log_sens_mc v. err','Z_tau1','p-value1','tau - log_sens_2 v. err','Z_tau2','p-value2','Pearson - log10(log_sens_mc) v. log10(err)','T_3','p-value3','Pearson - log10(log_sens_2) v. log10(err)','T_4','p-value4','tau - rim v. err','Z_tau5','p-value5','Pearson - log10(rim) v. log10(err)','T_6','p-value6'};

% Create table and save data 
corr_data_1 = array2table(correl_1);
corr_data_1.Properties.RowNames = rowname;
corr_data_1.Properties.VariableNames = headings_1;
writetable(corr_data_1,'../results/corr-ring-log_sens_rim.csv','WriteRowNames',true);

corr_data_2 = array2table(correl_2);
corr_data_2.Properties.RowNames = rowname;
corr_data_2.Properties.VariableNames = headings_2;
writetable(corr_data_2,'../results/corr-ring-log_sens-rim-err.csv','WriteRowNames',true);

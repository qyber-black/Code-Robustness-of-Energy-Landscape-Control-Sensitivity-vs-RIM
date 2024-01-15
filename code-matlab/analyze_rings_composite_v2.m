%% Analyze sensitivity and adjusted RIM trends for rings (analyze_rings_composite_v2.m)

function [] = analyze_rings_composite_v2()

% SPDX-FileCopyrightText: Copyright (C) 2023 Sean Patrick O'Neil <seanonei@usc.edu>
% SPDX-FileCopyrightText: Copyright (C) 2023 Frank C Langbein <frank@langbein.org>
% SPDX-FileCopyrightText: Copyright (C) 2023 SM Shermer <lw1660@gmail.com>
% SPDX-License-Identifier: CC-BY-SA-4.0  

% This script takes as input the log-sensitivity data (analytic and kde) in the 
% ../results/log_sens_results/ directory and RIM data from the 
% ../results/rim_data_mat/ directory and calcuates the correlations between
% the adjusted RIM, sensitivity, and both versus the fidelity error. 
% The output is saved as a  spreadsheet in the ../results/ directory

% On output
% spreadsheet of correlation data organized in two sheets with the
% following columns:
%   "sens_v_rim":
%   1 - Kendall tau for sensitivity measures
%   2 - Kendall tau test statistic for (1)
%   3 - p-value for test statistic in (3) 
%   4 - Pearson r for sensitivity measures
%   5 - Pearson r test statistic for (4)
%   6 - p-value for test statistic in (5) 
%   7 - Kendall tau for sensitivity (analyitc) and RIM 
%   8 - Kendall tau test statistic for (7)
%   9 - p-value for test statistic in (8) 
%   10 - Pearson r for sensitivity (analytic) and RIM
%   11 - Pearson r test statistic for (10)
%   12 - p-value for test statistic in (11) 
%   "sens_and_rim_v_error"
%   1 - Kendall tau for sensitivity (kde) vs. error
%   2 - Kendall tau test statistic for (1)
%   3 - p-value for test statistic in (3) 
%   4 - Kendall tau for sensitivity (analytic) vs. error
%   5 - Kendall tau test statistic for (4)
%   6 - p-value for test statistic in (5) 
%   7 - Pearson r for sensitivity (kde) vs. error (log-log)
%   8 - Pearson r test statistic for (7)
%   9 - p-value for test statistic in (8) 
%   10 - Pearson r for sensitivity (analytic) vs. error (log-log)
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

    %load log-sensitivity data 
    tag1 = sprintf('../results/log_sens-rings/log_sens_%s_%d-ring_1-%d.mat',opt,N,out);
    tag2 = sprintf('../results/kde-rings/log_kde_%s_%d-ring_1-%d.mat',opt,N,out);
    tag3 = sprintf('../results/controllers-rings/%s_%d-ring_1-%d.mat',opt,N,out);
    load(tag1); % load analytical results 
    load(tag2); % load data set 2 density data  
    load(tag3); % load data files for controller
    
    err_1 = arrayfun(@(n) 1-sys{n}.fidelity,1:100)';  
    test = arrayfun(@(n) density{n}.mean(1,1),1:100)';  
    log_sens_1 = arrayfun(@(n) density{n}.sensitivity,1:100)';
   
    % sort error and log-sensitivity saved in Controller and Density
    Z = [err_1 log_sens_1 test];
    Z = sortrows(Z);
    err_1 = Z(:,1);
    log_sens_1 = Z(:,2);
    test = Z(:,3);
    log_sens_2 = log_sens(:,num+1);

    % load RIM data 
    tag4 = sprintf('../results/rim-rings/rim_%s_%d-ring_1-%d.mat',opt,N,out);
    load(tag4);
    
    % calcuate sensitiviyt and adjusted RIM
    sens_analytic = log_sens_2.*err;
    sens_kde = log_sens_2.*err;
    rim = controller_rim-err;
    
    % Calculate Kendall tau statistics for sens agreement
    correl_1(count,1) = corr(sens_kde,sens_analytic,'type','kendall');
    n = length(err);
    sigma_k = sqrt(2*(2*n+5)/(9*n*(n-1)));
    correl_1(count,2) = correl_1(count,1)/sigma_k;
    
    if correl_1(count,2) > 0
        correl_1(count,3) = 1*(1-normcdf(correl_1(count,2)));
    else
        correl_1(count,3) = 1*normcdf(correl_1(count,2));
    end
    
    % Calculate Pearson r statistics for sens agreement on linear scale
    correl_1(count,4) = corr((sens_analytic),(sens_kde),'type','pearson');
    sigma_p = sqrt(n-2);
    correl_1(count,5) = correl_1(count,4)*sigma_p/sqrt(1 - correl_1(count,4)^2);
    
    if correl_1(count,5) > 0
        correl_1(count,6) = 1*(1-tcdf(correl_1(count,5),n-2));
    else
        correl_1(count,6) = 1*tcdf(correl_1(count,5),n-2);
    end
    
    % Calculate Kendall tau statistics for sens and rim
    correl_1(count,7) = corr(sens_analytic,rim,'type','kendall');
    n = length(err);
    sigma_k = sqrt(2*(2*n+5)/(9*n*(n-1)));
    correl_1(count,8) = correl_1(count,7)/sigma_k;
    
    if correl_1(count,8) > 0
        correl_1(count,9) = 1*(1-normcdf(correl_1(count,8)));
    else
        correl_1(count,9) = 1*normcdf(correl_1(count,8));
    end
    
    % Calculate Pearson r statistics for sens and rim on linear scale
    correl_1(count,10) = corr((sens_analytic),(rim),'type','pearson');
    sigma_p = sqrt(n-2);
    correl_1(count,11) = correl_1(count,10)*sigma_p/sqrt(1 - correl_1(count,10)^2);
    
    if correl_1(count,11) > 0
        correl_1(count,12) = 1*(1-tcdf(correl_1(count,11),n-2));
    else
        correl_1(count,12) = 1*tcdf(correl_1(count,11),n-2);
    end

    % Calculate Kendall tau statistics for sens_kde versus error
    correl_2(count,1) = corr(err,sens_kde,'type','kendall');
    sigma_k = sqrt(2*(2*n+5)/(9*n*(n-1)));
    correl_2(count,2) = correl_2(count,1)/sigma_k;
    
    if correl_2(count,2) > 0
        correl_2(count,3) = 1*(1-normcdf(correl_2(count,2)));
    else
        correl_2(count,3) = 1*normcdf(correl_2(count,2));
    end
    
    % Calculate Kendall tau statistics for sens_analytic versus error
    correl_2(count,4) = corr(err,sens_analytic,'type','kendall');
    sigma_k = sqrt(2*(2*n+5)/(9*n*(n-1)));
    correl_2(count,5) = correl_2(count,4)/sigma_k;
    
    if correl_2(count,5) > 0
        correl_2(count,6) = 1*(1-normcdf(correl_2(count,5)));
    else
        correl_2(count,6) = 1*normcdf(correl_2(count,5));
    end

    % Calculate Pearson r statistics for sens_kde versus err on
    % log-log scale
    correl_2(count,7) = corr(log10(err),log10(sens_kde),'type','pearson');
    sigma_p = sqrt(n-2);
    correl_2(count,8) = correl_2(count,7)*sigma_p/sqrt(1 - correl_2(count,7)^2);
    
    if correl_2(count,8) > 0
        correl_2(count,9) = 1*(1-tcdf(correl_2(count,8),n-2));
    else
        correl_2(count,9) = 1*tcdf(correl_2(count,8),n-2);
    end

    % Calculate Pearson r statistics for sens_analytic versus err on
    % log-log scale
    correl_2(count,10) = corr(log10(err),log10(sens_analytic),'type','pearson');
    sigma_p = sqrt(n-2);
    correl_2(count,11) = correl_2(count,10)*sigma_p/sqrt(1 - correl_2(count,10)^2);
    
    if correl_2(count,11) > 0
        correl_2(count,12) = 1*(1-tcdf(correl_2(count,11),n-2));
    else
        correl_2(count,12) = 1*tcdf(correl_2(count,11),n-2);
    end

    % Calculate Kendall tau statistics for rim v. error
    correl_2(count,13) = corr(err,rim,'type','kendall');
    sigma_k = sqrt(2*(2*n+5)/(9*n*(n-1)));
    correl_2(count,14) = correl_2(count,13)/sigma_k;
    
    if correl_2(count,14) > 0
        correl_2(count,15) = 1*(1-normcdf(correl_2(count,14)));
    else
        correl_2(count,15) = 1*normcdf(correl_2(count,14));
    end

    % Calculate Pearson r statistics for rim v. error
    correl_2(count,16) = corr(log10(err),log10(rim),'type','pearson');
    sigma_p = sqrt(n-2);
    correl_2(count,17) = correl_2(count,16)*sigma_p/sqrt(1 - correl_2(count,16)^2);
    
    if correl_2(count,17) > 0
        correl_2(count,18) = 1*(1-tcdf(correl_2(count,17),n-2));
    else
        correl_2(count,18) = 1*tcdf(correl_2(count,17),n-2);
    end
    rowtag = sprintf('N=%d out=%d %s',N,out,opt);
    rowname{count,1} = rowtag;

    clear log_sens;
    close;

end
end
end

headings_1 = {'tau - sens_kde v. _sens_analytic','Z_tau1','p-value1','Pearson - sens_kde v. sens_analytic','T_2','p-value2','tau - sens v. rim','Z_tau3','p-value3','Pearson - sens v. rim','T_4','p-value4'};
headings_2 = {'tau - sens_kde v. err','Z_tau1','p-value1','tau - sens_analytic v. err','Z_tau2','p-value2','Pearson - log10(sens_kde) v. log10(err)','T_3','p-value3','Pearson - log10(sens_analytic) v. log10(err)','T_4','p-value4','tau - rim v. err','Z_tau5','p-value5','Pearson - log10(rim) v. log10(err)','T_6','p-value6'};

% Create table and save data 
corr_data_1 = array2table(correl_1);
corr_data_1.Properties.RowNames = rowname;
corr_data_1.Properties.VariableNames = headings_1;
writetable(corr_data_1,'../results/corr-ring-sens_adj_rim.csv','WriteRowNames',true);

corr_data_2 = array2table(correl_2);
corr_data_2.Properties.RowNames = rowname;
corr_data_2.Properties.VariableNames = headings_2;
writetable(corr_data_2,'../results/corr-ring-sens-adj_rim-err.csv','WriteRowNames',true);


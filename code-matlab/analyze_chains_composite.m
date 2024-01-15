%% Analyze log-sensitivity and RIM trends for chains -  (analyze_chains_composite.m)

function [] = analyze_chains_composite(noise_level)

% SPDX-FileCopyrightText: Copyright (C) 2023 Sean Patrick O'Neil <seanonei@usc.edu>
% SPDX-FileCopyrightText: Copyright (C) 2023 Frank C Langbein <frank@langbein.org>
% SPDX-FileCopyrightText: Copyright (C) 2023 SM Shermer <lw1660@gmail.com>
% SPDX-License-Identifier: CC-BY-SA-4.0  

% This script takes as input the log-sensitivity data (analytic and kde) in 
% ../results/log_sens-chains/ and ../results/kde-chains and RIM data from the 
% ../results/rim-chains/ directory and calcuates the correlations between
% the RIM, log-sensitivity, and both versus the fidelity error. The output is saved as a  
% spreadsheet in the ../results/ directory

% Input: noise level from the following list of options
% {'0.0','0.01','0.02','0.03','0.04','0.05','0.1'};
% default if not specified '0.1'

% On input:
% N   - chain size 
% out - target spin
% opt - optimization option (nmplus, ppo, or snob)

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
%   1 - Kendall tau for log-sensitivity (analytic) vs. error
%   2 - Kendall tau test statistic for (1)
%   3 - p-value for test statistic in (3) 
%   4 - Kendall tau for log-sensitivity (kde) vs. error
%   5 - Kendall tau test statistic for (4)
%   6 - p-value for test statistic in (5) 
%   7 - Pearson r for log-sensitivity (analytic) vs. error (log-log)
%   8 - Pearson r test statistic for (7)
%   9 - p-value for test statistic in (8) 
%   10 - Pearson r for log-sensitivity (kde) vs. error (log-log)
%   11 - Pearson r test statistic for (10)
%   12 - p-value for test statistic in (11) 
%   13 - Kendall tau for RIM vs. error
%   14 - Kendall tau test statistic for (13)
%   15 - p-value for test statistic in (14) 
%   16 - Pearson r for RIM vs. error (log-log)
%   17 - Pearson r test statistic for (16)
%   18 - p-value for test statistic in (17)


noise={'0.0','0.01','0.02','0.03','0.04','0.05','0.1'};
if ~exist('noise_level','var') || all(cellfun(@(x) strcmp(x,noise_level), noise))
    noise_level = '0.1'
end

disp(sprintf('noise level selected %s',noise_level))

if strcmp(noise_level,"0.0")
    q0 = 1; % include l-bfgs-controllers
else
    q0 = 2; % no l-bfgs controllers
end

c   = 1;
num = 100;
in  = 1;
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

        for q = q0:length(type)
            opt = type{q};

            disp(sprintf('running noise level %s %s controllers for chain N=%d target %d',noise_level,opt,N,out))
  
            % load data for correlation calculations 
            loadtag1 = sprintf('../results/log_sens-chains/log_sens_%s_%s_%d-chain_1-%d.mat',noise_level,opt,N,out);
            load(loadtag1);
            loadtag2 = sprintf('../results/kde-chains/log_kde_%s_%s_%d-chain_1-%d.mat',noise_level,opt,N,out);
            load(loadtag2);
            loadtag3 = sprintf('../results/rim-chains/rim_noise_%s_%s_%d-chain_1-%d.mat',noise_level,opt,N,out);
            load(loadtag3);
            
            temp = arrayfun(@(n) density{n}.sensitivity(1),1:num,'UniformOutput',false);
            log_sens_kde = cell2mat(temp');
            log_sens_analytic = mean_log_sens;
           
            % Calculate Kendall tau statistics for log_sens agreement
            correl_1(c,1) = corr(log_sens_analytic,log_sens_kde,'type','kendall');
            n = length(err);
            sigma_k = sqrt(2*(2*n+5)/(9*n*(n-1)));
            correl_1(c,2) = correl_1(c,1)/sigma_k;
    
            if correl_1(c,2) > 0
              correl_1(c,3) = 1*(1-normcdf(correl_1(c,2)));
            else
              correl_1(c,3) = 1*normcdf(correl_1(c,2));
            end
    
            % Calculate Pearson r statistics for log_sens agreement on linear scale
            correl_1(c,4) = corr((log_sens_kde),(log_sens_analytic),'type','pearson');
            sigma_p = sqrt(n-2);
            correl_1(c,5) = correl_1(c,4)*sigma_p/sqrt(1 - correl_1(c,4)^2);
    
            if correl_1(c,5) > 0
                correl_1(c,6) = 1*(1-tcdf(correl_1(c,5),n-2));
            else
                correl_1(c,6) = 1*tcdf(correl_1(c,5),n-2);
            end

            % Calculate Kendall tau statistics for log_sens and rim
            correl_1(c,7) = corr(log_sens_analytic,controller_rim,'type','kendall');
            n = length(err);
            sigma_k = sqrt(2*(2*n+5)/(9*n*(n-1)));
            correl_1(c,8) = correl_1(c,7)/sigma_k;
    
            if correl_1(c,8) > 0
                correl_1(c,9) = 1*(1-normcdf(correl_1(c,8)));
            else
                correl_1(c,9) = 1*normcdf(correl_1(c,8));
            end

            % Calculate Pearson r statistics for log_sens ard rim on linear scale
            correl_1(c,10) = corr((controller_rim),(log_sens_analytic),'type','pearson');
            sigma_p = sqrt(n-2);
            correl_1(c,11) = correl_1(c,10)*sigma_p/sqrt(1 - correl_1(c,10)^2);
    
            if correl_1(c,11) > 0
                correl_1(c,12) = 1*(1-tcdf(correl_1(c,11),n-2));
            else
                correl_1(c,12) = 1*tcdf(correl_1(c,11),n-2);
            end
    
            % Calculate Kendall tau statistics for log_sens_analytic versus error
            correl_2(c,1) = corr(err,log_sens_analytic,'type','kendall');
            sigma_k = sqrt(2*(2*n+5)/(9*n*(n-1)));
            correl_2(c,2) = correl_2(c,1)/sigma_k;
    
            if correl_2(c,2) > 0
                correl_2(c,3) = 1*(1-normcdf(correl_2(c,2)));
            else
                correl_2(c,3) = 1*normcdf(correl_2(c,2));
            end
    
            % Calculate Kendall tau statistics for log_sens_kde versus error
            correl_2(c,4) = corr(err,log_sens_kde,'type','kendall');
            sigma_k = sqrt(2*(2*n+5)/(9*n*(n-1)));
            correl_2(c,5) = correl_2(c,4)/sigma_k;
    
            if correl_2(c,5) > 0
                correl_2(c,6) = 1*(1-normcdf(correl_2(c,5)));
            else
                correl_2(c,6) = 1*normcdf(correl_2(c,5));
            end

            % Calculate Pearson r statistics for log_sens_analytic versus err on
            % log-log scale
            correl_2(c,7) = corr(log10(err),log10(log_sens_analytic),'type','pearson');
            sigma_p = sqrt(n-2);
            correl_2(c,8) = correl_2(c,7)*sigma_p/sqrt(1 - correl_2(c,7)^2);
    
            if correl_2(c,8) > 0
                correl_2(c,9) = 1*(1-tcdf(correl_2(c,8),n-2));
            else
                correl_2(c,9) = 1*tcdf(correl_2(c,8),n-2);
            end

            % Calculate Pearson r statistics for log_sens_kde versus err on
            % log-log scale
            correl_2(c,10) = corr(log10(err),log10(log_sens_kde),'type','pearson');
            sigma_p = sqrt(n-2);
            correl_2(c,11) = correl_2(c,10)*sigma_p/sqrt(1 - correl_2(c,10)^2);
    
            if correl_2(c,11) > 0
                correl_2(c,12) = 1*(1-tcdf(correl_2(c,11),n-2));
            else
                correl_2(c,12) = 1*tcdf(correl_2(c,11),n-2);
            end

            % Calculate Kendall tau statistics for rim versus error
            correl_2(c,13) = corr(err,controller_rim,'type','kendall');
            sigma_k = sqrt(2*(2*n+5)/(9*n*(n-1)));
            correl_2(c,14) = correl_2(c,13)/sigma_k;
            if correl_2(c,14) > 0
                correl_2(c,15) = 1*(1-normcdf(correl_2(c,14)));
            else
                correl_2(c,15) = 1*normcdf(correl_2(c,14));
            end

            % Calculate Pearson r statistics for rim vs. error
            correl_2(c,16) = corr(log10(err),log10(controller_rim),'type','pearson');
            sigma_p = sqrt(n-2);
            correl_2(c,17) = correl_2(c,16)*sigma_p/sqrt(1 - correl_2(c,16)^2);
    
            if correl_2(c,17) > 0
                correl_2(c,18) = 1*(1-tcdf(correl_2(c,17),n-2));
            else
                correl_2(c,18) = 1*tcdf(correl_2(c,17),n-2);
            end
            rowtag = sprintf('N=%d out=%d %s',N,out,opt);
            rowname{c,1} = rowtag;

            clear log_sens;
            close;
            c = c+1;
        end
    end
end

%% Generate Spreadsheets
headings_1 = {'tau - log_sens_kde v. _log_sens_analytic','Z_tau1','p-value1','Pearson - log_sens_kde v. log_sens_analytic','T_2','p-value2','tau - log_sens v. rim','Z_tau3','p-value3','Pearson - log_sens v. rim','T_4','p-value4'};
headings_2 = {'tau - log_sens_kde v. err','Z_tau1','p-value1','tau - log_sens_analytic v. err','Z_tau2','p-value2','Pearson - log10(log_sens_kde) v. log10(err)','T_3','p-value3','Pearson - log10(log_sens_analytic) v. log10(err)','T_4','p-value4','tau - rim v. err','Z_t','p-value5','Pearson - rim vs error','T_5','p-value6'}

% Create table and save data 
corr_data_1 = array2table(correl_1);
corr_data_1.Properties.RowNames = rowname;
corr_data_1.Properties.VariableNames = headings_1;
writetable(corr_data_1,sprintf('../results/corr-chain-noise_%s-log_sens-rim.csv',noise_level),'WriteRowNames',true);

corr_data_2 = array2table(correl_2);
corr_data_2.Properties.RowNames = rowname;
corr_data_2.Properties.VariableNames = headings_2;
writetable(corr_data_2,sprintf('../results/corr-chain-noise_%s-log_sens-rim-err.csv',noise_level),'WriteRowNames',true);

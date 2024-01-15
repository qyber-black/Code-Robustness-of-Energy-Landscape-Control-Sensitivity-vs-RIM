%% Analyze sensitivity and adjusted RIM trends for chains (analyze_chains_composite_v2.m)

function [] = analyze_chains_composite_v2(noise_level)

% SPDX-FileCopyrightText: Copyright (C) 2023 Sean Patrick O'Neil <seanonei@usc.edu>
% SPDX-FileCopyrightText: Copyright (C) 2023 Frank C Langbein <frank@langbein.org>
% SPDX-FileCopyrightText: Copyright (C) 2023 SM Shermer <lw1660@gmail.com>
% SPDX-License-Identifier: CC-BY-SA-4.0  

% This script takes as input the log-sensitivity data (analytic and kde) in 
% ../results/log_sens-chains/ and ../results/kde-chains and RIM data from 
% ../results/rim-chains/ directory and calcuates the correlations between
% the adjusted RIM, sensitivity, and the fidelity error. The output is saved 
% as a spreadsheet in the ../results/ directory

% Input: noise level from the following list of options
% {'0.0','0.01','0.02','0.03','0.04','0.05','0.1'};
% default if not specified '0.1'

% On Input
% N - size of spin-ring
% out - target spin for transfer 
% option - optimization option (nmplus, ppo, or snob)

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
%   7 - Kendall tau for sensitivity (analyitc) and adjusted RIM 
%   8 - Kendall tau test statistic for (7)
%   9 - p-value for test statistic in (8) 
%   10 - Pearson r for sensitivity (analytic) and adjusted RIM
%   11 - Pearson r test statistic for (10)
%   12 - p-value for test statistic in (11) 
%   "sens_and_rim_v_error"
%   1 - Kendall tau for sensitivity (analytic) vs. error
%   2 - Kendall tau test statistic for (1)
%   3 - p-value for test statistic in (3) 
%   4 - Kendall tau for sensitivity (kde) vs. error
%   5 - Kendall tau test statistic for (4)
%   6 - p-value for test statistic in (5) 
%   7 - Pearson r for sensitivity (analytic) vs. error (log-log)
%   8 - Pearson r test statistic for (7)
%   9 - p-value for test statistic in (8) 
%   10 - Pearson r for sensitivity (kde) vs. error (log-log)
%   11 - Pearson r test statistic for (10)
%   12 - p-value for test statistic in (11) 
%   13 - Kendall tau for adjusted RIM vs. error
%   14 - Kendall tau test statistic for (13)
%   15 - p-value for test statistic in (14) 
%   16 - Pearson r for adjusted RIM vs. error (log-log)
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

c     = 1;
num   = 100;
in    = 1;
noise = {'0.1'};
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
    
            noise_level = noise{1};
            disp(sprintf('running noise level %s %s controllers for chain N=%d target %d',noise_level,opt,N,out))
            loadtag1 = sprintf('../results/log_sens-chains/log_sens_%s_%s_%d-chain_1-%d.mat',noise_level,opt,N,out);
            load(loadtag1);
            loadtag2 = sprintf('../results/kde-chains/log_kde_%s_%s_%d-chain_1-%d.mat',noise_level,opt,N,out);
            load(loadtag2);
            loadtag3 = sprintf('../results/rim-chains/rim_noise_%s_%s_%d-chain_1-%d.mat',noise_level,opt,N,out);
            load(loadtag3);

            % extract sensitivity data from log-sensitivity
            temp = arrayfun(@(n) density{n}.sensitivity(1),1:num,'UniformOutput',false);
            log_sens_kde = cell2mat(temp');
            log_sens_analytic = mean_log_sens;
            sens_kde = log_sens_kde.*err;
            sens_analytic = log_sens_analytic.*err;
            % calcuate adjusted RIM
            rim = controller_rim-err;
           
            % Calculate Kendall tau statistics for sens agreement
            correl_1(c,1) = corr(sens_analytic,sens_kde,'type','kendall');
            n = length(err);
            sigma_k = sqrt(2*(2*n+5)/(9*n*(n-1)));
            correl_1(c,2) = correl_1(c,1)/sigma_k;
    
            if correl_1(c,2) > 0
              correl_1(c,3) = 1*(1-normcdf(correl_1(c,2)));
            else
              correl_1(c,3) = 1*normcdf(correl_1(c,2));
            end
    
            % Calculate Pearson r statistics for sens agreement on linear scale
            correl_1(c,4) = corr((sens_kde),(sens_analytic),'type','pearson');
            sigma_p = sqrt(n-2);
            correl_1(c,5) = correl_1(c,4)*sigma_p/sqrt(1 - correl_1(c,4)^2);
    
            if correl_1(c,5) > 0
               correl_1(c,6) = 1*(1-tcdf(correl_1(c,5),n-2));
            else
            correl_1(c,6) = 1*tcdf(correl_1(c,5),n-2);
        end

        % Calculate Kendall tau statistics for sens and rim
        correl_1(c,7) = corr(sens_analytic,rim,'type','kendall');
        n = length(err);
        sigma_k = sqrt(2*(2*n+5)/(9*n*(n-1)));
        correl_1(c,8) = correl_1(c,7)/sigma_k;
    
        if correl_1(c,8) > 0
           correl_1(c,9) = 1*(1-normcdf(correl_1(c,8)));
        else
           correl_1(c,9) = 1*normcdf(correl_1(c,8));
        end

        % Calculate Pearson r statistics for sens ard rim on linear scale
        correl_1(c,10) = corr((rim),(sens_analytic),'type','pearson');
        sigma_p = sqrt(n-2);
        correl_1(c,11) = correl_1(c,10)*sigma_p/sqrt(1 - correl_1(c,10)^2);
    
        if correl_1(c,11) > 0
           correl_1(c,12) = 1*(1-tcdf(correl_1(c,11),n-2));
        else
           correl_1(c,12) = 1*tcdf(correl_1(c,11),n-2);
        end
    
        % Calculate Kendall tau statistics for sens_analytic versus error
        correl_2(c,1) = corr(err,sens_analytic,'type','kendall');
        sigma_k = sqrt(2*(2*n+5)/(9*n*(n-1)));
        correl_2(c,2) = correl_2(c,1)/sigma_k;
    
        if correl_2(c,2) > 0
           correl_2(c,3) = 1*(1-normcdf(correl_2(c,2)));
        else
           correl_2(c,3) = 1*normcdf(correl_2(c,2));
        end
    
        % Calculate Kendall tau statistics for sens_kde versus error
        correl_2(c,4) = corr(err,sens_kde,'type','kendall');
        sigma_k = sqrt(2*(2*n+5)/(9*n*(n-1)));
        correl_2(c,5) = correl_2(c,4)/sigma_k;
    
        if correl_2(c,5) > 0
           correl_2(c,6) = 1*(1-normcdf(correl_2(c,5)));
        else
          correl_2(c,6) = 1*normcdf(correl_2(c,5));
        end

        % Calculate Pearson r statistics for sens_analytic versus err on
        % log-log scale
        correl_2(c,7) = corr(log10(err),log10(sens_analytic),'type','pearson');
        sigma_p = sqrt(n-2);
        correl_2(c,8) = correl_2(c,7)*sigma_p/sqrt(1 - correl_2(c,7)^2);
    
        if correl_2(c,8) > 0
           correl_2(c,9) = 1*(1-tcdf(correl_2(c,8),n-2));
        else
           correl_2(c,9) = 1*tcdf(correl_2(c,8),n-2);
        end

        % Calculate Pearson r statistics for sens_kde versus err on
        % log-log scale
        correl_2(c,10) = corr(log10(err),log10(sens_kde),'type','pearson');
        sigma_p = sqrt(n-2);
        correl_2(c,11) = correl_2(c,10)*sigma_p/sqrt(1 - correl_2(c,10)^2);
    
        if correl_2(c,11) > 0
          correl_2(c,12) = 1*(1-tcdf(correl_2(c,11),n-2));
        else
          correl_2(c,12) = 1*tcdf(correl_2(c,11),n-2);
        end

        % Calculate Kendall tau statistics for rim versus error
        correl_2(c,13) = corr(err,rim,'type','kendall');
        sigma_k = sqrt(2*(2*n+5)/(9*n*(n-1)));
        correl_2(c,14) = correl_2(c,13)/sigma_k;
        if correl_2(c,14) > 0
          correl_2(c,15) = 1*(1-normcdf(correl_2(c,14)));
        else
          correl_2(c,15) = 1*normcdf(correl_2(c,14));
        end

        % Calculate Pearson r statistics for rim vs. error
        correl_2(c,16) = corr(log10(err),log10(rim),'type','pearson');
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

%% Generate spreadsheets
headings_1 = {'tau - sens_analytic v. sens_kde','Z_tau1','p-value1','Pearson - sens_analytic v. sens_kde','T_2','p-value2','tau - sens_analytic v. rim','Z_tau3','p-value3','Pearson - sens_analytic v. rim','T_3','p-value4'};
headings_2 = {'tau - sens_analytic v. err','Z_tau1','p-value1','tau - log_sens_kde v. err','Z_tau2','p-value2','Pearson - log10(sens_analytic) v. log10(err)','T_3','p-value3','Pearson - log10(sens_kde) v. log10(err)','T_4','p-value4','tau - rim v. err','Z_tau5','p-value5','Pearson - log10(rim) v. log10(err)','T_6','p-value6'};

% Create table and save data 
corr_data_1 = array2table(correl_1);
corr_data_1.Properties.RowNames = rowname;
corr_data_1.Properties.VariableNames = headings_1;
writetable(corr_data_1,sprintf('../results/corr-chain-noise_%s-sens-adj_rim.csv',noise_level),'WriteRowNames',true);

corr_data_2 = array2table(correl_2);
corr_data_2.Properties.RowNames = rowname;
corr_data_2.Properties.VariableNames = headings_2;
writetable(corr_data_2,sprintf('../results/corr-chain-noise_%s-sens-adj_rim-err.csv',noise_level),'WriteRowNames',true);


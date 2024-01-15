% SPDX-FileCopyrightText: Copyright (C) 2022 Frank C Langbein <frank@langbein.org>, Cardiff University
% SPDX-FileCopyrightText: Copyright (C) 2022 S Shermer <lw1660@gmail.com>, Swansea University
% SPDX-FileCopyrightText: Copyright (C) 2022 Sean P O'Neil <seanonei@usc.edu>, University of Southern California
% 
% SPDX-License-Identifier: CC-BY-SA-4.0

function density = AnalyzeErrorDensity(Err,Ctrl)

% Analyze error distributions 
%
% Input
%   Err - array of errors for a given controller produced by the
%   AnalyseRobustnessError.m script
%   Ctrl - controller data 
%   Err(d,s) - error for fixed controller with uncertainty d of strength s
%
% Output
%   density - cell array with statiscal information (mean and variance) of
%   the distribution for the error for each controller - also includes a
%   best fit curve 
%   plots - plots depicting the error density and saved in ../figures
%   directory 

% General

[d_n,s_n] = size(Err);
min_err = min(min(Err));
max_err = max(max(Err));

% Mean and variance per strength

Mean = mean(Err);
Variance = var(Err);

% Kernel density estimation for each column (over decoherence operators at fixed strength)

E = reshape(Err,numel(Err),1);
bin_size = 3.5 * std(E) * numel(E)^(-1/3); % scott
bin_size = max([bin_size,(max_err-min_err)/1000,0.001]);
XI = min_err:bin_size:max_err;
bins = numel(XI);

Density = zeros(s_n,bins);
for s = 1:s_n
  Density(s,:) = ksdensity(Err(:,s),XI,'Bandwidth',bin_size,'Support',[min_err-bin_size*0.001 max_err+bin_size*0.001],'BoundaryCorrection','reflection');
end

% Sensitivity calculation - adjusted for maximum \delta of 0.1 
X = [0:s_n-1]/(s_n-1)*0.1;
[mean_fit,mean_gof] = fit(X',Mean','smoothingspline');
sens(1,1) = differentiate(mean_fit,0) / Mean(1);
%sens(2,1) = differentiate(mean_fit,0.001) / Mean(2);
%sens(3,1) = differentiate(mean_fit,0.01) / Mean(11);
%sens(4,1) = differentiate(mean_fit,0.1) / Mean(101);
%sens(5,1) = differentiate(mean_fit,0.5) / Mean(501);

err_kde(1,1) = Mean(1);
%err_kde(2,1) = Mean(2);
%err_kde(3,1) = Mean(11);
%err_kde(4,1) = Mean(101);
%err_kde(5,1) = Mean(501);


% Total distribution

[TotalDensity TotalXI] = ksdensity(E,'Bandwidth',bin_size,'Support',[min_err-bin_size*0.001 max_err+bin_size*0.001],'BoundaryCorrection','reflection');
TotalMean = mean(E);
TotalVariance = var(E);

density.mean = Mean;
density.var = Variance;
density.max_err = max_err;
density.min_err = min_err;
density.delta_kde = Density;
density.delta_xi = XI;
density.total_kde = TotalDensity;
density.total_xi = TotalXI;
density.total_mean = TotalMean;
density.total_var = TotalVariance;
density.sensitivity = sens;
density.err_kde = err_kde;


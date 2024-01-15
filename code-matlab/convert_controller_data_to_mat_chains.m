%% Convert rim controller data to .mat files for chains

function [] = convert_controller_data_to_mat_chains(noise_level)
%
% SPDX-FileCopyrightText: Copyright (C) 2023 Sean Patrick O'Neil <seanonei@usc.edu>
% SPDX-FileCopyrightText: Copyright (C) 2024 S M Shermer <lw1660@gmail.com>
% SPDX-License-Identifier: CC-BY-SA-4.0 
%
% This function takes the .csv data for the RIM controllers for chains from 
% ../data-raw/controller-chains/ and prodcues .mat files for use with the
% other routines in this data set.  RIM controller data is saved in the
% ../results/controllers-chains/ directory 
%
% Input: noise level from the following list of options
% {'0.0','0.01','0.02','0.03','0.04','0.05','0.1'};
% default if not specified '0.1'

if ~exist('../results/controllers-chains','dir')
    mkdir('../results/controllers-chains');
end

in    = 1;
noise = {'0.0','0.01','0.02','0.03','0.04','0.05','0.1'};
type  = {'lbfgs','nmplus','ppo','snob'};
x = 1;

if ~exist('noise_level','var') || all(cellfun(@(x) strcmp(x,noise_level), noise))
    noise_level = '0.1'
end

for N = 5:6
    if N == 5
        target = [3 5];
    else
        target = [4 6];
    end
    for x = 1:2
	    out = target(x);
    	if strcmp(noise_level,'0.0')
    	    q = 1:4;
        else    
            q = 2:4;
	        
        end    
        for q = min(q):max(q)
            opt = type{q};
	        tag = sprintf('../data-raw/controllers-chains/train_noise_%s/%s_Nspin_%d_outspin_%d.csv',noise_level,opt,N,out-1);
			controller_data = readmatrix(tag);
			controller_data(:,1) = [];
			controller_data(1,:) = [];
			controller_data_sorted = sortrows(controller_data,N+2,'descend');
			sys = convert_controller(controller_data_sorted,in,out,opt,N,x);
			savetag = sprintf('../results/controllers-chains/noise_%s_%s_%d-chain_1-%d.mat',noise_level,opt,N,out);
			save(savetag,'sys');
    	end
	end
end


%% Convert rim and error data to .mat files for chains
%
function [] = convert_rim_raw_data_to_mat_chains()
%
% SPDX-FileCopyrightText: Copyright (C) 2023 Sean Patrick O'Neil <seanonei@usc.edu>
% SPDX-FileCopyrightText: Copyright (C) 2024 S M Shermer <lw1660@gmail.com>
% SPDX-License-Identifier: CC-BY-SA-4.0 

% This function takes the .csv data from ../data-raw/rim-chains/ directory 
% and prodcues .mat files for use with the other routines in this data set. 
% The data is saved in the ../results/rim-chains/


if ~exist('../results/rim-chains','dir')
    mkdir('../results/rim-chains');
end

in    = 1;
noise = {'0.0','0.01','0.02','0.03','0.04','0.05','0.1'};
type  = {'lbfgs','nmplus','ppo','snob'};


for N = 5:6
    if N == 5
        target = [3 5];
    else
        target = [4 6];
    end
    for x = 1:2
        out = target(x);
	    for y = 1:length(noise)
	        noise_level = noise{y};
	        if y == 1
	            q = 1:4;
	        else
	            q = 2:4;
	        end
	        for q = min(q):max(q)
                opt = type{q};

				tag = sprintf('../data-raw/rim-chains/tn_%s_%s_dephasing_rim_1s_trace_chain_%d_%d.csv',noise_level,opt,N,out);
                disp(tag)
				rim_data = readmatrix(tag);
				controller_rim = rim_data(51,:)';

				savetag = sprintf('../results/rim-chains/rim_noise_%s_%s_%d-chain_1-%d.mat',noise_level,opt,N,out);
				save(savetag,'controller_rim','rim_data');
				clear controller_rim rim_data;
	        end
	    end
	end
end


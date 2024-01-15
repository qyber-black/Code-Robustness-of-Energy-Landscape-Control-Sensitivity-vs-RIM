%% Convert rim and error data to .mat files for rings
%
function [] = convert_rim_raw_data_to_mat_rings()
%
% SPDX-FileCopyrightText: Copyright (C) 2023 Sean Patrick O'Neil <seanonei@usc.edu>
% SPDX-FileCopyrightText: Copyright (C) 2024 S M Shermer <lw1660@gmail.com>
% SPDX-License-Identifier: CC-BY-SA-4.0 

% This function takes the .csv data from ../data-raw/rim-rings/ directory 
% and prodcues .mat files for use with the other routines in this data set. 
% The data is saved in the ../results/rim-rings/


if ~exist('../results/rim-rings','dir')
    mkdir('../results/rim-rings');
end

num = 1000
option = {'dephasing';'fidelity';'overlap'};
for x = 1:3
    opt = option{x};
	for N = 5:6
		for out = 2:floor(N/2)+1
        	tag = sprintf('../data-raw/rim-rings/%s_dephasing_rim_1s_trace_ring_%d_%d.csv',opt,N,out);
            disp(tag)
        	rim_data = readmatrix(tag);
        	controller_rim = rim_data(51,:)';
        	savetag = sprintf('../results/rim-rings/rim_%s_%d-ring_1-%d.mat',opt,N,out);
        	save(savetag,'controller_rim','rim_data');
        	clear controller_rim rim_data;
		end
	end
end

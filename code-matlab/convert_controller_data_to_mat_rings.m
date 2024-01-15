%% Convert rim controller data to .mat files for rings
%
function [] = convert_controller_data_to_mat_rings
%
% SPDX-FileCopyrightText: Copyright (C) 2023 Sean Patrick O'Neil <seanonei@usc.edu>
% SPDX-FileCopyrightText: Copyright (C) 2024 S M Shermer <lw1660@gmail.com>
% SPDX-License-Identifier: CC-BY-SA-4.0 
%
% This function takes the .csv data for the ring controllers from 
% ../data-raw/controllers-rings/ and prodcues .mat files for use with the
% other routines in this data set.  Controller data is saved in the
% ../results/controllers-rings/ directory 
%

if ~exist('../results/controllers-rings','dir')
    mkdir('../results/controllers-rings');
end

in    = 1;
type  = {'fidelity','dephasing','overlap'};
x = 2;


for N = 5:6
    if N == 5
        target = [2 3];
    else
        target = [2 3 4];
    end
    for out = target	      
        for q = 1:length(type)
            opt = type{q};
	        tag = sprintf('../data-raw/controllers-rings/%s_%d-ring_1-%d.csv',opt,N,out);
			controller_data = readmatrix(tag);
		    sys = convert_controller(controller_data,in,out,opt,N,x);
			savetag = sprintf('../results/controllers-rings/%s_%d-ring_1-%d.mat',opt,N,out);
			save(savetag,'sys');
    	end
	end
end


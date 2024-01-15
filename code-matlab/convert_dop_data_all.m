%% Translate dephasing operator raw data to .mat files

function [] = convert_dephasing_op_all()

% SPDX-FileCopyrightText: Copyright (C) 2023 Sean Patrick O'Neil <seanonei@usc.edu>
% SPDX-FileCopyrightText: Copyright (C) 2024 S M Shermer <lw1660@gmail.com>
% SPDX-License-Identifier: CC-BY-SA-4.0 

% This function takes the .csv data for dephasing operators from the
% ../data-raw/dephasing_op/ directory and prodcues .mat files for use with
% other routines in this data set.  dephasing operators in mat format are
% saved in the ../results/dephasing_op/ directory. 

if ~exist('../results/dephasing_op','dir')
    mkdir('../results/dephasing_op');
end

for N = 5:6
	loadtag = sprintf('../data-raw/dephasing_op/dephasingop_ld_%d_100000.csv',N);
	array = readmatrix(loadtag);	
	DephasingOp = convert_dop(array,N,1);
    savetag=(sprintf('../results/dephasing_op/dephasingop_ld_ring-xx-%d',N));
    save(savetag,'DephasingOp')
end


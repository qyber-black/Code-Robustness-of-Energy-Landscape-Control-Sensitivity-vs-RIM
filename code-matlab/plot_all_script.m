%% script to run all plot routines (plot_all_script.m)

% SPDX-FileCopyrightText: Copyright (C) 2024 Sean Patrick O'Neil <seanonei@usc.edu>
% SPDX-FileCopyrightText: Copyright (C) 2024 S M Shermer <lw1660@gmail.com>
% SPDX-License-Identifier: CC-BY-SA-4.0 

% set noise level for RIM controllers 
noise_choice = 7;
noise = ["0.0" "0.01" "0.02" "0.03" "0.04" "0.05" "0.1"];

disp('Running chain plot routines for log-sensitivity and RIM. . .');
plot_all_chain(noise(noise_choice));
disp('Running chain plot routines for sensitivity and adjusted RIM. . .');
plot_all_chain_v2(noise(noise_choice));
disp('Running ring plot routines for log-sensitivity and RIM. . .');
plot_all_ring;
disp('Running ring plot routines for sensitivity and adjusted RIM. . .');
plot_all_ring_v2;
disp('Running heatmap plot routines for chains. . .');
plot_all_heatmaps_chains(noise(noise_choice))
disp('Running heatmap plot routines for rings. . .');
plot_all_heatmaps_rings
plot_all_rim_sens_heatmap_chains(noise(noise_choice))
plot_all_rim_sens_heatmap_rings
close all;



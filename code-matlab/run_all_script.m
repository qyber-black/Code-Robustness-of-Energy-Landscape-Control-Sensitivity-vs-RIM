%% script to run all conversion and calculation routines (run_all_script.m)

% SPDX-FileCopyrightText: Copyright (C) 2024 Sean Patrick O'Neil <seanonei@usc.edu>
% SPDX-FileCopyrightText: Copyright (C) 2024 S M Shermer <lw1660@gmail.com>
% SPDX-License-Identifier: CC-BY-SA-4.0 

% set noise level for RIM controllers 
noise_choice = 7;
noise = ["0.0" "0.01" "0.02" "0.03" "0.04" "0.05" "0.1"];


disp('convert_dop_all.m');
convert_dop_data_all

disp('convert_controller_data');
convert_rim_data_to_mat_chains
convert_rim_data_to_mat_rings
convert_controller_data_to_mat_chains(noise(noise_choice));
convert_controller_data_to_mat_rings
calc_log_sens_chain(noise(noise_choice))
calc_log_sens_chain_kde(noise(noise_choice))
calc_log_sens_ring
calc_log_sens_ring_kde
disp('sensitivity versus adjusted RIM calculations: ');
calc_sens_vs_adjusted_rim(noise(noise_choice))
analyze_chains_composite(noise(noise_choice))
analyze_chains_composite_v2(noise(noise_choice))
analyze_rings_composite
analyze_rings_composite_v2




# Code - Robustness of Energy Landscape Control - Sensitivity vs RIM

> SPDX-FileCopyrightText: Copyright (C) 2023-2024 Sean Patrick O'Neil <seanonei@usc.edu>\
> SPDX-FileCopyrightText: Copyright (C) 2023 Irtaza Khalid\
> SPDX-FileCopyrightText: Copyright (C) 2023-2024 SM Shermer <lw1660@gmail.com>\
> SPDX-FileCopyrightText: Copyright (C) 2023-2024 Frank C Langbein <frank@langbein.org>
> 
> SPDX-License-Identifier: CC-BY-SA-4.0

Recent achievements in quantum control have resulted in advanced techniques for designing controllers for applications in quantum communication, computing, and sensing. However, the susceptibility of such systems to noise and uncertainties necessitates robust controllers that perform effectively under uncertainty conditions to realize the full potential of quantum devices. The time-domain log-sensitivity and a recently introduced robustness infidelity measure (RIM) are two means to quantify controller robustness in quantum systems. The code and data in this data set allows for comparison between the log-sensitivity and RIM for controllers designed to achieve high-fidelity of excitation transfer in spin rings and chains. This code and data was originally developed for:

[1] S. P. O'Neil, _et al_., "Analyzing and Unifying Robustness Measures for Excitation Transfer Control in Spin Networks," in _IEEE Control Systems Letters_, vol. 7, pp. 1783-1788, 2023, [[DOI:10.1109/LCSYS.2023.3279797]](https://doi.org/10.1109/LCSYS.2023.3279797) [[arXiv:2303.09518]](https://arxiv.org/abs/2303.09518) [[PDF]](https://spinnet.qyber.dev/robustness/paper-cdc2023/cdc2023.pdf).

The results of running the code is available in the `results` submodule that can be cloned via `git submodule update --init --recursive` or directly at [https://qyber.black/spinnet/results-robustness-of-energy-landscape-control-sensitivity-vs-rim/](https://qyber.black/spinnet/results-robustness-of-energy-landscape-control-sensitivity-vs-rim/).

This repository is developed on [qyber/black](https://qyber.black/) at https://qyber.black/spinnet/code-robustness-of-energy-landscape-control-sensitivity-vs-rim and is mirrored on [github](https://github.com/qyber-black/code-robustness-of-energy-landscape-control-sensitivity-vs-rim).

## Versions

**Version 1.0.0**: initial release for the paper

## Notation

The files and figures in this project use the following notation:

```
N   - size (number of spins) for a given ring or chain; takes values of 5 or 6
out - target spin for excitation transfer 
opt - optimizer/optimization condition descriptor
      for chains takes values in {lbgfs,ppo,nmplus,snob}
      for rings takes values in {fidelity,dephasing,overlap}
nl  - (noise level) size of the training noise for RIM data; takes values in [0.0,0.01,0.02,0.03,0.04,0.05,0.1] with default    
      value 0.1 
```

## Data Description

The data intended for input to the computation and analysis code is saved in the `data-raw` directory. Upon downloading the code and data, this directory should be saved directly subordinate to the root directory. The established hierarchy uses the following sub-directories :   
  * `controllers-chains` - All controller data for chain problems is saved in this directory as a `.csv` file with columns `2` through `N+1` containing the bias field values, column `N+2` containing the read-out time, and column `N+3` containing the nominal fidelity.
  * `controllers-rings` - All controller data for ring problems is saved in this directory as a `.csv` file with the first `N` columns containing the bias fields values, column `N+1` containing the read-out time, and column `N+2` containing the nominal fidelity. 
  * `dephasing_op` - Each set of `10^6` dephasing operators for rings and chains of size `N=5` and `N=6` is saved as a `.csv` file in this directory. Each row of the file corresponds to one dephasing operator.
  * `rim-chains` - Spreadsheet of RIM data for chain problems saved as `tn_nl_opt_dephasing_rim_1s_trace_chain_N_out.csv`. 
  * `rim-rings` - Spreadsheet of RIM data for ring problems saved as `opt_dephasing_rim_1s_trace_ring_N_out.csv`.

## Computation and Analysis Code Description

The code to compute and analyze the data based on the inputs form the `data-raw` directory are located in the `code-matlab` directory. The `code-matlab` directory should be saved intact and directly subordiante to the root directory and all routines are designed to be run from the `code-matlab` directory. 

| Routine                                   | Description                                                  |
| ----------------------------------------- | ------------------------------------------------------------ |
| `run_all_script.m`                        | When run, this script calls on each routine below to compute and generate all results saved in the `results` directory. The routines are designed to create each necessary sub-directory within `results`. |
| `convert_rnim_raw_data_to_mat_chains.m`    | Converts relevant RIM data from the `.csv` files in `data-raw/rim-chains` and saves as `.mat` files in `results/rim-chains` as `rim_noise_nl_opt_N-chain_1-out.mat`. |
| `convert_rim_raw_data_to_mat_rings.m`     | Converts relevant RIM data from the `.csv` files in `data-raw/rim-rings` and saves as `.mat` files in `results/rim-rings as `rim_opt_N-ring_1-out.mat`. |
| `convert_dop_data_all.m`                  | Takes the `.csv` files for the dephasing processes from `data-raw/dephasing_op` and translates them into `.mat` files saved in `results/dephasing_op` as `dephasingop_1d_ring-xx-N.mat`. Data saved as a cell array of 100,000 upper triangular matrices. Calls on `convert_dop.m` and `InvertDephasing.m`. |
| `convert_dop.m`                           | Executes conversion of dephasing operators in `.csv` spreadsheet to an array for use with `MATLAB` routines. |
| `InvertDephasing.m`                       | Tests each dephasing operator for complete positivity.       |
| `convert_controller_data_to_mat_chains.m` | Converts the controller data in `data-raw/controllers-chains` to usable `.mat` files and saves as `noise_nl_opt_N-chain_1-out.mat` in the `results/controllers-chains` directory. Calls on `convert_controller.m`. |
| `convert_controller_data_to_mat_rings.m`  | Converts the controller data in `data-raw/controllers-rings` to usable `.mat` files and saves as `opt_N-ring_1-out.mat` in the `results/controllers-rings` directory. Calls on `convert_controller.m`. |
| `convert_controller.m`                    | Reads in `.csv` file from `data-raw/controllers-{rings,chains}` and converts data to a `MATLAB` structure. |
| `calc_log_sens_chain.m`                   | Computes the analytic log-sensitivity for chain problems and saves as `log_sens_nl_opt_N-chain_1-out.mat` in the `results/log_sens-chains` directory. Calls on `bloch_basis.m`. |
| `calc_log_sens_chain_kde.m`               | Computes the KDE-based log-sensitivity for chain problems and saves as `log_kde_nl_opt_N-chain_1-out.mat` in the `results/kde-chains` directory. Calls on `AnalyzeRobustnessDephasing.m` and `AnalyzeErrorDensity.m`. |
| `calc_log_sens_ring.m`                    | Computes the analytic log-sensitivity for ring problems and saves as `log_sens_opt_N-ring_1-out.mat` in the `results/log_sens-rings` directory . Calls on `bloch_basis.m`. |
| `calc_log_sens_ring_kde.m`                | Computes the KDE-based log-sensitivity for ring problems and saves as `log_kde_opt_N-ring_1-out.mat` in the `results/kde-rings` directory. Calls on `AnalyzeRobustnessDephasing.m` and `AnalyzeErrorDensity.m`. |
| `bloch_basis.m`                           | Computes basis matrices for computation of the analytic log-sensitivity. |
| `AnalyzeRobustnessDephasing.m`            | Produces kernel of error values for input to KDE.            |
| `AnalyzeErrorDensity.m`                   | Produces kernel density estimate of mean fidelity error.     |
| `calc_sens_vs_adjusted_rim.m`             | Computes the relative difference between the analytic differential sensitivity and the adjusted RIM for all ring and chain controllers and saves the result as `results/sens_v_rim_sens.mat`. |
| `analyze_chains_composite.m`              | Computes correlation between analytic and KDE-based log-sensitivity versus RIM as well as the correlation between both log-sensitivity measures and RIM versus the nominal fidelity error and saves the results as `results/corr-chain-noise_nl-log_sens-rim.csv` and `results/corr-chain-noise_nl-log_sens-rim-err.csv`. |
| `analyze_chains_composite_v2.m`           | Computes correlation between both differential sensitivity measures (analytic and KDE) versus RIM as well as the correlation between both sensitivity measures and RIM versus the nominal fidelity error and saves results as `results/corr-chain-noise_nl_-sens-adj_rim.csv` and `results/corr-chain-noise_nl-sens-adj_rim-err.csv` |
| `analyze_rings_composite.m`               | Computes correlation between analytic and KDE-based log-sensitivity versus RIM as well as the correlation between both log-sensitivity measures and RIM versus the nominal fidelity error and saves the results as `results/corr-ring-log_sens_rim.csv` and `results/corr-ring-log_sens-rim-err.csv`. |
| `analyze_rings_composite_v2.m`            | Computes correlation between both differential sensitivity measures (analytic and KDE) versus RIM as well as the correlation between both sensitivity measures and RIM versus the nominal fidelity error and saves results as `results/corr-ring-sens_adj_rim.csv` and `results/corr-ring-sens-adj_rim-err.csv`. |

The python files to process the RIM results are in `code-python` and the dependencies can be installed via `pip3 install -r code-python/requirements.txt`. The files are:

| File | Description |
| ---- | ----------- |
| `rim_vs_logsensitivity.py` | RIM vs. log-sensitivity for chains and rings; creats the `results/figures/rim_vs_logsensitivity-{rings,chains}` plots. |
| `generate_rims_heatmap.py` | RIM heatmaps and correlations with log sensitivity; computes the `rim_heat-{chains,rings}` files in `results/results` and corresponding plots in `results/figures/rim_heat-{chains,rings}`. |
| `plot_all_adj_rim_vs_kdtau.py` | Plot Kendall taus for RIM vs. differential sensitivity. |
| `rim_vs_diffsens_dists_chains.py` | Kendall-tau between RIM and log sensitivity for chains (for `plot_all_adj_rim_vs_kdtau.py`). |
| `rim_vs_diffsens_dists_rings.py` | Kendall-tau between RIM and log sensitivity for rings (for `plot_all_adj_rim_vs_kdtau.py`). |
| `dephasing_fid.py` | Compute infidelities with dephasing operators. |
| `noise_model.py` |Perturbation models. |
| `utilities.py` | Helper functions. |

## Plot Routines

The following plot routines are designed to be run from the `code-matlab` directory following completion of the computation and analysis routines above. 

| Routine                              | Description                                                                                                                                                                                                                                                                                                                                                                        |
| ------------------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `plot_all_script.m`                  | This script executes all plot routines below and requires the user to input the `noise_choice` variable to select the proper RIM data.                                                                                                                                                                                                                                             |
| `plot_all_chain.m`                   | Constructs plots of both log-sensitivity measures, RIM, and nominal fidelity error versus controller index saved as `figures/composite-chains/comparison_nl_opt_N-chain_1-out.{fig.png}` and scatter plots of both log-sensitivity measures and RIM versus nominal fidelity error saved as `figures/composite-chains/scatter_nl_opt_N-chain_1-out.{fig,png}`.                      |
| `plot_all_chain_v2.m`                | Constructs plots of both sensitivity measures, adjusted RIM, and nominal fidelity error versus controller index saved as `figures/composite-chains-v2/comparison_v2_nl_opt_N-chain_1-out.{fig.png}` and scatter plots of both sensitivity measures and adjusted RIM versus nominal fidelity error saved as `figures/composite-chains-v2/scatter_v2_nl_opt_N-chain_1-out.{fig,png}`. |
| `plot_all_ring.m`                    | Constructs plots of both log-sensitivity measures, RIM, and nominal fidelity error versus controller index saved as `figures/composite-rings/comparison_opt_N-ring_1-out.{fig.png}` and scatter plots of log-sensitivity measures and RIM versus nominal fidelity error saved as `figures/composite-rings/scatter_opt_N-ring_1-out.{fig,png}`.                                     |
| `plot_all_ring_v2.m`                 | Constructs plots of both sensitivity measures, adjusted RIM, and nominal fidelity error versus controller index saved as `figures/composite-rings-v2/comparison_v2_opt_N-ring_1-out.{fig.png}` and scatter plot of both sensitivity measures and adjusted RIM versus nominal fidelity error saved as `figures/composite-rings-v2/scatter_v2_opt_N-ring_1-out.{fig,png}`.            |
| `plot_all_heatmaps_chains.m`         | Produces a multi-pane figure with a heatmap of the fidelity error distribution, bias fields, initial density matrix, and density matrix at read-out time and saves as `figures/heatmaps-chains/noise_nl_opt_N-chain_1-out_ctrl_#.{fig.png}`. User may set number of controllers to plot in the routine under `index` variable. Calls on `PlotHeatMaps.m`.                          |
| `plot_all_heatmaps_rings.m`          | Produces a multi-pane figure with a heatmap of the fidelity error distribution, bias fields, initial density matrix, and density matrix at read-out time and saves as `figures/heatmaps-rings/opt_N-ring_1-out_ctrl_#.{fig.png}`. User may set number of controllers to plot in the routine under `index` variable. Calls on `PlotHeatMaps.m`.                                     |
| `plot_all_rim_sens_heatmap_chains.m` | Produces a heatmap of RIM for increasing perturbation strength versus the differential sensitivity and saves as `figures/heat_maps/heat_map_chain_opt_N_out.{fig,png}`.                                                                                                                                                                                                            |
| `plot_all_rim_sens_heatmap_rings.m`  | Produces a heatmap of RIM for increasing perturbation strength versus the differential sensitivity and saves as `figures/heat_maps/heat_map_ring_opt_N_out.{fig,png}`.                                                                                                                                                                                                             |
| `PlotHeatMaps.m`                     | Routine that produces the heat map for `plot_all_heatmaps_{ring/chains}.m`.                                                                                                                                                                                                                                                                                                        |
| `heatmaptext.m`                      | Helper script for producing heat maps called on by `PlotHeatMaps.m`.                                                                                                                                                                                                                                                                                                               |
| `Inferno.m`                          | Helper script for coloring heat maps called on by `PlotHeatMaps.m`.                                                                                                                                                                                                                                                                                    |

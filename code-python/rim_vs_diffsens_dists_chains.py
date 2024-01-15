#!/usr/bin/env python3
#
# SPDX-FileCopyrightText: Copyright (C) 2023 Irtaza Khalid
# SPDX-FileCopyrightText: Copyright (C) 2024 Frank C Langbein <frank@langbein.org>
# SPDX-License-Identifier: CC-BY-SA-4.0

if __name__ == "__main__":
    import os
    from scipy.io import loadmat
    import numpy as np
    from scipy.stats import kendalltau
    from utilities import get_ranks
    import pickle

    CHAIN_PATH_logs = os.path.join("results", "results", "log_sens-chains")
    CHAIN_PATH_rim = os.path.join("data-raw", "rim-chains")
    kendall_tau_dir = os.path.join("results", "results", "kendall_tau_rims_vs_diffsens")
    corr_fname = "chains.pickle"

    os.makedirs(kendall_tau_dir, exist_ok=True)
    kendall_taus_0p05 = []
    kendall_taus_0p001 = []
    kendall_taus_0p001_p_errs = []
    kendall_taus_0p05_logsens = []
    kendall_taus_0p001_logsens = []
    kendall_taus_0p05_p_err = []
    kendall_taus_0p05_p_err = []
    algos = ["snob", "nmplus", "ppo"]
    nspins = [5, 6]
    outspins = {5: [3, 5], 6: [4, 6]}
    for algo in algos:
        for nspin in nspins:
            for outspin in outspins[nspin]:
                # file paths
                rim_path = f"tn_0.1_{algo}_dephasing_rim_1s_trace_chain_{nspin}_{outspin}.csv"
                log_sens_path = f"log_sens_0.1_{algo}_{nspin}-chain_1-{outspin}.mat"
                # load data
                # analytical log sensitivity
                load_dict = loadmat(os.path.join(CHAIN_PATH_logs, log_sens_path))
                rims = np.loadtxt(os.path.join(CHAIN_PATH_rim, rim_path), delimiter=",")
                # preprocess
                mean_log_sens = load_dict["log_sens"].mean(axis=-1)
                rims_sorted = rims  # sort_2d_array_by_zero_index_array(rims)
                infids_mat = load_dict["err"].reshape(-1)
                errs = rims_sorted[0]
                # differential sensitivity
                diffsens = errs * mean_log_sens
                # adjusted rim
                adjusted_rims_0p05 = rims_sorted[50] - errs
                adjusted_rims_0p001 = rims_sorted[1] - errs
                # get ranks
                adjusted_rims_0p05_ranks = get_ranks(
                    adjusted_rims_0p05,
                    cluster=False,
                )[0][:]
                adjusted_rims_0p001_ranks = get_ranks(
                    adjusted_rims_0p001, cluster=False
                )[0][:]
                adjusted_rims_0p001_p_err_ranks = get_ranks(
                    adjusted_rims_0p001 + errs, cluster=False
                )[0][:]
                adjusted_rims_0p05_p_err_ranks = get_ranks(
                    adjusted_rims_0p05 + errs, cluster=False
                )[0][:]
                rims_0p001_ranks = get_ranks(rims_sorted[1], cluster=False)[0][:]
                rims_0p05_ranks = get_ranks(rims_sorted[50], cluster=False)[0][:]
                diffsens_ranks = get_ranks(diffsens, cluster=False)[0][:]
                logs_ranks = get_ranks(mean_log_sens, cluster=False)[0][:]
                # store various kendall taus
                kendall_taus_0p05.append(
                    kendalltau(adjusted_rims_0p05, diffsens_ranks)[0]
                )
                kendall_taus_0p001.append(
                    kendalltau(adjusted_rims_0p001, diffsens_ranks)[0]
                )
                kendall_taus_0p001_p_errs.append(
                    kendalltau(adjusted_rims_0p001_p_err_ranks, diffsens_ranks)[0]
                )
                kendall_taus_0p05_logsens.append(
                    kendalltau(rims_0p05_ranks, logs_ranks)[0]
                )
                kendall_taus_0p001_logsens.append(
                    kendalltau(rims_0p001_ranks, logs_ranks)[0]
                )
                kendall_taus_0p05_p_err.append(
                    kendalltau(adjusted_rims_0p05_p_err_ranks, diffsens_ranks)[0]
                )
    pickle.dump(
        [
            kendall_taus_0p05,
            kendall_taus_0p001,
            kendall_taus_0p001_p_errs,
            kendall_taus_0p05_p_err,
            kendall_taus_0p05_logsens,
            kendall_taus_0p001_logsens,
        ],
        open(os.path.join(kendall_tau_dir, corr_fname), "wb"),
    )

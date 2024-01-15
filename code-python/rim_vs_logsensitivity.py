#!/usr/bin/env python3
#
# SPDX-FileCopyrightText: Copyright (C) 2023 Irtaza Khalid
# SPDX-FileCopyrightText: Copyright (C) 2024 Frank C Langbein <frank@langbein.org>
# SPDX-License-Identifier: CC-BY-SA-4.0

from utilities import get_ranks, sort_2d_array_by_zero_index_array
import os
from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import kendalltau


def get_kde_based_log_sens(file):
    load_dict = loadmat(file)
    kde_dens = load_dict["density"]
    kde_dens_clean = {"mean": [], "sensitivity": [], "err_kde": []}
    for controller in range(len(kde_dens[0])):
        kde_dens_clean["mean"].append(kde_dens[0][controller][0][0][0])
        kde_dens_clean["sensitivity"].append(kde_dens[0][controller][0][0][8])
        kde_dens_clean["err_kde"].append(kde_dens[0][controller][0][0][9])
    return np.array(kde_dens_clean["sensitivity"])[:, 0].squeeze(-1)


def do_chains():
    CHAIN_PATH_logs = os.path.join("results", "results", "log_sens-chains")
    CHAIN_PATH_kde_logs = os.path.join("results", "results", "kde-chains")
    CHAIN_PATH_rim = os.path.join("data-raw", "rim-chains")
    plots_dir = os.path.join("results", "figures", "rim_vs_logsensitivity-chains")
    if not os.path.exists(plots_dir):
        os.makedirs(plots_dir, exist_ok=True)
    algos = ["snob", "nmplus", "ppo"]
    nspins = [5, 6]
    outspins = {5: [3, 5], 6: [4, 6]}
    for algo in algos:
        for nspin in nspins:
            for outspin in outspins[nspin]:
                # file paths
                rim_path = f"tn_0.1_{algo}_dephasing_rim_1s_trace_chain_{nspin}_{outspin}.csv"
                log_sens_path = f"log_sens_0.1_{algo}_{nspin}-chain_1-{outspin}.mat"
                kde_log_sens_path = f"log_kde_0.1_{algo}_{nspin}-chain_1-{outspin}.mat"
                # load data
                # analytical log sensitivity
                load_dict = loadmat(os.path.join(CHAIN_PATH_logs, log_sens_path))
                # kde based log sensitivity
                kde_log_sens = get_kde_based_log_sens(os.path.join(CHAIN_PATH_kde_logs, kde_log_sens_path))
                rims = np.loadtxt(os.path.join(CHAIN_PATH_rim, rim_path), delimiter=",")
                # preprocess
                mean_log_sens = load_dict["log_sens"].mean(axis=-1)
                if not np.allclose(mean_log_sens - kde_log_sens, 0, atol=1e-4):
                    err = np.abs(mean_log_sens - kde_log_sens).max()
                    print(
                        "WARNING: kde log sens and analytical log sens average mismatch by: ",
                        err,
                    )
                rims_sorted = rims  # sort_2d_array_by_zero_index_array(rims)
                assert np.allclose(
                    1 - rims_sorted[0] - load_dict["fid_rim"].reshape(-1), 0
                ), "rim and log sens rank order doesn't match!"

                # plot measures
                fig, ax = plt.subplots()
                ax2 = ax.twinx()
                ax.plot(
                    range(len(mean_log_sens)),
                    mean_log_sens,
                    label="analytic mean log sens",
                    c="blue",
                )
                ax.plot(
                    range(len(mean_log_sens)),
                    kde_log_sens,
                    label="kde mean log sens",
                    c="green",
                )
                ax2.plot(
                    range(len(rims_sorted[50])),
                    rims_sorted[-1],
                    label="rim at 0.05",
                    c="red",
                )
                ax.set_xlabel("controller index")
                ax.set_ylabel("Log sensitivity")
                ax2.set_ylabel("RIM at dephasing strength 0.05")
                fig.savefig(os.path.join(plots_dir, f"tn_0.1_{algo}_{nspin}-{outspin}.png"), dpi=600)

                # plot measure in ranks
                logs_ranks = get_ranks(
                    mean_log_sens,
                    cluster=False,
                )[
                    0
                ][:]
                kde_logs_ranks = get_ranks(
                    kde_log_sens,
                    cluster=False,
                )[
                    0
                ][:]
                rim_ranks = get_ranks(rims_sorted[50], cluster=False)[0][:]
                fig2, ax3 = plt.subplots()
                ax3.plot(
                    range(len(logs_ranks)),
                    logs_ranks,
                    label="rank analytic mean log sens",
                )
                ax3.plot(
                    range(len(logs_ranks)),
                    kde_logs_ranks,
                    label="rank kde mean log sens",
                )
                ax3.plot(range(len(logs_ranks)), rim_ranks, label="rank RIM at 0.05")
                ax3.set_title(
                    "kendall tau: "
                    + str(np.round(kendalltau(logs_ranks, rim_ranks)[0], 3))
                )
                ax3.set_xlabel("controller index")
                ax3.set_ylabel("rank w.r.t. measure")
                ax3.legend()
                fig2.savefig(os.path.join(plots_dir, f"tn_0.1_{algo}_{nspin}-{outspin}_ranks.png"), dpi=600)

                # correlation as a function dephasing strength
                get_taus_at_rim = lambda rims_sorted_i: kendalltau(
                    logs_ranks, get_ranks(rims_sorted_i, cluster=False)[0]
                )[0]
                get_taus_at_rim_kde = lambda rims_sorted_i: kendalltau(
                    kde_logs_ranks, get_ranks(rims_sorted_i, cluster=False)[0]
                )[0]
                taus_vs_strength = list(map(get_taus_at_rim, rims_sorted))
                taus_vs_strength_kde = list(map(get_taus_at_rim_kde, rims_sorted))
                fig3, ax4 = plt.subplots()
                ax4.plot(
                    np.linspace(0, 0.1, 101),
                    taus_vs_strength,
                    label="analytic mean log sensitivity",
                    alpha=0.5,
                )
                ax4.plot(
                    np.linspace(0, 0.1, 101),
                    taus_vs_strength_kde,
                    label="kde mean log sensitivity",
                    alpha=0.5,
                )
                ax4.set_xlabel("dephasing strength")
                ax4.set_ylabel(
                    "Correlation between analytic mean log sensitivity and dephasing RIM"
                )
                fig3.savefig(os.path.join(plots_dir, f"tn_0.1_{algo}_{nspin}-{outspin}_taus_vs_strength.png"), dpi=600)

                # middle of the pack correlation as a function dephasing strength
                logs_ranks = logs_ranks[40:80]
                kde_logs_ranks = kde_logs_ranks[40:80]
                get_taus_at_rim = lambda rims_sorted_i: kendalltau(
                    logs_ranks, get_ranks(rims_sorted_i, cluster=False)[0][40:80]
                )[0]
                get_taus_at_rim_kde = lambda rims_sorted_i: kendalltau(
                    kde_logs_ranks, get_ranks(rims_sorted_i, cluster=False)[0][40:80]
                )[0]
                taus_vs_strength = list(map(get_taus_at_rim, rims_sorted))
                taus_vs_strength_kde = list(map(get_taus_at_rim_kde, rims_sorted))
                fig4, ax5 = plt.subplots()
                ax5.plot(
                    np.linspace(0, 0.1, 101),
                    taus_vs_strength,
                    label="analytic mean log sensitivity",
                    alpha=0.5,
                )
                ax5.plot(
                    np.linspace(0, 0.1, 101),
                    taus_vs_strength_kde,
                    label="kde mean log sensitivity",
                    alpha=0.4,
                )
                ax5.set_xlabel("dephasing strength")
                ax5.set_ylabel(
                    "Correlation between analytic mean log sensitivity and dephasing RIM"
                )
                fig4.savefig(os.path.join(plots_dir, f"tn_0.1_{algo}_{nspin}-{outspin}_taus_vs_strength_middle.png"), dpi=600)

                plt.close("all")

def do_rings():
    RING_LOGS_PATH = os.path.join("results", "results", "log_sens-rings")
    CHAIN_KDE_LOGS_PATH = os.path.join("results", "results", "kde-rings")
    RING_rim_path = os.path.join("data-raw", "rim-rings")
    plots_dir = os.path.join("results", "figures", "rim_vs_logsensitivity-rings")
    if not os.path.exists(plots_dir):
        os.makedirs(plots_dir, exist_ok=True)
    opts = ["fidelity", "dephasing", "overlap"]
    nspins = [5, 6]
    outspins = {5: [2, 3], 6: [2, 3, 4]}
    for opt in opts:
        for nspin in nspins:
            for outspin in outspins[nspin]:
                # file paths
                rim_path = f"{opt}_dephasing_rim_1s_trace_ring_{nspin}_{outspin}.csv"
                log_sens_path = f"log_sens_{opt}_{nspin}-ring_1-{outspin}.mat"

                # load data
                # analytical log sensitivity
                load_dict = loadmat(os.path.join(RING_LOGS_PATH, log_sens_path))

                rims = np.loadtxt(os.path.join(RING_rim_path, rim_path), delimiter=",")
                # preprocess
                mean_log_sens = load_dict["log_sens"].mean(axis=-1)
                rims_sorted = rims  # sort_2d_array_by_zero_index_array(rims)
                infids_mat = load_dict["err"].reshape(-1)
                assert np.allclose(
                    rims_sorted[0] - infids_mat, 0
                ), "rim and log sens rank order doesn't match!"

                # plot measures
                fig, ax = plt.subplots()
                ax2 = ax.twinx()
                ax.plot(
                    range(len(mean_log_sens)),
                    mean_log_sens,
                    label="analytic mean log sens",
                    c="blue",
                )

                ax2.plot(
                    range(len(rims_sorted[50])),
                    rims_sorted[-1],
                    label="rim at 0.05",
                    c="red",
                )
                ax.set_xlabel("controller index")
                ax.set_ylabel("Log sensitivity")
                ax2.set_ylabel("RIM at dephasing strength 0.05")
                fig.savefig(os.path.join(plots_dir, f"ring_{opt}_{nspin}-{outspin}.png"), dpi=600)

                # plot measure in ranks
                logs_ranks = get_ranks(
                    mean_log_sens,
                    cluster=False,
                )[
                    0
                ][:]
                rim_ranks = get_ranks(rims_sorted[50], cluster=False)[0][:]
                fig2, ax3 = plt.subplots()
                ax3.plot(
                    range(len(logs_ranks)),
                    logs_ranks,
                    label="rank analytic mean log sens",
                )
                ax3.plot(range(len(logs_ranks)), rim_ranks, label="rank RIM at 0.05")
                ax3.set_title(
                    "kendall tau: "
                    + str(np.round(kendalltau(logs_ranks, rim_ranks)[0], 3))
                )
                ax3.set_xlabel("controller index")
                ax3.set_ylabel("rank w.r.t. measure")
                ax3.legend()
                fig2.savefig(os.path.join(plots_dir, f"ring_{opt}_{nspin}-{outspin}_ranks.png"), dpi=600)

                # correlation as a function dephasing strength
                get_taus_at_rim = lambda rims_sorted_i: kendalltau(
                    logs_ranks, get_ranks(rims_sorted_i, cluster=False)[0]
                )[0]
                taus_vs_strength = list(map(get_taus_at_rim, rims_sorted))

                fig3, ax4 = plt.subplots()
                ax4.plot(
                    np.linspace(0, 0.1, 101),
                    taus_vs_strength,
                    label="analytic mean log sensitivity",
                    alpha=0.5,
                )

                ax4.set_xlabel("dephasing strength")
                ax4.set_ylabel(
                    "Correlation between analytic mean log sensitivity and dephasing RIM"
                )
                fig3.savefig(os.path.join(plots_dir, f"ring_{opt}_{nspin}-{outspin}_taus_vs_strength.png"), dpi=600)

                # middle of the pack correlation as a function dephasing strength
                logs_ranks = logs_ranks[40:80]
                get_taus_at_rim = lambda rims_sorted_i: kendalltau(
                    logs_ranks, get_ranks(rims_sorted_i, cluster=False)[0][40:80]
                )[0]
                taus_vs_strength = list(map(get_taus_at_rim, rims_sorted))

                fig4, ax5 = plt.subplots()
                ax5.plot(
                    np.linspace(0, 0.1, 101),
                    taus_vs_strength,
                    label="analytic mean log sensitivity",
                    alpha=0.5,
                )

                ax5.set_xlabel("dephasing strength")
                ax5.set_ylabel(
                    "Correlation between analytic mean log sensitivity and dephasing RIM"
                )
                fig4.savefig(os.path.join(plots_dir, f"ring_{opt}_{nspin}-{outspin}_taus_vs_strength_middle.png"), dpi=600)

                plt.close("all")

if __name__ == "__main__":
    do_chains()
    do_rings()

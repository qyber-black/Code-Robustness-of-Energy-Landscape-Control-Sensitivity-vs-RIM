#!/usr/bin/env python3
#
# SPDX-FileCopyrightText: Copyright (C) 2023 Irtaza Khalid
# SPDX-FileCopyrightText: Copyright (C) 2024 Frank C Langbein <frank@langbein.org>
# SPDX-License-Identifier: CC-BY-SA-4.0

import matplotlib.pyplot as plt
import numpy as np
import pickle
import os

kendall_tau_dir = os.path.join("results", "results", "kendall_tau_rims_vs_diffsens")
plots_dir = os.path.join("results", "figures", "kendall_tau_rims_vs_diffsens")
os.makedirs(plots_dir, exist_ok=True)
corr_fnames = ["rings", "chains"]
alpha = 0.5
linewidth = 1
fig, ax = plt.subplots(nrows=2, figsize=(4, 7))
for i, corr_fname in enumerate(corr_fnames):
    (
        kendall_taus_0p05,
        kendall_taus_0p001,
        kendall_taus_0p001_p_errs,
        kendall_taus_0p05_p_err,
        kendall_taus_0p05_logsens,
        kendall_taus_0p001_logsens,
    ) = pickle.load(open(os.path.join(kendall_tau_dir, corr_fname + ".pickle"), "rb"))

    bins = 10 ** np.linspace(-10, 1, 20)
    # plt.figure()
    bins = 80
    fontsize = 13
    ax[i].hist(
        kendall_taus_0p05,
        bins=bins,
        range=(-1.0, 1),
        alpha=alpha,
        edgecolor="black",
        linewidth=linewidth,
        label=r"$\widetilde{{\rm RIM}}_1(0.05)$"
        + " vs. "
        + r"$E_{\mathbf{P}(\mathbf{S})}\left[\zeta\right]$",
    )
    # ax[i].hist(1-np.array(kendall_taus_0p05), bins=bins, alpha=.5, label="adjusted rim(0.05) vs diffsens tau")
    ax[i].hist(
        kendall_taus_0p001,
        bins=bins,
        range=(-1.0, 1),
        alpha=alpha,
        edgecolor="black",
        linewidth=linewidth,
        label=r"$\widetilde{{\rm RIM}}_1(0.001)$"
        + " vs. "
        + r"$E_{\mathbf{P}(\mathbf{S})}\left[\zeta\right]$",
    )
    ax[i].hist(
        kendall_taus_0p001_p_errs,
        bins=bins,
        range=(-1.0, 1),
        alpha=alpha,
        edgecolor="black",
        linewidth=linewidth,
        label=r"${\rm RIM}_1(0.001)$"
        + " vs. "
        + r"$E_{\mathbf{P}(\mathbf{S})}\left[\zeta\right]$",
    )
    ax[i].hist(
        kendall_taus_0p05_p_err,
        bins=bins,
        range=(-1.0, 1),
        edgecolor="black",
        linewidth=linewidth,
        alpha=alpha,
        label=r"${\rm RIM}_1(0.05)$"
        + " vs. "
        + r"$E_{\mathbf{P}(\mathbf{S})}\left[\zeta\right]$",
    )
    ax[i].hist(
        kendall_taus_0p001_logsens,
        bins=bins,
        range=(-1.0, 1),
        edgecolor="black",
        linewidth=linewidth,
        alpha=alpha,
        label=r"$\widetilde{{\rm RIM}}_1(0.001)$"
        + " vs. "
        + r"$E_{\mathbf{P}(\mathbf{S})}\left[s\right]$",
    )
    # ax[i].hist(kendall_taus_0p05_logsens, bins=bins, range=(-1., 1), alpha=alpha,
    #         edgecolor="black",
    #     linewidth=linewidth,
    # label=r"${\rm RIM}_1(0.05)$"+" vs. "+r"$E_{\mathbf{P}(\mathbf{S})}\left[s\right]$")
    # ax[i].hist(1-np.array(kendall_taus_0p001), bins=bins, alpha=.5, label="adjusted rim(0.001) vs diffsens tau")
    # ax[i].hist(kendall_taus_0p001, bins=20, range=(0.5, 1), alpha=.5, label="rim(0.001) vs diffsens tau")
    # ax[i].xscale('log')
    if i != 0:
        ax[i].set_xlabel(r"$\tau$", fontsize=fontsize)
    ax[i].set_ylabel("count", fontsize=fontsize)
    ax[i].set_title(f"{corr_fname}", fontsize=fontsize)
    if i == 0:
        ax[i].legend(fontsize=fontsize - 3)
    # plt.show()
fig.savefig(os.path.join(plots_dir, f"kendalltau_dist_rings_p_chains.png"), dpi=600)

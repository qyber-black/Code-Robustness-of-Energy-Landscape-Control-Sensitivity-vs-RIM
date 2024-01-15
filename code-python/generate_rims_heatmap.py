#!/usr/bin/env python3
#
# SPDX-FileCopyrightText: Copyright (C) 2023 Irtaza Khalid
# SPDX-FileCopyrightText: Copyright (C) 2024 Frank C Langbein <frank@langbein.org>
# SPDX-License-Identifier: CC-BY-SA-4.0

import os
from noise_model import Dephasing_with_perturbation, unstructured_perturbation
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool
from typing import Callable, List, Any
import seaborn as sns
from utilities import get_top_k_by_fid, sort_2d_array_by_zero_index_array, get_pairwise_kendall_tau_ranks, pcolorwrm

def make_dephasing_rim_1s(
    Nspin=5,
    outspin=1,
    topo: str = "ring",
    controller_path: str = None,
    ham_noisy=False,
    num_dephasing_ops=1000,
    opt="fidelity",
    DATA_PATH=None,
    FIG_PATH=None,
    dep_strength=0.1,
    dephasing_res=101,
    plot=False,
    skip_first_col=False,
    thanos=False,
    fixed_noise=False,
) -> None:
    """
    Make dephasing RIMs for a given set of arbitrary time-independent spin chain
    controllers.

    Parameters
    ----------
    Nspin : int, optional
        Length of the chain, by default 5
    outspin : int, optional
        Final transition spin atom (by default all chains are initialized
        in the first atom as a convention), by default 1
    topo : str, optional
        chain topology can be either "chain" or "ring", by default "ring"
    controller_path : str, optional
        file location of the controllers (assumed csv), by default None
    ham_noisy : bool, optional
        Add structured perturbations, by default False
    num_dephasing_ops : int, optional
        Number of dephasing jump operators, by default 1000
    opt : str, optional
        prefix to make the controllers more unique, by default "fidelity"
    DATA_PATH : str, optional
        roughly a global location of all the `controller paths`,
        by default None
    FIG_PATH : str, optional
        roughly a global location of where the figures are stored,
        by default None
    dep_strength : float, optional
        deterministic scalar \delta, the strength of the dephasing,
        by default 0.1
    dephasing_res : int, optional
        number of incremental steps to to up to `dep_strength` starting from 0,
        by default 101
    plot : bool, optional
        generate plots again?, by default False
    skip_first_col : bool, optional
        remove the first column in the controllers (could be an index as an artifact of
        the initial storing mechanism), by default False
    fixed_noise : float, optional
        if not None, then use this fixed noise strength instead of the dep_strength
    """
    # load the dephasing operators and the noise model
    # conform to the 1 index notation with the rest of the codebase so that its easy to read
    # and write about later
    
    if fixed_noise:
        RIM_save_path_1 += f"_fixed_noise_{fixed_noise}"
    else:
        RIM_save_path_1 = f"{opt}_dephasing_rim_1s_trace_{topo}_{Nspin}_{outspin+1}"
    if not os.path.exists(DATA_PATH + RIM_save_path_1 + ".csv"):
        nmodel_dep = Dephasing_with_perturbation(
            Nspin=Nspin,
            outspin=outspin,
            num_dephasing_ops=num_dephasing_ops,
            dephasing_op_fname=os.path.join("data-raw", "dephasing_op", f"dephasingop_ld_{Nspin}_100000.csv"),
            topo=topo,
            dephasing_res=dephasing_res,
            dep_strength=dep_strength,
        )
        # load the controllers
        if topo == "chain":
            controllers = pd.read_csv(controller_path).to_numpy()
        else:
            controllers = pd.read_csv(controller_path, header=None if not thanos else 0).to_numpy()  # don't ignore the first row with header=None

        if skip_first_col:
            controllers = controllers[:, 1:]
        if not fixed_noise:
            # evaluate the fidelity of the controllers
            deph_fid_map_all = lambda controller: nmodel_dep.evaluate_noisy_fidelity(controller, ham_noisy=ham_noisy)
        else:
            deph_fid_map_all = lambda controller: nmodel_dep._evaluate_noisy_fidelity_at_fixed_dep_strength(controller,
                                                                                                            ham_noisy=ham_noisy,
                                                                                                            dep_strength=fixed_noise)
        # shape (num_controllers, dephasing_ops, 2, dephasing_steps)
        all_fids = np.array(list(map(deph_fid_map_all, controllers)))
        num_controllers, *_ = all_fids.shape
        RIMs = 1 - all_fids.mean(axis=1)
        if not fixed_noise:
            RIMs_1, RIMs_2 = (RIMs[:, 0, :], RIMs[:, 1, :])  # two definitions of the fidelities (trace, norm)
        else:
            RIMs_1 = RIMs
        # get top noiseless fidelity -- dummy for now
        RIMs_1_topk = get_top_k_by_fid(RIMs_1.T, topk=num_controllers)
        # sort the RIMs by the topk noiseless fidelities
        RIMs_1_sorted = sort_2d_array_by_zero_index_array(RIMs_1_topk)

        np.savetxt(os.path.join(DATA_PATH, RIM_save_path_1 + ".csv"), RIMs_1_sorted, delimiter=",")

    else:
        # save RIMs sorted by fidelity
        RIMs_1_sorted = np.loadtxt(os.path.join(DATA_PATH, RIM_save_path_1 + ".csv"), delimiter=",")

    if plot:
        # plot the RIMs
        fig, ax = plt.subplots(figsize=(10, 5))
        pcolorwrm(
            RIMs_1_sorted,
            f"{topo} N={Nspin}, O={outspin+1}",
            fig,
            ax,
            pltcolbar=True,
            sigma_sims=np.linspace(0, dep_strength, dephasing_res),
        )
        ax.set_xlabel("controller index", fontsize=20)
        ax.set_ylabel(r"$\delta$", fontsize=20)
        fig.savefig(os.path.join(FIG_PATH, RIM_save_path_1 + ".png"), dpi=600, bbox_inches="tight")
        plt.close(fig)
        # kendall tau heatmaps
        fig2, ax2 = plt.subplots(figsize=(5, 5))
        kendalltaus = get_pairwise_kendall_tau_ranks(RIMs_1_sorted)
        sns.heatmap(kendalltaus, ax=ax2, cbar_kws={'label': r'$\tau(\delta_1, \delta_2)$'})
        ax2.set_yticks(np.linspace(0, 100, 11))
        ax2.set_yticklabels(np.linspace(0, 0.1, 11))
        ax2.set_xticks(np.linspace(0, 100, 11))
        ax2.set_xticklabels(np.linspace(0, 0.1, 11))
        ax2.set_xlabel(r"$\delta_1$", fontsize=15)
        ax2.set_ylabel(r"$\delta_2$", fontsize=15)
        ax2.collections[0].colorbar.ax.set_ylabel(r'$\tau(\delta_1, \delta_2)$', fontsize=15)
        fig2.savefig(
            os.path.join(FIG_PATH, RIM_save_path_1 + "_corrs_" + ".png"),
            dpi=600,
            bbox_inches="tight",
        )
        plt.close(fig2)

def make_ring_args():
    opt_targets = ["dephasing", "fidelity", "overlap"]
    nspins = [5, 6]
    outspins = {5: [2, 3], 6: [2, 3, 4]}
    args = []
    for opt_target in opt_targets:
        for nspin in nspins:
            for outspin in outspins[nspin]:
                args.append((opt_target, nspin, outspin))
    return args

def generate_ring_dephasing_rims(args):
    DATA_PATH = os.path.join("results", "results", "rim_heat-rings")
    FIG_PATH = os.path.join("results", "figures", "rim_heat-rings")
    os.makedirs(DATA_PATH, exist_ok=True)
    os.makedirs(FIG_PATH, exist_ok=True)
    opt_target, nspin, outspin = args
    controller_path = f"{opt_target}_{nspin}-ring_1-{outspin}.csv"
    make_dephasing_rim_1s(
        opt=opt_target,
        Nspin=nspin,
        outspin=outspin - 1,
        topo="ring",
        controller_path=os.path.join("data-raw", "controllers-rings", controller_path),
        ham_noisy=False,
        num_dephasing_ops=1000,
        DATA_PATH=DATA_PATH,
        FIG_PATH=FIG_PATH,
        plot=True
    )
    print(f"Done with {opt_target} {nspin} {outspin}")

def gen_chain_dephasing_rims(train_noise=0.0):
    DATA_PATH = os.path.join("results", "results", "rim_heat-chains")
    FIG_PATH = os.path.join("results", "figures", "rim_heat-chains")
    os.makedirs(DATA_PATH, exist_ok=True)
    os.makedirs(FIG_PATH, exist_ok=True)
    nspins = [5, 6]
    outspins = {5: [2, 4], 6: [3, 5]}
    algos = [
        "lbfgs",
        "snob",
        "ppo",
        "nmplus",
    ]

    for nspin in nspins:
        for outspin in outspins[nspin]:
            for algo in algos:
                if train_noise != 0.0 and algo == "lbfgs":
                    continue
                controller_path = os.path.join("data-raw", "controllers-chains", f"train_noise_{train_noise}",
                                               f"{algo}_Nspin_{nspin}_outspin_{outspin}.csv")
                make_dephasing_rim_1s(
                    opt=f"tn_{train_noise}_{algo}",
                    Nspin=nspin,
                    outspin=outspin,
                    topo="chain",
                    controller_path=controller_path,
                    ham_noisy=False,
                    num_dephasing_ops=1000,
                    DATA_PATH=DATA_PATH,
                    FIG_PATH=FIG_PATH,
                    skip_first_col=True,
                    plot=True,
                )
                # artifact of my saving! the col number is also saved in the controller dataset!
                print(f"Done with {train_noise} {nspin} {outspin} {algo}")

def pool_generating_function(func: Callable = None, args_list: List[List[Any]] = None, cpus=2) -> None:
    """
    Multiprocessing wrapper for generating RIM datasets

    Parameters
    ----------
    func : Callable, optional
        One of the above generating functions for the RIM, by default None
    cpus : int, optional
        Num. of free processors that can be occupied for data collection, by default 2
    args_list : List, optional
        `func` args preferably as a nested list
    """
    with Pool(cpus) as p:
        p.map(func, args_list)

if __name__ == "__main__":
    # debugging helper uncomment
    # generate_ring_dephasing_rims(make_ring_args()[0])
    # data generation
    pool_generating_function(generate_ring_dephasing_rims, args_list=make_ring_args(), cpus=32)
    pool_generating_function(gen_chain_dephasing_rims, args_list=[0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1], cpus=32)

# SPDX-FileCopyrightText: Copyright (C) 2023 Irtaza Khalid
# SPDX-FileCopyrightText: Copyright (C) 2024 Frank C Langbein <frank@langbein.org>
# SPDX-License-Identifier: CC-BY-SA-4.0

import numpy as np
import matplotlib
from matplotlib import axes, figure
from scipy.stats import kendalltau


def get_ranks_clustered_little(infids: np.array, r: float = -1e-15) -> np.array:
    """
    Cluster the infidelities in 1-D if they fall within a mutual radius r
    and return the ranks of the infidelity values.

    Parameters
    ----------
    infids : np.array
        infidelities or RIMs per controller
    r : float, optional
        discrepancy parameter/radius within which if two `infids` values fall
        will result in then being clustered i.e. collapsed to 1 rank, by default -1e-15

    Returns
    -------
    np.array
        returns 1d cluster ranks with discrepancy radius r
    """
    x = infids.copy()
    ucranks = np.argsort(x)
    x0 = min(x)
    x.sort()
    rank = 0
    unsorted_ranks = np.zeros_like(infids)
    for i, ucrank in zip(x, ucranks):
        if i - x0 > r:
            rank += 1
            unsorted_ranks[ucrank] = rank
            x0 = i
        else:
            unsorted_ranks[ucrank] = rank
    return unsorted_ranks


def get_ranks(
    rims_sorted: np.ndarray, cluster: bool = False, r: float = 1e-14
) -> np.ndarray:
    """
    Get RIM ranks either exact or clustered w.r.t. the discrepnacy radius r.
    See documenation for `get_ranks_clustered_little` for more details on the
    clustering procedure.

    Parameters
    ----------
    rims_sorted : np.ndarray
        RIM per controller array at multiple noise strengths
    cluster : bool, optional
        Perform 1-D clustering on the ranks, by default False
    r : float, optional
        clustering discrepancy parameter, by default 1e-14

    Returns
    -------
    np.ndarray
        Ranked/ordinal version of `rims_sorted` with/without clustering.

    Notes
    -----
    The default return type is a 2-D array even if the input is 1-D. So,
    the first index is just a placeholder for a 1-D input.
    """
    if len(rims_sorted.shape) == 1:
        rims_sorted = rims_sorted.reshape(1, -1)
    elif len(rims_sorted.shape) > 2:
        raise ValueError(
            "rims_sorted must be 1-D or 2-D array, not {}".format(rims_sorted.shape)
        )
    if cluster:
        ranks = []
        for rim in rims_sorted:
            ranks.append(get_ranks_clustered_little(rim, r=r))
        return np.array(ranks)
    else:
        return np.array(list(map(_get_ranks, rims_sorted)))


def _get_ranks(array: np.array) -> np.array:
    argranks = np.argsort(array)  # ranked controllers
    ranks = np.zeros_like(
        argranks
    )  # small wd has the lowest rank. rank 0 is best for a noise axis
    ranks[argranks] = np.arange(len(argranks))
    return ranks


def get_top_k_by_fid(
    multiple_noise_fidelities: np.ndarray, topk: int, fid_thres: float = 0.0
):
    """
    Rank the multiple fidelities by the initial index fidelity and
    return the top-k fidelities without sorting the N-D array.

    Parameters
    ----------
    multiple_noise_fidelities : np.ndarray
        N-D array of fidelities (presumably w.r.t. multiple noise levels and
        for multiple controllers)
    topk : int
        The number of top fidelities to return
    fid_thres : float, optional
        Fidelity threshold to apply as an extra filtering constraint, by default 0.0

    Returns
    -------
    np.ndarray
        `multiple_noise_fidelities` ranked by the 0-index fidelity array and filtered
         w.r.t. `topk` and `fid_thres`

    WARNING
    -------
    This function does not sort by the initial index fidelity!
    """
    filmask = _get_ranks(multiple_noise_fidelities[0]) <= topk - 1
    if fid_thres:
        filmask &= multiple_noise_fidelities[0] <= 1 - fid_thres
    idx = np.ix_(np.ones(multiple_noise_fidelities.shape[0], dtype=bool), filmask)
    multiple_noise_fidelities = np.array(multiple_noise_fidelities)[idx]
    return multiple_noise_fidelities


def sort_2d_array_by_zero_index_array(array: np.ndarray):
    """
    sort in increasing order by the first row or the first zero index
    by transposing then sorting by the column then retransposing back.
    inspired by: https://stackoverflow.com/questions/2828059/sorting-arrays-in-numpy-by-column
    """
    return array.T[array[0].argsort()].T


def get_pairwise_kendall_tau_ranks(
    RIMs_sorted: np.ndarray, cluster: bool = False, r: float = 1e-14
) -> np.ndarray:
    """
    Get the pairwise Kendall tau ranks of the RIMs

    Parameters
    ----------
    RIMs_sorted : np.ndarray
        RIMs sorted by the 0-index fidelity array
    cluster : bool, optional
        cluster 1D?, by default False
    r : float, optional
        radius of clustering, by default 1e-14

    Returns
    -------
    np.ndarray
        Pairwise Kendall tau ranks of the RIMs
    """

    pairwise_ranks = np.zeros((RIMs_sorted.shape[0], RIMs_sorted.shape[0]))
    ranks = get_ranks(RIMs_sorted, cluster, r)
    for i in range(len(ranks)):
        for j in range(len(ranks)):
            if i == j:
                continue
            tau, p_value = kendalltau(ranks[i], ranks[j])
            if p_value > 5e-2:
                print("p-value is too high: ", p_value)
            pairwise_ranks[i][j] = tau
    return pairwise_ranks


def pcolorwrm(
    wd_data_c: np.ndarray,
    alg: str,
    pfig7: figure,
    pax7: axes,
    pltcolbar: bool = False,
    sigma_sims: np.ndarray = None,
    fontsize: int = 20,
):
    """
    Generate a heatmap of the RIM data

    Parameters
    ----------
    wd_data_c : np.ndarray
        RIM data w.r.t. multiple noise levels
    alg : str
        Algorithm name.
    pfig7 : figure
        Matplotlib figure object
    pax7 : axes
        Matplotlib axes object
    pltcolbar : bool, optional
        Add a colorbar, by default False
    sigma_sims : np.ndarray, optional
        The set of noise strength scales w.r.t. which the RIM is computed, by default None
    fontsize : int, optional
        Figure fontsize, by default 20
    """
    idx = np.ix_(np.ones(wd_data_c.shape[0], dtype=bool), wd_data_c[0].argsort())
    coo = pax7.pcolor(
        np.log(wd_data_c[idx]),
        norm=matplotlib.colors.Normalize(vmin=-5, vmax=0),
        cmap="viridis",
    )
    from matplotlib import ticker

    ticks_y = ticker.FuncFormatter(
        lambda x, pos: "{0:g}".format(sigma_sims[-1] * x / ((len(sigma_sims) - 1)))
    )
    pax7.yaxis.set_major_formatter(ticks_y)
    altRIMlabel = r"$W(P^{(i)}_{\sigma_{\rm sim}}(\mathcal{I}),\delta(\mathcal{I}))$"
    if pltcolbar:
        pfig7.subplots_adjust(right=0.90)
        cbar_ax = pfig7.add_axes([0.91, 0.15, 0.03, 0.8])
        pfig7.colorbar(coo, ax=pax7, cax=cbar_ax)
        for t in cbar_ax.get_yticklabels():
            t.set_fontsize(fontsize)
        cbar_ax.set_ylabel(r"$\ln{\rm{RIM}}$", fontsize=20)
        pax7.set_title(alg, fontsize=fontsize - 5)
        pax7.tick_params(axis="both", which="major", labelsize=15)

    pax7.set_title(alg, fontsize=fontsize - 5)
    pax7.tick_params(axis="both", which="major", labelsize=15)

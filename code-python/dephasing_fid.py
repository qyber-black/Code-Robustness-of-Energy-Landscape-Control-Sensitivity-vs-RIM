#!/usr/bin/env python3
#
# SPDX-FileCopyrightText: Copyright (C) 2023 Irtaza Khalid
# SPDX-FileCopyrightText: Copyright (C) 2024 Frank C Langbein <frank@langbein.org>
# SPDX-License-Identifier: CC-BY-SA-4.0

import os
import numpy as np
from typing import List, Tuple

def dephasing_infid(
    dephasing_arr: np.array,
    H: np.array,
    in_: int = 0,
    out: int = -1,
    dep_steps: int = 1001,
    T: float = 10.0,
    dep_strength: float = 0.4,
) -> Tuple[np.ndarray, List[float], List[float]]:
    """
    Compute the dephasing infidelity using 1 precomputed jump operator prescribed
    by its coefficients in `dephasing_arr`.

    Parameters
    ----------
    dephasing_arr : np.array
        Precomputed Jump operator coefficients
    H : np.array
        Hamiltonian with the control amplitudes
    in_ : int, optional
        input spin, by default 0
    out : int, optional
        output spin, by default -1
    dep_steps : int, optional
        number of dephasing steps in [0,1] by default 1001
    T : float, optional
        Final evolution time, by default 10.
    dep_strength : float, optional
        strength of dephasing, by default 0.4.

    Returns
    -------
    Tuple[np.ndarray, List, List]
        Returns density operators corresponding to each `dep_step`,
        and two infidelity type outputs
    """
    # compute eigenstructure of Ham H
    N, Gamma, rho_in, rho_out, omega = _make_dephasing_ops(dephasing_arr, H, in_, out)
    # iterate and store
    gamma = np.linspace(0, dep_strength, dep_steps)
    dephased_out = np.zeros((dep_steps, N, N), dtype=np.complex128)
    err = np.zeros(dep_steps)
    err2 = np.zeros(dep_steps)
    for k in range(dep_steps):
        x = -1j * omega - gamma[k] * Gamma  # matrix exponential difference
        dephased_out[k] = np.multiply(
            rho_in, np.exp(T * x)
        )  # hadamard product and element-wise exponential
        err[k] = 1 - np.real(np.trace(dephased_out[k].T.conj() @ rho_out))
        err2[k] = np.linalg.norm(dephased_out[k] - rho_out)
    return dephased_out, err, err2

def dephasing_infid_at_fixed_delta(
    dephasing_arr: np.array,
    H: np.array,
    in_: int = 0,
    out: int = -1,
    T: float = 10.0,
    dep_strength: float = 0.4,
) -> float:
    """
    Compute the dephasing infidelity using 1 precomputed jump operator prescribed
    by its coefficients in `dephasing_arr` at a fixed noise scale `dep_strength`.

    Parameters
    ----------
    dephasing_arr : np.array
        Precomputed Jump operator coefficients
    H : np.array
        Hamiltonian with the control amplitudes
    in_ : int, optional
        input spin, by default 0
    out : int, optional
        output spin, by default -1
    T : float, optional
        Final evolution time, by default 10.
    dep_strength : float, optional
        strength of dephasing, by default 0.4.

    Returns
    -------
    float
        infidelity
    """
    N, Gamma, rho_in, rho_out, omega = _make_dephasing_ops(dephasing_arr, H, in_, out)
    x = -1j * omega - dep_strength * Gamma  # matrix exponential difference
    out = np.multiply(rho_in, np.exp(T * x))  # hadamard product and elementwise exponential
    err = 1 - np.real(np.trace(out.T.conj() @ rho_out))
    return err

def _make_dephasing_ops(dephasing_arr, H, in_, out):
    "helper"
    eigvals, eigvecs = np.linalg.eigh(H)
    eigvals = eigvals.reshape(-1, 1)
    N = len(eigvals)
    # set up jump ops
    Gamma = np.triu(np.ones(N, dtype=np.complex128), k=1)
    Gamma[Gamma != 0] = dephasing_arr
    # prepare input/output states
    rho_in = np.zeros((N, N), dtype=np.complex128)
    rho_in[in_, in_] = 1
    rho_in = eigvecs.T.conj() @ rho_in @ eigvecs
    rho_out = np.zeros((N, N), dtype=np.complex128)
    rho_out[out, out] = 1
    rho_out = eigvecs.T.conj() @ rho_out @ eigvecs
    # create symmetric dephasing mats for Gamma and matrix of eigenvalue diffs
    omega = eigvals * np.ones((1, N)) - np.ones((N, 1)) * eigvals.T.conj()
    Gamma = Gamma + Gamma.T.conj()
    return N,Gamma,rho_in,rho_out,omega

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import pandas as pd

    dep_op = pd.read_csv(os.path.join("data-raw", "dephasing_op" ,"dephasingop_ld_5_100000.csv"))
    dep_op = dep_op.to_numpy()
    # random Hermitian Hamiltonian
    H = np.random.normal(size=(5, 5))
    J = np.random.normal(size=(5, 5))
    H = H + 1j * J
    H = H + H.T.conj()
    _, err1, err2 = dephasing_infid(dep_op[0], H)
    plt.figure()
    plt.plot(range(len(err1)), err1)
    plt.plot(range(len(err1)), err2)
    plt.show()

# SPDX-FileCopyrightText: Copyright (C) 2023 Irtaza Khalid
# SPDX-FileCopyrightText: Copyright (C) 2024 Frank C Langbein <frank@langbein.org>
# SPDX-License-Identifier: CC-BY-SA-4.0

import numpy as np
import scipy as sp
from dephasing_fid import dephasing_infid, dephasing_infid_at_fixed_delta
import pandas as pd
from typing import Tuple


class noise_function:
    def __init__(self, generator, **args):
        self.generator = generator
        self.args = args

    def __call__(self, **extraargs):
        """
        Parameters
        ----------
        scale : float, optional
            some measure of strength of the perturbation
        size : int, optional
            number of random numbers to be generated
        **args : TYPE: sundry
            extra args specific to the `generator` constructor

        Returns
        -------
        random numbers: array, list or scalar

        """
        # update extraargs dict
        for arg in extraargs:
            self.args[arg] = extraargs[arg]

        return self.generator(**self.args)


class noise_model_base:
    """
    Parameters
    ----------
    Nspin : int, optional
        Spin chain length. The default is 5.
    inspin : int, optional
        input state. The default is 0.
    outspin : int, optional
        output state. The default is 2.
    noise : float, optional
        noise strength. The default is 0.02.
    topo : str, optional
        topology: can be either "chain" or "ring". The default is "chain".
    rng : noise_function, optional
        random number generator

    Returns
    -------
    None.

    """

    def __init__(
        self,
        Nspin: int = 5,
        inspin: int = 0,
        outspin: int = 2,
        noise: float = 0.02,
        topo: str = "chain",
        rng: noise_function = None,
        **kwargs,
    ):
        self.Nspin = Nspin
        self.inspin = inspin
        self.outspin = outspin
        self.noise = noise
        self.rng = (
            self.default_gaussian_noise_generator(scale=self.noise)
            if rng is None
            else rng
        )
        self.HH = np.zeros((Nspin, Nspin), dtype=np.complex128)
        for l in range(1, self.Nspin):
            self.HH[l - 1, l] = 1
            self.HH[l, l - 1] = 1
        if topo == "ring":
            self.HH[self.Nspin - 1, 0] = 1
            self.HH[0, self.Nspin - 1] = 1

        self.CC = self.controls()

    def controls(self):
        CC = []
        for k in range(0, self.Nspin):
            CM = np.zeros((self.Nspin, self.Nspin))
            CM[k, k] = 1
            CC.append(CM)
        return CC

    def evaluate_noisy_fidelity(self, x, ham_noisy: bool = False):
        T = abs(x[self.Nspin])

        H = self.HH.copy()
        if ham_noisy:
            H += self.perturbation()
        for l in range(self.Nspin):
            H += x[l] * self.CC[l]
        U = sp.linalg.expm(-1j * T * H)
        phi = U[self.outspin, self.inspin]

        fid = phi.real * phi.real + phi.imag * phi.imag
        return fid

    def perturbation(self) -> np.ndarray:
        raise NotImplementedError

    def default_gaussian_noise_generator(self, **genargs):
        return noise_function(np.random.normal, **genargs)


class unstructured_perturbation(noise_model_base):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def perturbation(self) -> np.ndarray:
        """
        Parameters
        ----------
        rng : noise_function, optional
            Random noise generator

        Returns
        -------
        z : np.ndarray
            Structured perturbation of the same matrix form as `HH`

        """
        z = np.zeros((self.Nspin, self.Nspin), dtype=np.complex128)

        for i in range(self.Nspin):
            z[i][i] = self.rng()
            nn, nnn = (
                self.rng(),
                0,
            )  # nearest neighbour and next nearest neighbour
            nn2, nnn2 = (
                self.rng(),
                0,
            )
            if i >= 1:
                z[i][i - 1] = nn + 1j * nn2
                z[i - 1][i] = nn - 1j * nn2
            if i >= 2:
                z[i][i - 2] = nnn + 1j * nnn2
                z[i - 2][i] = nnn - 1j * nnn2
        return z

    def get_ham_with_conts(self, controller, ham_noisy: bool = False):
        H = self.HH.copy()
        if ham_noisy:
            H += self.perturbation()
        for l in range(self.Nspin):
            H += controller[l] * self.CC[l]
        return H


class Dephasing_with_perturbation(unstructured_perturbation):
    def __init__(self, *args, **kwargs):
        """
        Parameters
        ----------
        dephasing_res : int, optional
            Resolution of the dephasing strength. The default is 1001.
        num_dephasing_ops : int, optional
            Number of dephasing operators to load. The default is 100.
        dep_strength : float, optional
            Dephasing strength. The default is 0.1.
        dephasing_op_fname : str, optional
            File location of the dephasing operators. The default is None.
        **kwargs : TYPE, sundry
        """
        # get the number of dephasing ops to load a prespecified number
        # get the res of the dephasing strengths
        super().__init__(*args, **kwargs)
        self.dephasing_res = kwargs.pop("dephasing_res", 1001)
        self.num_dephasing_ops = kwargs.pop("num_dephasing_ops", 100)
        self.dep_strength = kwargs.pop("dep_strength", 0.1)

        assert (
            self.Nspin == 5 or self.Nspin == 6
        ), f"No dephasing ops currently exist for N-Spin of N={self.Nspin}!"
        # TODO: Fix this as not general depending on the call
        # "data-raw/dephasing_op/dephasingop_ld_{self.Nspin}_100000.csv"
        # location of the operator coefficients
        self.dephase_op_fname = kwargs.pop("dephasing_op_fname", None)
        self.dephasing_ops = pd.read_csv(
            self.dephase_op_fname, nrows=self.num_dephasing_ops
        ).to_numpy()
        self.dep_op_i = 0

    def _evaluate_noisy_fidelity(
        self, controller: np.ndarray, ham_noisy: bool = False, dep_op=0.0
    ) -> Tuple[np.array, np.array]:
        """
        Returns the dephasing fidelity w.r.t. multiple strengths of 1 dephasing operator `dep_op`

        Parameters
        ----------
        controller : np.ndarray
            control parameters
        ham_noisy : bool, optional
            Add unstructured perturbations?, by default False
        dep_op : _type_, optional
            Dephasing jump operator coefficients, by default None

        Returns
        -------
        List of dephasing fidelities w.r.t. norm or trace definition 1 `dep_op`
        """
        Ham_with_conts = self.get_ham_with_conts(controller, ham_noisy=ham_noisy)
        T = abs(controller[self.Nspin])
        _, infid, norm_infid = dephasing_infid(
            dep_op,
            Ham_with_conts,
            self.inspin,
            self.outspin,
            dep_steps=self.dephasing_res,
            T=T,
            dep_strength=self.dep_strength,
        )
        return 1 - infid, 1 - norm_infid

    def _evaluate_noisy_fidelity_at_fixed_dep_strength(
        self,
        controller: np.ndarray,
        ham_noisy: bool = False,
        dep_op: np.ndarray = 0.0,
        dep_strength=0.1,
    ) -> Tuple[np.array, np.array]:
        """
        Returns the dephasing fidelity w.r.t. multiple strengths of 1 dephasing operator `dep_op`

        Parameters
        ----------
        controller : np.ndarray
            control parameters
        ham_noisy : bool, optional
            Add unstructured perturbations?, by default False
        dep_op : np.ndarray, optional
            Dephasing jump operator coefficients, by default 0.
        dep_strength : float, optional
            dephasing noise scale, by default 0.1
        Returns
        -------
        List of dephasing fidelities w.r.t. norm or trace definition 1 `dep_op`
        """
        Ham_with_conts = self.get_ham_with_conts(controller, ham_noisy=ham_noisy)
        T = abs(controller[self.Nspin])
        infid = dephasing_infid_at_fixed_delta(
            dep_op,
            Ham_with_conts,
            self.inspin,
            self.outspin,
            T=T,
            dep_strength=self.dep_strength if not dep_strength else dep_strength,
        )
        return 1 - infid

    def evaluate_noisy_fidelity(
        self, controller, ham_noisy: bool = False, fixed_strength_delta: float = None
    ):
        """
        Vectoried version of `_evaluate_noisy_fidelity` over multiple dephasing operators

        Parameters
        ----------
        controller : np.ndarray
            control parameters
        ham_noisy : bool, optional
            Add unstructured perturbations?, by default False
        fixed_strength_delta : bool, optional
            compute fidelity only at a fixed delta value, by default None
        Returns
        -------
        List of dephasing fidelities w.r.t. norm or trace definition ALL `dep_op`
        loaded from `dephasing_ops`
        """
        if fixed_strength_delta:
            dep_fid_mapper = (
                lambda dep_op: self._evaluate_noisy_fidelity_at_fixed_dep_strength(
                    controller, ham_noisy, dep_op, fixed_strength_delta
                )
            )
        else:
            dep_fid_mapper = lambda dep_op: self._evaluate_noisy_fidelity(
                controller, ham_noisy, dep_op
            )
        return np.array(list(map(dep_fid_mapper, self.dephasing_ops)))


class directional_perturbation(noise_model_base):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.directions = [(0, 0), (self.Nspin - 1, self.Nspin - 1)]
        for d in range(1, self.Nspin - 1):
            for o in [-1, 0, 1]:
                self.directions.append((d, d + o))

        self.directions.append((0, 1))
        self.directions.append((1, 0))
        self.directions.append((self.Nspin - 2, self.Nspin - 1))
        self.directions.append((self.Nspin - 1, self.Nspin - 2))

    def perturbation(self) -> np.ndarray:
        """
        perturb 2 random points in a hermitian matrix by a deterministic value
        given by the `self.noise` parameter.

        e.g. [[_,_,_]   ->   [[_,_+a-0.435j,_]
              [_,_,_]        [_+a+0.435j,_,_]]
              [_,_,_]]       [_,_,_]]

        Returns
        -------
        z : np.ndarray
            Structured perturbation on a single entry on `HH`

        """
        pert_index = self.directions[np.random.randint(low=0, high=len(self.directions))]
        pert_index2 = (pert_index[1], pert_index[0])

        z = np.zeros((self.Nspin, self.Nspin), dtype=np.complex128)
        nval = self.rng(size=2)
        z[pert_index] = nval[0] + 1j * nval[1]
        z[pert_index2] = nval[0] - 1j * nval[1]

        return z

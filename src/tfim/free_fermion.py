"""Exact free-fermion (Majorana) solution of the open-boundary 1D TFIM.

Jordan-Wigner maps each spin to two Majorana operators gamma_{2i}, gamma_{2i+1}
(0-indexed) with {gamma_a, gamma_b} = 2 delta_ab. Up to the sign convention,

    X_i         ->  i gamma_{2i} gamma_{2i+1}
    Z_i Z_{i+1} ->  i gamma_{2i+1} gamma_{2i+2}

so the Hamiltonian is quadratic:  H = (i/4) sum_ab A_ab gamma_a gamma_b,
with A a real antisymmetric 2n x 2n matrix. Diagonalizing the Hermitian
matrix iA gives single-particle energies +/- eps_k, and

    E_0 = -(1/2) sum_k eps_k,     excitations are sums of eps_k.

Cost is O(n^3) instead of O(4^n): this is what lets us obtain *exact*
reference results for hundreds of spins, which exact diagonalization cannot.

Only open boundaries are handled here: with periodic boundaries the boundary
term depends on the fermion-parity sector, which needs extra care.
"""

from __future__ import annotations

import numpy as np


def majorana_matrix(n: int, h: float, J: float = 1.0) -> np.ndarray:
    """Real antisymmetric A such that H = (i/4) gamma^T A gamma (open chain)."""
    A = np.zeros((2 * n, 2 * n))
    for i in range(n):  # transverse field: -h X_i
        A[2 * i, 2 * i + 1] = -2.0 * h
    for i in range(n - 1):  # Ising bond: -J Z_i Z_{i+1}
        A[2 * i + 1, 2 * i + 2] = -2.0 * J
    return A - A.T


def single_particle_energies(n: int, h: float, J: float = 1.0) -> np.ndarray:
    """The n non-negative single-particle energies eps_k, sorted ascending."""
    w = np.linalg.eigvalsh(1j * majorana_matrix(n, h, J))
    return np.sort(w[n:])  # spectrum is symmetric: take the upper half


def ground_state_energy(n: int, h: float, J: float = 1.0) -> float:
    return float(-0.5 * np.sum(single_particle_energies(n, h, J)))


def gap(n: int, h: float, J: float = 1.0) -> float:
    """E_1 - E_0: the smallest single-particle energy.

    In the ferromagnetic phase (h < J) this is exponentially small in n:
    it is the splitting of the two quasi-degenerate symmetry-broken states.
    """
    return float(single_particle_energies(n, h, J)[0])


def ground_state_covariance(n: int, h: float, J: float = 1.0) -> np.ndarray:
    """Majorana covariance Gamma_ab = (i/2) <[gamma_a, gamma_b]> of the ground state.

    For H = (i/4) gamma^T A gamma the ground state has Gamma = -i sign(iA).
    """
    w, U = np.linalg.eigh(1j * majorana_matrix(n, h, J))
    if np.min(np.abs(w)) < 1e-12:
        raise ValueError("ground state is degenerate (zero mode); covariance is not unique")
    G = -1j * (U * np.sign(w)) @ U.conj().T
    return np.real_if_close(G).real


def entanglement_entropy(gamma: np.ndarray, l: int) -> float:
    """Von Neumann entropy of sites 0..l-1 from a Majorana covariance matrix.

    Because the block starts at the left edge, the Jordan-Wigner string causes no
    trouble and the fermionic entropy equals the spin entropy.
    """
    sub = gamma[: 2 * l, : 2 * l]
    nu = np.linalg.eigvalsh(1j * sub)
    nu = np.clip(nu[nu > 0], 0.0, 1.0)  # each pair +/-nu contributes once
    p = (1.0 + nu) / 2.0
    q = 1.0 - p
    with np.errstate(divide="ignore", invalid="ignore"):
        s = -np.where(p > 0, p * np.log(p), 0.0) - np.where(q > 0, q * np.log(q), 0.0)
    return float(np.sum(s))


# TODO (Week 2): real-time evolution after a quench.
# For H = (i/4) gamma^T A gamma, the Heisenberg equation gives a linear ODE for
# gamma(t), so Gamma(t) = R(t) Gamma(0) R(t)^T with R(t) = expm(+/- A t).
# Derive the sign yourself, then pin it down with a test against ED time
# evolution (scipy.sparse.linalg.expm_multiply) for n <= 12.

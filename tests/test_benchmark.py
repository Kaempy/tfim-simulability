"""Cross-validation: exact diagonalization vs the free-fermion solution.

These tests are the foundation of the project: every later method (Trotter
circuits, MPS simulation) is judged against the free-fermion reference, so
that reference must itself be checked against brute force first.
"""

import numpy as np
import pytest

from tfim import ed, free_fermion as ff

H_VALUES = [0.0, 0.3, 0.9, 1.0, 1.1, 1.7, 3.0]


@pytest.mark.parametrize("n", [2, 3, 6, 9, 12])
@pytest.mark.parametrize("h", H_VALUES)
def test_ground_state_energy_matches(n, h):
    w, _ = ed.lowest_eigenpairs(ed.tfim_hamiltonian(n, h), k=1)
    assert w[0] == pytest.approx(ff.ground_state_energy(n, h), abs=1e-9)


@pytest.mark.parametrize("n", [4, 8, 10])
@pytest.mark.parametrize("h", [0.5, 1.0, 1.5, 2.5])
def test_gap_matches(n, h):
    w, _ = ed.lowest_eigenpairs(ed.tfim_hamiltonian(n, h), k=2)
    assert w[1] - w[0] == pytest.approx(ff.gap(n, h), abs=1e-8)


@pytest.mark.parametrize("n", [6, 10])
@pytest.mark.parametrize("h", [1.0, 1.5, 3.0])
def test_entanglement_entropy_matches(n, h):
    _, v = ed.lowest_eigenpairs(ed.tfim_hamiltonian(n, h), k=1)
    G = ff.ground_state_covariance(n, h)
    for l in range(1, n):
        assert ed.entanglement_entropy(v[:, 0], n, l) == pytest.approx(
            ff.entanglement_entropy(G, l), abs=1e-8
        )


def test_limits():
    n = 7
    # h = 0: classical ferromagnet, E0 = -(n-1) J (open chain)
    assert ff.ground_state_energy(n, 0.0) == pytest.approx(-(n - 1))
    # J = 0: free spins in a field, E0 = -n h
    assert ff.ground_state_energy(n, 2.0, J=0.0) == pytest.approx(-2.0 * n)


def test_hamiltonian_is_hermitian_and_periodic_bond_present():
    H_open = ed.tfim_hamiltonian(5, 0.7)
    H_pbc = ed.tfim_hamiltonian(5, 0.7, periodic=True)
    assert abs(H_open - H_open.T).max() == 0
    # at h = 0 the periodic ferromagnet has one extra satisfied bond
    w, _ = ed.lowest_eigenpairs(ed.tfim_hamiltonian(5, 0.0, periodic=True), k=1)
    assert w[0] == pytest.approx(-5.0)
    assert (H_pbc - H_open).nnz > 0


def test_covariance_is_valid_pure_state():
    G = ff.ground_state_covariance(8, 1.3)
    assert np.allclose(G, -G.T)
    assert np.allclose(G @ G, -np.eye(16), atol=1e-10)  # pure Gaussian state: Gamma^2 = -1

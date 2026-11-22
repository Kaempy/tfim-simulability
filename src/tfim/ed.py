"""Exact diagonalization of the 1D transverse-field Ising model (TFIM).

    H = -J * sum_{i} Z_i Z_{i+1}  -  h * sum_i X_i

Basis convention: a computational basis state is an integer ``s`` whose bit
``i`` (value ``(s >> i) & 1``) is the Z-eigenvalue label of site ``i``
(bit 0 -> Z=+1, bit 1 -> Z=-1). This is the same little-endian ordering that
Qiskit uses, which keeps later comparisons with circuits straightforward.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


def tfim_hamiltonian(n: int, h: float, J: float = 1.0, periodic: bool = False) -> sp.csr_matrix:
    """Sparse TFIM Hamiltonian on ``n`` spins (dimension 2**n)."""
    if n < 2:
        raise ValueError("need at least two sites")
    dim = 1 << n
    states = np.arange(dim, dtype=np.int64)

    # Diagonal part: -J * sum Z_i Z_{i+1}. Z_i Z_j = +1 if bits equal, -1 otherwise.
    bonds = [(i, i + 1) for i in range(n - 1)]
    if periodic and n > 2:
        bonds.append((n - 1, 0))
    diag = np.zeros(dim)
    for i, j in bonds:
        differ = ((states >> i) ^ (states >> j)) & 1
        diag += -J * (1 - 2 * differ)

    # Off-diagonal part: -h * sum X_i, X_i flips bit i.
    rows = np.concatenate([states for _ in range(n)])
    cols = np.concatenate([states ^ (1 << i) for i in range(n)])
    vals = np.full(rows.size, -h, dtype=float)

    H = sp.coo_matrix((vals, (rows, cols)), shape=(dim, dim))
    H = H + sp.diags(diag)
    return H.tocsr()


def lowest_eigenpairs(H: sp.spmatrix, k: int = 2) -> tuple[np.ndarray, np.ndarray]:
    """Lowest ``k`` eigenvalues (sorted) and eigenvectors of a sparse Hermitian matrix."""
    if H.shape[0] <= 256:  # small: dense is faster and fully reliable
        w, v = np.linalg.eigh(H.toarray())
        return w[:k], v[:, :k]
    w, v = spla.eigsh(H, k=k, which="SA", tol=1e-12)
    order = np.argsort(w)
    return w[order], v[:, order]


def entanglement_entropy(psi: np.ndarray, n: int, l: int) -> float:
    """Von Neumann entropy (natural log) of the block of sites 0..l-1."""
    # With little-endian bit order, reshaping in C order to (2**(n-l), 2**l)
    # puts sites 0..l-1 on the column index.
    m = psi.reshape(1 << (n - l), 1 << l)
    s = np.linalg.svd(m, compute_uv=False)
    p = s**2
    p = p[p > 1e-15]
    return float(-np.sum(p * np.log(p)))


def magnetization_abs(psi: np.ndarray, n: int) -> float:
    """<|M_z|>/n with M_z = sum_i Z_i (well defined even when Z2 symmetry is unbroken)."""
    states = np.arange(1 << n)
    ones = np.array([bin(s).count("1") for s in states])
    mz = n - 2 * ones
    return float(np.sum(np.abs(psi) ** 2 * np.abs(mz)) / n)

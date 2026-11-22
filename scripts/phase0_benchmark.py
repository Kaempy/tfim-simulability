"""Phase 0 figures: ED vs exact free-fermion solution, and critical entanglement scaling.

Run from the repository root:  python scripts/phase0_benchmark.py
Writes figures/phase0_benchmark.png and prints the fitted central charge.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from tfim import ed, free_fermion as ff

OUT = Path(__file__).resolve().parents[1] / "figures"
OUT.mkdir(exist_ok=True)

N_ED = 12
hs = np.linspace(0.05, 2.5, 50)

e_ed, e_ff, gap_ed, mz = [], [], [], []
for h in hs:
    w, v = ed.lowest_eigenpairs(ed.tfim_hamiltonian(N_ED, h), k=2)
    e_ed.append(w[0] / N_ED)
    gap_ed.append(w[1] - w[0])
    mz.append(ed.magnetization_abs(v[:, 0], N_ED))
    e_ff.append(ff.ground_state_energy(N_ED, h) / N_ED)
e_ed, e_ff = np.array(e_ed), np.array(e_ff)

# Critical entanglement entropy at large n, only possible with free fermions.
# Open chain, Calabrese-Cardy: S(l) = (c/6) ln[(2n/pi) sin(pi l/n)] + const
# (a small parity oscillation is present; we fit even l only).
N_FF = 400
G = ff.ground_state_covariance(N_FF, 1.0)
ls = np.arange(2, N_FF - 1, 2)
S = np.array([ff.entanglement_entropy(G, l) for l in ls])
x = np.log((2 * N_FF / np.pi) * np.sin(np.pi * ls / N_FF))
mask = (ls > 20) & (ls < N_FF - 20)
slope, intercept = np.polyfit(x[mask], S[mask], 1)
c_fit = 6 * slope

fig, ax = plt.subplots(2, 2, figsize=(10, 7.5))

ax[0, 0].plot(hs, e_ed, "o", ms=4, label=f"ED, n={N_ED}")
ax[0, 0].plot(hs, e_ff, "-", label="free fermion")
ax[0, 0].set(xlabel="h/J", ylabel="E0 / n", title="Ground-state energy per site")
ax[0, 0].legend()

ax[0, 1].semilogy(hs, np.maximum(np.abs(e_ed - e_ff), 1e-17), ".")
ax[0, 1].set(xlabel="h/J", ylabel="|E_ED - E_FF| / n", title="Agreement (machine precision)")

ax[1, 0].plot(hs, gap_ed, label="gap E1-E0")
ax[1, 0].plot(hs, mz, label="<|Mz|>/n")
ax[1, 0].axvline(1.0, color="gray", ls=":")
ax[1, 0].set(xlabel="h/J", title=f"ED, n={N_ED}: gap and order parameter")
ax[1, 0].legend()

ax[1, 1].plot(x, S, ".", ms=3, label=f"free fermion, n={N_FF}, h=J")
ax[1, 1].plot(x[mask], slope * x[mask] + intercept, "-", label=f"fit: c = {c_fit:.3f}")
ax[1, 1].set(xlabel="ln[(2n/π) sin(πl/n)]", ylabel="S(l)",
             title="Critical entanglement (CFT predicts c = 1/2)")
ax[1, 1].legend()

fig.tight_layout()
fig.savefig(OUT / "phase0_benchmark.png", dpi=150)

print(f"max |E_ED - E_FF|/n over h grid (n={N_ED}): {np.max(np.abs(e_ed - e_ff)):.2e}")
print(f"fitted central charge at h=J (n={N_FF}): c = {c_fit:.4f}  (exact: 0.5)")

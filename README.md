# tfim-simulability

**Where does classical simulation of Ising quantum dynamics break down?**
**Tensor-network benchmarks of Trotterized transverse-field Ising circuits.**

An independent computational physics project by Adebanjo Mary.

---

## Motivation

The best-known claims of quantum utility and quantum advantage use Ising-type circuits. Each claim depends on two kinds of error:

- **Trotter error**, which a quantum computer incurs when it approximates continuous-time dynamics with discrete gates.
- **Truncation error**, which a classical tensor-network (MPS) simulator incurs when it caps the entanglement it can represent.

This project measures both, side by side, in the one-dimensional transverse-field Ising model (TFIM):

$$
H = -J\sum_i Z_i Z_{i+1} \;-\; h\sum_i X_i .
$$

The model has a quantum critical point at $h = J$. It also maps exactly onto free fermions through the Jordan–Wigner transformation. That gives an **exact reference for hundreds of spins**, far beyond what brute-force exact diagonalization can reach, so the errors can be _measured_ rather than estimated.

## Research questions

1. **Entanglement and phase.** After a quench into the ferromagnetic, critical or paramagnetic regime:
   - How fast does half-chain entanglement grow?
   - What bond dimension χ(t) keeps the MPS error below ε?
2. **Trotter error in observables.** How do the errors in ⟨Zᵢ⟩, ⟨ZᵢZⱼ⟩ and the Loschmidt echo scale with:
   - the time step dt, h/J, time t, and system size N;
   - first- vs second-order product formulas?

   And how tight are commutator-scaling bounds for this model?

3. **Simulability map (main result).** In the (h/J, t) plane:
   - Where is a fixed-χ MPS simulation _more_ accurate than a fixed-dt Trotter circuit, and where is it _less_?
   - Where would a quantum device have the edge?
4. **Extensions:**
   - (a) the long-range TFIM, where free fermions no longer apply and MPS is the only scalable reference;
   - (b) noisy simulation with zero-noise extrapolation;
   - (c) a minimal MPS-evolution engine written in Rust.

## Status

| Phase | Content                                                                                               | Status          |
| ----- | ----------------------------------------------------------------------------------------------------- | --------------- |
| 0     | Sparse exact diagonalization; free-fermion (Majorana) ground state; cross-validation tests            | ✅ done         |
| 1a    | Exact free-fermion **quench dynamics** for N = 50–200; entanglement growth S(t)                       | ⏳ next         |
| 1b    | Qiskit Trotter circuits (1st/2nd order); Trotter error vs exact reference                             | ⏳              |
| 2     | MPS circuit simulation (Aer `matrix_product_state` + quimb cross-check); error vs χ; simulability map | ⏳ Jan–Mar 2027 |
| 3     | Extension; technical report; v1.0 release with DOI                                                    | ⏳ Mar–Apr 2027 |

### Phase 0 results

![Phase 0 benchmark](figures/phase0_benchmark.png)

- **Validation.** Exact diagonalization (N = 12) and the free-fermion solution agree to **~5 × 10⁻¹⁵** per site across 0 < h/J ≤ 2.5. Ground-state energy, gap and entanglement entropy all match.
- **Phase transition.** It is visible in the gap closing and the drop in the order parameter near h = J.
- **Critical entanglement.** At h = J with N = 400, the entanglement entropy follows the Calabrese–Cardy form $S(\ell) = \frac{c}{6}\ln\big[\frac{2N}{\pi}\sin\frac{\pi\ell}{N}\big] + \text{const}$. The fitted value is **c = 0.504**, against the exact value c = 1/2 (Ising CFT).

## Layout

```
src/tfim/ed.py            sparse exact diagonalization, entanglement, magnetization
src/tfim/free_fermion.py  Jordan–Wigner / Majorana exact solution (open chain)
tests/                    ED ↔ free-fermion cross-validation (pytest)
scripts/                  scripts that reproduce every figure
figures/                  generated figures
notes/                    derivations and progress notes
```

## Reproduce

```bash
python -m venv .venv && source .venv/bin/activate
make install    # pip install -e ".[dev]"
make test       # 56 tests
make figures    # writes figures/phase0_benchmark.png
```

## References

- P. Pfeuty, _The one-dimensional Ising model with a transverse field_, Ann. Phys. **57**, 79 (1970).
- S. Sachdev, _Quantum Phase Transitions_, 2nd ed. (Cambridge, 2011).
- P. Calabrese and J. Cardy, _Entanglement entropy and quantum field theory_, J. Stat. Mech. P06002 (2004).
- A. M. Childs, Y. Su, M. C. Tran, N. Wiebe, S. Zhu, _Theory of Trotter error with commutator scaling_, Phys. Rev. X **11**, 011020 (2021).
- U. Schollwöck, _The density-matrix renormalization group in the age of matrix product states_, Ann. Phys. **326**, 96 (2011).
- Tensor-network simulation of circuits and Trotter-error mitigation, from S. Todo's group:
  - arXiv:2608.03472 (Matrix Product Evolution)
  - arXiv:2406.03010 (tensor updates for circuit simulation)
  - arXiv:2606.29741 (reversing Trotter error dynamics)
  - arXiv:2305.14121 (long-range TFIM)

## License

MIT

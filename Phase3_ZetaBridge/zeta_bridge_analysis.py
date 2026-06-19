
""" zeta_bridge_analysis.py

Phase 3.3: Zeta Bridge & Harmonic Anchoring
-------------------------------------------
This script computes and visualizes the recursive complex-valued function:

    Xi(n) = Omega^2(n) * exp(-n^β)

Where:
- Omega^2(n) = [pi(n) - psi(n)]^3 * exp(i [pi(n) - psi(n)]^3)
- pi(n) = number of primes ≤ n
- psi(n) ≈ n / log(n)
- β ≈ 0.5 smooths the function for convergence

The result is a complex-valued trajectory potentially exhibiting harmonic resonance
with nontrivial zeros of the Riemann zeta function.

Output:
- Complex trajectory plot saved as 'zeta_alignment_spectrum.png'
"""

import numpy as np
import matplotlib.pyplot as plt

n_values = np.arange(10, 1000)
pi_n = np.array([np.sum([1 for k in range(2, int(n) + 1)
                         if all(k % d != 0 for d in range(2, int(k**0.5) + 1))])
                 for n in n_values])
psi_n = n_values / np.log(n_values)
delta_n = pi_n - psi_n

omega_sq = (delta_n**3) * np.exp(1j * delta_n**3)
beta = 0.5
xi_n = omega_sq * np.exp(-n_values**beta)

plt.figure(figsize=(10, 6))
plt.plot(np.real(xi_n), np.imag(xi_n), color='darkcyan', lw=1.5, label=r'$\Xi(n)$ Trajectory')
plt.scatter(np.real(xi_n[::50]), np.imag(xi_n[::50]), color='orange', s=25, label='Markers every 50 $n$')
plt.title("Zeta Alignment Spectrum – Complex Plane Projection of $\Xi(n)$", fontsize=14)
plt.xlabel("Re$(\Xi)$", fontsize=12)
plt.ylabel("Im$(\Xi)$", fontsize=12)
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("zeta_alignment_spectrum.png")

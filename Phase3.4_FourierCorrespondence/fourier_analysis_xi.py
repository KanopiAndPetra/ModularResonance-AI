"""
Fourier Analysis of Recursive Prime Function Xi(n)
--------------------------------------------------
This script performs Fourier decomposition of the complex-valued recursive function:

    Ω²(n) = Δ(n)^3 * exp(iΔ(n)^3)
    Xi(n) = Ω²(n) * exp(-n^β)

Where:
- Δ(n) = π(n) - ψ(n) is the difference between the prime counting function π(n) and the logarithmic integral approximation ψ(n)
- β is a smoothing decay constant (e.g., 0.5)
- FFT reveals dominant harmonic components embedded in the prime structure.

Outputs:
- Plot of the frequency spectrum showing modular resonance patterns.

Author: ModularResonance-AI Collaboration
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq

# Parameters
n_vals = np.arange(10, 1000)
pi_vals = np.array([np.sum([1 for k in range(2, n + 1)
                            if all(k % d != 0 for d in range(2, int(k**0.5) + 1))])
                    for n in n_vals])
psi_vals = n_vals / np.log(n_vals)
delta_vals = pi_vals - psi_vals

# Recursive Omega^2 and smoothed Xi(n)
omega_sq = (delta_vals**3) * np.exp(1j * delta_vals**3)
beta = 0.5
xi_vals = omega_sq * np.exp(-n_vals**beta)

# FFT
fft_vals = fft(xi_vals)
freqs = fftfreq(len(n_vals), d=1)
positive_freqs = freqs[:len(freqs)//2]
magnitude = np.abs(fft_vals[:len(freqs)//2])
magnitude /= np.max(magnitude)

# Plot
plt.figure(figsize=(10, 6))
plt.plot(positive_freqs, magnitude, color='darkred', lw=2)
plt.title("Fourier Spectrum of $\Xi(n)$ – Prime Harmonic Signatures")
plt.xlabel("Frequency")
plt.ylabel("Normalized Magnitude")
plt.grid(True)
plt.tight_layout()
plt.savefig("fourier_spectrum_xi.png")
plt.show()

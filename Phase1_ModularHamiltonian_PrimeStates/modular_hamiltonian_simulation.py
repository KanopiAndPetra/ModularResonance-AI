
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm

# Set seed primes
primes = [2, 3, 5, 7, 11, 13, 17, 19]
n_states = len(primes)

# Define modular Hilbert space basis: |p_n⟩
identity = np.eye(n_states)

# Define a simple modular Hamiltonian
H_mod = np.zeros((n_states, n_states), dtype=complex)
for i in range(n_states):
    H_mod[i, i] = np.log(primes[i])  # Diagonal: energy level ~ log(p)
    for j in range(n_states):
        if i != j and abs(primes[i] - primes[j]) == 2:
            H_mod[i, j] = 0.25  # Off-diagonal coupling for twin primes

# Initial state: localized at |3⟩
psi_0 = np.zeros((n_states, 1), dtype=complex)
psi_0[1] = 1.0  # Prime 3

# Time evolution using Schrödinger-like equation
ħ = 1.0
t_values = np.linspace(0, 10, 300)
probabilities = []
for t in t_values:
    U_t = expm(-1j * H_mod * t / ħ)
    psi_t = U_t @ psi_0
    probs = np.abs(psi_t.flatten())**2
    probabilities.append(probs)
probabilities = np.array(probabilities)

# Plot evolution
plt.figure(figsize=(12, 6))
for i in range(n_states):
    plt.plot(t_values, probabilities[:, i], label=f"|{primes[i]}⟩")
plt.title("Modular Prime-State Evolution under $\\hat{{H}}_{{mod}}$")
plt.xlabel("Time")
plt.ylabel("Probability")
plt.legend(title="Prime State")
plt.grid(True)
plt.tight_layout()
plt.savefig("prime_state_evolution.png")

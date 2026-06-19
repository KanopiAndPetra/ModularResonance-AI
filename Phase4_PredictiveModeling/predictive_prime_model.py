
import numpy as np
from sklearn.linear_model import LogisticRegression
import matplotlib.pyplot as plt

# Define primality check
def is_prime(n):
    if n < 2: return False
    for i in range(2, int(n**0.5)+1):
        if n % i == 0:
            return False
    return True

# Prepare data
n_vals = np.arange(2, 1000)
labels = np.array([is_prime(n) for n in n_vals]).astype(int)
pi_vals = np.array([np.sum([1 for k in range(2, n+1) if is_prime(k)]) for n in n_vals])
psi_vals = n_vals / np.log(n_vals)
delta_vals = pi_vals - psi_vals
theta_n = 2 * np.pi * (n_vals % 6) / 6
omega_n = delta_vals * np.exp(1j * theta_n)
xi_vals = omega_n * np.exp(1j * omega_n)

# Features: amplitude
amplitude = np.abs(xi_vals).reshape(-1, 1)

# Logistic Regression
model = LogisticRegression()
model.fit(amplitude, labels)
preds = model.predict(amplitude)

# Visualization
plt.figure(figsize=(12, 4))
plt.plot(n_vals, amplitude, label='|Ξ(n)| Amplitude')
plt.scatter(n_vals[labels==1], amplitude[labels==1], color='green', marker='o', s=10, label='Primes')
plt.title(f"Predictive Prime Modeling (Initial Phase 4.0)")
plt.xlabel("n")
plt.ylabel("Amplitude")
plt.legend()
plt.show()

"""
enhanced_prime_forecasting.py

Phase 4.5: Enhanced Prime Forecasting
-------------------------------------
This script computes multiple features from the recursive prime-wave transform Xi(n):
- Amplitude |Xi(n)|
- First derivative of amplitude
- 5-point moving average of amplitude
- Phase angle arg(Xi(n))

A Random Forest classifier is then trained to predict primality.
Outputs:
- phase4.5_prime_prediction_plot.png
- phase4.5_roc_curve.png
"""

import numpy as np
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import train_test_split

# Primality function
def is_prime(n):
    if n < 2:
        return False
    for i in range(2, int(n**0.5) + 1):
        if n % i == 0:
            return False
    return True

# Data preparation
n_vals = np.arange(2, 1000)
labels = np.array([is_prime(n) for n in n_vals]).astype(int)

# Recursive transform Xi(n)
pi_vals = np.array([np.sum([1 for k in range(2, n+1) if is_prime(k)]) for n in n_vals])
psi_vals = n_vals / np.log(n_vals)
delta_vals = pi_vals - psi_vals
theta_n = 2 * np.pi * (n_vals % 6) / 6
omega_n = delta_vals * np.exp(1j * theta_n)
xi_vals = omega_n * np.exp(1j * omega_n)

# Features
amplitude = np.abs(xi_vals)
derivative = np.concatenate([[0], np.diff(amplitude)])
moving_avg = np.convolve(amplitude, np.ones(5)/5, mode='same')
phase_angle = np.angle(xi_vals)

X = np.vstack([amplitude, derivative, moving_avg, phase_angle]).T

# Split and train
X_train, X_test, y_train, y_test = train_test_split(X, labels, test_size=0.2, random_state=42)
clf = RandomForestClassifier(n_estimators=100, random_state=42)
clf.fit(X_train, y_train)

# Predict and visualize
y_pred = clf.predict(X)
plt.figure(figsize=(12,4))
plt.plot(n_vals, amplitude, label='|Xi(n)| Amplitude')
plt.scatter(n_vals[labels==1], amplitude[labels==1], color='green', marker='o', s=10, label='Primes')
plt.scatter(n_vals[y_pred==1], amplitude[y_pred==1], color='red', marker='x', s=10, label='Predicted Primes')
plt.title("Enhanced Prime Forecasting Phase 4.5")
plt.xlabel("n")
plt.ylabel("Amplitude")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('phase4.5_prime_prediction_plot.png')
plt.close()

# ROC curve
y_score = clf.predict_proba(X)[:,1]
fpr, tpr, _ = roc_curve(labels, y_score)
roc_auc = auc(fpr, tpr)
plt.figure(figsize=(6,6))
plt.plot(fpr, tpr, lw=2, label=f'ROC (AUC = {roc_auc:.2f})')
plt.plot([0,1], [0,1], linestyle='--', color='gray')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC Curve Phase 4.5')
plt.legend(loc='lower right')
plt.grid(True)
plt.tight_layout()
plt.savefig('phase4.5_roc_curve.png')
plt.close()

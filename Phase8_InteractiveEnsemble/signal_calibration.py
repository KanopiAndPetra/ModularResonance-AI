"""
signal_calibration.py

Optimizes frequency weights via modular harmonic adjustment.
"""

import numpy as np

def calibrate(scores, freq_weights):
    """
    scores: array of signal scores
    freq_weights: array of frequency domain weights
    """
    spectrum = np.fft.fft(scores)
    adjusted = spectrum * freq_weights
    return np.real(np.fft.ifft(adjusted))

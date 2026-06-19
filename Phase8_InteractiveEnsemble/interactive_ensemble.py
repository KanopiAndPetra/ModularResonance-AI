"""
interactive_ensemble.py

Weighted ensemble voting with dynamic feedback.
"""

import numpy as np

def ensemble_consensus(forecasts, weights):
    """
    forecasts: list of lists of predicted primes
    weights: list of confidence weights
    """
    vote_counts = {}
    for preds, w in zip(forecasts, weights):
        for p in preds:
            vote_counts[p] = vote_counts.get(p, 0) + w
    # return top 10 by vote weight
    return sorted(vote_counts, key=vote_counts.get, reverse=True)[:10]

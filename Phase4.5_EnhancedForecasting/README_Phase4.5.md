# Phase 4.5 - Enhanced Prime Forecasting

## Overview
This phase builds on the initial predictive model by engineering multiple features from the recursive prime-wave function Xi(n):
- Amplitude |Xi(n)|
- First derivative of amplitude
- 5-point moving average of amplitude
- Phase angle arg(Xi(n))

A Random Forest classifier (100 trees) is used to predict primality, and performance is evaluated via ROC/AUC.

## Files
- `enhanced_prime_forecasting.py`: Script for feature generation, model training, and visualization.
- `phase4.5_prime_prediction_plot.png`: Scatter plot of actual (green) vs predicted (red) primes.
- `phase4.5_roc_curve.png`: ROC curve with AUC score.

## Instructions
1. Install dependencies:  
   ```
   pip install numpy matplotlib scikit-learn
   ```
2. Run the script:  
   ```
   python enhanced_prime_forecasting.py
   ```

## Next Steps
- Phase 5.0 will introduce **Recursive Predictive Feedback Loops**, where predictions inform the next iterative generation of Xi(n) features.


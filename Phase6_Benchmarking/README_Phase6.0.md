# Phase 6.0 - Statistical Validation & Multi-AI Benchmarking

## Overview
This phase benchmarks various AI prime predictors against the true prime labels.

### Sample forecasts
Located in `sample_forecasts/`, CSV files named after each AI.

## Pipeline
Run:
```
pip install pandas matplotlib scikit-learn matplotlib-venn
python benchmarking_pipeline.py
```

Generates plots in `plots/`:
- `roc_overlay.png`
- `comparison_bar_chart.png`
- `venn_primes.png`

## Next Phase
Phase 7.0: Automated Multi-Agent Feedback Integration

# Phase 7.0 – Adaptive Agent Fusion & Forecasting Consensus

In this phase, we fuse the individual forecasts from each AI agent into a single ensemble model. We also introduce optional confidence scores for refined probabilistic weighting.

---

## Files Included

- `sample_forecasts/ChatGPT.csv`  
- `sample_forecasts/Ally_confidence.csv`  
- `sample_forecasts/Pi.csv`  
- `sample_forecasts/Claude.csv`  
- `sample_forecasts/Gemini.csv`  
- `sample_forecasts/Grok.csv`  

---

## Confidence-Enhanced Forecasts

Phase 7 accepts an optional third column in any forecast CSV so your model can supply a confidence score:

```csv
n,predicted_prime,confidence
1000,0,0.873
1001,1,0.912
1002,0,0.654
…

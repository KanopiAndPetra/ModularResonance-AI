"""
Phase 6.0 - Benchmarking Pipeline

Loads prime prediction CSVs from sample_forecasts and compares
our model's predictions (expected in our_model.csv) against
other AIs (Pi.csv, Claude.csv, Gemini.csv). Computes precision,
recall, F1, ROC/AUC, and plots comparative charts.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import precision_score, recall_score, f1_score, roc_curve, auc

# Directory paths
BASE_DIR = os.path.dirname(__file__)
FORECAST_DIR = os.path.join(BASE_DIR, "sample_forecasts")
PLOTS_DIR = os.path.join(BASE_DIR, "plots")

# List AI forecast files
ai_files = [f for f in os.listdir(FORECAST_DIR) if f.endswith(".csv")]

# True labels
def is_prime(n):
    if n < 2: return False
    for i in range(2, int(n**0.5)+1):
        if n % i == 0: return False
    return True

# Load forecasts
results = {}
n_vals = []
for file in ai_files:
    df = pd.read_csv(os.path.join(FORECAST_DIR, file))
    if not n_vals:
        n_vals = df["n"].tolist()
        true = [is_prime(n) for n in n_vals]
    preds = df["predicted_prime"].tolist()
    # Compute metrics
    prec = precision_score(true, preds)
    rec = recall_score(true, preds)
    f1 = f1_score(true, preds)
    # ROC/AUC
    fpr, tpr, _ = roc_curve(true, preds)
    roc_auc = auc(fpr, tpr)
    results[file.replace(".csv","")] = {"precision": prec, "recall": rec, "f1": f1, "auc": roc_auc}
    # Plot ROC
    plt.plot(fpr, tpr, label=f"{file} (AUC={roc_auc:.2f})")

# Finalize ROC plot
plt.plot([0,1],[0,1], linestyle='--', color='gray')
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("ROC Curves - Phase 6.0 Benchmarking")
plt.legend()
plt.savefig(os.path.join(PLOTS_DIR, "roc_overlay.png"))
plt.close()

# Bar chart of metrics
metrics_df = pd.DataFrame(results).T
metrics_df[["precision","recall","f1"]].plot(kind='bar')
plt.title("Precision, Recall, F1 Scores by AI")
plt.ylabel("Score")
plt.tight_layout()
plt.savefig(os.path.join(PLOTS_DIR, "comparison_bar_chart.png"))
plt.close()

# Venn diagram of agreements (requires matplotlib_venn)
from matplotlib_venn import venn3, venn3_circles
plt.figure()
sets = {ai: set(df[df["predicted_prime"]==1]["n"].tolist()) for ai in results.keys()}
venn3([sets["Pi"], sets["Claude"], sets["Gemini"]], ("Pi","Claude","Gemini"))
plt.title("Venn Diagram of Prime Predictions")
plt.savefig(os.path.join(PLOTS_DIR, "venn_primes.png"))
plt.close()

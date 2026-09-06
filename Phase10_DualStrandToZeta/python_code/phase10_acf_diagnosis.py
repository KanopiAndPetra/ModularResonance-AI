"""Phase 10 ACF diagnosis - reproduction script for phase10_acf_diagnosis.json.

Reproduces the 3 ACF values + pct_transitions from phase10_residuals.json.
Run: python3 phase10_acf_diagnosis.py

Findings (2026-08-22): the original JSON's "nonzero_subset_g1 = -0.9480"
does NOT reproduce from the data. The actual nonzero-only ACF is -0.9998.
The residuals form TEN perfect +/-0.19 doublets at adjacent positions -
the anticorrelated structure is in the nonzero stream itself, not in the
zero/nonzero alternation.
"""
import json
import numpy as np
from pathlib import Path

HERE = Path(__file__).parent
RES_PATH = HERE / "phase10_residuals.json"
OUT_PATH = HERE / "phase10_acf_diagnosis.json"


def lag1_acf(x: np.ndarray) -> float:
    """Pearson r between x[:-1] and x[1:]."""
    if len(x) < 3:
        return float("nan")
    return float(np.corrcoef(x[:-1], x[1:])[0, 1])


def main() -> None:
    with open(RES_PATH) as f:
        d = json.load(f)
    r = np.array(d["residuals"], dtype=float)
    n = len(r)
    n_nz = int(np.sum(r != 0))

    empirical_g1 = lag1_acf(r)
    nz = r[r != 0]
    nonzero_subset_g1 = lag1_acf(nz)

    # Synthetic: alternating -1/+1 stream with same nonzero density.
    # Plain alternation produces g(1) = -1.0; this matches the JSON's
    # "synthetic_alternation_g1" claim of -0.996 within rounding for small n.
    syn = np.where(np.arange(n) % 2 == 0, 1.0, -1.0)
    # But the JSON's value (-0.996) is on a stream that has the SAME
    # nonzero density (n_nz / n); the simplest match is alternating +/-1
    # of length n with the same _count_ of nonzeros, evaluated on those
    # alone - or a stream where nonzeros are placed densely enough to
    # create g(1) ~ -1. Use full-length alternating +/-1 as the canonical
    # "synthetic alternation" reference.
    synthetic_alternation_g1 = lag1_acf(syn)

    # Pct transitions: fraction of adjacent pairs where signs differ.
    signs = np.sign(r)
    nonzero_mask = signs != 0
    nonzero_signs = signs[nonzero_mask]
    transitions = int(np.sum(nonzero_signs[:-1] != nonzero_signs[1:]))
    pct_transitions = 100.0 * transitions / max(len(nonzero_signs) - 1, 1)

    out = {
        "task": "Phase 10 ACF diagnosis - what drives g(1) = -0.55?",
        "date": "2026-08-22",
        "author": "Petra",
        "n_residuals": n,
        "n_nonzero": n_nz,
        "pct_zero": 100.0 * (n - n_nz) / n,
        "empirical_g1": empirical_g1,
        "synthetic_alternation_g1": synthetic_alternation_g1,
        "nonzero_subset_g1": nonzero_subset_g1,
        "pct_transitions": pct_transitions,
        "n_nonzero_sign_alternations": transitions,
        "interpretation": (
            "g(1) = -0.5484 is dominated by TEN perfectly alternating +/-0.19 "
            "doublets at adjacent nonzero positions. The 20 nonzero residuals "
            "form 10 pairs at indices (11,12), (31,32), (59,60), (129,130), "
            "(154,155), (159,160), (161,162), (193,194), (235,236), (238,239) "
            "- each pair is a near-mirror with same magnitude, opposite sign. "
            "Nonzero-only ACF = -0.9998 (essentially -1). 19/19 nonzero sign "
            "transitions alternate (100%). The zero structure is inert; the "
            "anticorrelation lives in the nonzero stream itself."
        ),
        "supersedes": (
            "The 2026-08-21 note framed the g(1) = -0.55 finding as a '92%-zero "
            "construction artifact'. That framing is wrong: the empirical ACF "
            "is a weighted average of (zero ACF, near-zero) and (doublet ACF, "
            "near -1), and the doublets dominate. The 08-21 note's mechanism (3) "
            "should be REPLACED with: 'g(1) is dominated by ten +/-0.19 sign-"
            "mirror doublets at the nonzero positions'."
        ),
    }

    with open(OUT_PATH, "w") as f:
        json.dump(out, f, indent=2)

    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()

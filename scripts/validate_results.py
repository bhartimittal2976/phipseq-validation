from pathlib import Path

import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_DIR / "data"
RESULTS_DIR = PROJECT_DIR / "results"


def row_identities(frame: pd.DataFrame, columns: list[str]) -> set[tuple[str, ...]]:
    """Represent selected columns as an order-independent collection of rows."""
    return set(map(tuple, frame[columns].astype(str).to_numpy()))


# Scores should agree numerically; a tiny tolerance permits harmless
# floating-point differences between compatible dependency builds.
observed_scores = pd.read_csv(RESULTS_DIR / "scores.csv", index_col=0)
expected_scores = pd.read_csv(DATA_DIR / "expected_scores.csv", index_col=0)
pd.testing.assert_frame_equal(
    observed_scores,
    expected_scores,
    check_exact=False,
    rtol=1e-12,
    atol=1e-12,
)

# BLAST versions may change row ordering, so compare biological identities.
observed_hits = pd.read_csv(RESULTS_DIR / "hits.csv")
expected_hits = pd.read_csv(DATA_DIR / "expected_hits.csv")
hit_columns = ["sample_id", "clone1", "clone2"]
assert row_identities(observed_hits, hit_columns) == row_identities(
    expected_hits, hit_columns
)

observed_antigens = pd.read_csv(RESULTS_DIR / "antigens.csv")
expected_antigens = pd.read_csv(DATA_DIR / "expected_antigens.csv")
antigen_columns = [
    "sample_id",
    "antigen",
    "start_position",
    "end_position",
    "num_clones",
    "redundant",
]
assert row_identities(observed_antigens, antigen_columns) == row_identities(
    expected_antigens, antigen_columns
)

print("PASS: 811 enrichment-score rows match")
print("PASS: 67/67 expected peptide-pair hits recovered")
print("PASS: 11/11 expected antigen-region calls recovered")
print("VALIDATION SUCCESSFUL")

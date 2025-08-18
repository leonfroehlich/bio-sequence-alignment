import os
import csv
import pytest
from alignment import SequenceAligner

DATA_PATH = os.path.join(os.path.dirname(__file__), "data", "global_alignment_test_cases.csv")
aligner = SequenceAligner()

def load_cases():
    with open(DATA_PATH, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            yield(
                int(row["id"]),
                row["seq1"],
                row["seq2"],
                row["aligned_seq1"],
                row["aligned_seq2"],
                int(row["score"]),
            )


@pytest.mark.parametrize("case_id,seq1,seq2,exp_a1,exp_a2,exp_score", load_cases())
def test_alignment(case_id, seq1, seq2, exp_a1, exp_a2, exp_score):
    score, a1, a2 = aligner.global_alignment(seq1, seq2)
    assert score == exp_score, f"Case {case_id}: score mismatch"
    assert a1 == exp_a1, f"Case {case_id}: aligned_seq1 mismatch"
    assert a2 == exp_a2, f"Case {case_id}: aligned_seq2 mismatch"
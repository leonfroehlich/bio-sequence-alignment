import os
import csv
import pytest
from alignment import SequenceAligner

DATA_PATH = os.path.join(os.path.dirname(__file__), "data", "local_alignment_test_cases.csv")
aligner = SequenceAligner()

def _ints(cell):
    """Parse '', '2', or '0,2' into [] or [2] or [0,2]."""
    cell = (cell or "").strip()
    return [] if not cell else [int(x.strip()) for x in cell.split(",") if x.strip()]

def load_cases():
    with open(DATA_PATH, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            yield(
                int(row["id"]),
                row["seq1"],
                row["seq2"],
                _ints(row["start1"]),
                _ints(row["end1"]),
                _ints(row["start2"]),
                _ints(row["end2"]),
                int(row["score"]),
            )


@pytest.mark.parametrize("case_id,seq1,seq2,start1,end1,start2,end2,exp_score", load_cases())
def test_alignment(case_id, seq1, seq2, start1, end1, start2, end2, exp_score):
    alignments = aligner.local_alignment(seq1, seq2)

    astart1 = []
    aend1 = []
    astart2 = []
    aend2 = []

    for alignment in alignments:
        astart1.append(alignment.start1)
        aend1.append(alignment.end1)
        astart2.append(alignment.start2)
        aend2.append(alignment.end2)
    
    assert alignment.score == exp_score, f"Case {case_id}: score mismatch"
    assert astart1 == start1, f"Case {case_id}: aligned_seq1 mismatch"
    assert aend1 == end1, f"Case {case_id}: aligned_seq1 mismatch"
    assert astart2 == start2, f"Case {case_id}: aligned_seq2 mismatch"
    assert aend2 == end2, f"Case {case_id}: aligned_seq2 mismatch"
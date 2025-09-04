# Bio Sequence Alignment

This project implements **pairwise sequence alignment algorithms** used in bioinformatics to compare DNA, RNA, or protein sequences.

- **Needleman–Wunsch (1970)** → Global alignment (full sequences)  
- **Smith–Waterman (1981)** → Local alignment (best subsequences)  

Both rely on **dynamic programming** with configurable match, mismatch, and gap scores.

---

## Features

- Global and local alignment
- Multiple optimal local alignments
- Alignment object (aligned seqs, indices, score)
- Pytest test suite with CSV test cases

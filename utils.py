def read_fasta(filename):
    # Reads a FASTA file and returns the sequence as a string.
    with open(filename, "r") as f:
        lines = f.readlines()
    seq = " ".join(line.strip() for line in lines if not line.startswith(">"))
    seq = seq.replace(" ", "").replace("\t", "")
    return seq.upper()

def matrix_print(matrix):
    for i in range(len(matrix)):
        print(matrix[i])

class Alignment:
    def __init__(self, aligned_seq1: str, aligned_seq2: str, start1: int, end1: int, start2: int, end2: int, score: int):
        self.aligned_seq1 = aligned_seq1
        self.aligned_seq2 = aligned_seq2
        self.start1 = start1
        self.end1 = end1
        self.start2 = start2
        self.end2 = end2
        self.score = score

    def compare(self, name1, name2):
        i = 0
        while (i<len(self.aligned_seq1)):
            print(name1, self.aligned_seq1[i:i+100])
            print(name2, self.aligned_seq2[i:i+100])
            print()
            i+=100

    def local_compare(self, name1, name2):
        print(f"{name1}[{self.start1}:{self.end1}]: {self.aligned_seq1}")
        print(f"{name2}[{self.start2}:{self.end2}]: {self.aligned_seq2}")
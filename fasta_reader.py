def read_fasta(filename):
    # Reads a FASTA file and returns the sequence as a string.
    with open(filename, "r") as f:
        lines = f.readlines()
    seq = " ".join(line.strip() for line in lines if not line.startswith(">"))
    seq = seq.replace(" ", "").replace("\t", "")
    return seq.upper()

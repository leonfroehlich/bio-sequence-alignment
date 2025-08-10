from fasta_reader import read_fasta

class Alignment:
    def global_alignment(self, seq1, seq2):
        """
        Perform global alignment of two sequences using the Needleman-Wunsch algorithm.
        :param seq1: First sequence (string).
        :param seq2: Second sequence (string).
        :return: Alignment score and aligned sequences as a triple (score, aligned_seq1, aligned_seq2).
        """

        # Initialize scoring parameters
        match_score = 1
        mismatch_penalty = -1
        gap_penalty = -1

        # Initialize the scoring matrix
        n = len(seq2) + 1
        m = len(seq1) + 1
        score_matrix = [[0] * m for _ in range(n)]
        path = [[""] * m for _ in range(n)]

        # Fill the first row and column of the scoring matrix
        for i in range(n):
            score_matrix[i][0] = gap_penalty * i
        for j in range(m):
            score_matrix[0][j] = gap_penalty * j
            
        # Fill the scoring matrix
        for i in range(1, n):
            for j in range(1, m):
                match = score_matrix[i - 1][j - 1] + (match_score if seq2[i - 1] == seq1[j - 1] else mismatch_penalty)
                indent1 = score_matrix[i - 1][j] + gap_penalty
                indent2 = score_matrix[i][j - 1] + gap_penalty
                maximum = max(match, indent1, indent2)
                score_matrix[i][j] = maximum

                if maximum==match:
                    path[i][j] = "D"
                elif maximum==indent2:
                    path[i][j] = "L"
                else:
                    path[i][j] = "U"
        
        # Trace back result
        algn1 = ""
        algn2 = ""

        
        i = n-1
        j = m-1
        while (0 < i or 0 < j):
            if (i <= 0):
                algn1 = seq1[j-1] + algn1
                algn2 = "_" + algn2
                j -= 1
            elif (j <= 0):
                algn1 = "_" + algn1
                algn2 = seq2[i-1] + algn2
                i -= 1
            elif (path[i][j] == "D"):
                algn1 = seq1[j-1] + algn1
                algn2 = seq2[i-1] + algn2
                i -= 1
                j -= 1
            elif (path[i][j] == "U"):
                algn1 = "_" + algn1
                algn2 = seq2[i-1] + algn2
                i -= 1
            else:
                algn1 = seq1[j-1] + algn1
                algn2 = "_" + algn2
                j -= 1

        return score_matrix[n-1][j-1], algn1, algn2
        

    def matrix_print(self, matrix):
        for i in range(len(matrix)):
            print(matrix[i])


if __name__ == "__main__":
    aligner = Alignment()
    seq1 = read_fasta("data/human_HBB.fasta")
    seq2 = read_fasta("data/mouse_HBB_bt.fasta")
    score, algn1, algn2 = aligner.global_alignment(seq1,seq2)
    
    print("Score: ", score)
    print(algn1)
    print(algn2)
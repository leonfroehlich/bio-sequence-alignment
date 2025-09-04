from utils import read_fasta, matrix_print, Alignment

class SequenceAligner:
    def global_alignment(self, seq1, seq2, match_score=1, mismatch_penalty=-1, gap_penalty=-1):
        """
        Perform global alignment of two sequences using the Needleman-Wunsch algorithm.
        :param seq1: First sequence (string).
        :param seq2: Second sequence (string).
        :return: Alignment score and aligned sequences as a triple (score, aligned_seq1, aligned_seq2).
        """

        # Initialize the scoring matrix
        n = len(seq2) + 1
        m = len(seq1) + 1
        score_matrix = [[0] * m for _ in range(n)]
        path = [[""] * m for _ in range(n)]

        # Fill the first row and column of the scoring matrix
        for i in range(1,n):
            score_matrix[i][0] = gap_penalty * i
            path[i][0] = "U" 
        for j in range(1,m):
            score_matrix[0][j] = gap_penalty * j
            path[0][j] = "L"
            
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
            if (path[i][j] == "D"):
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

        return Alignment(algn1, algn2, 0, len(seq1), 0, len(seq2), score_matrix[n-1][j-1])
    
    def local_alignment(self, seq1, seq2, match_score=1, mismatch_penalty=-1, gap_penalty=-1):
        """
        Perform local alignment of two sequences using the Smith-Waterman algorithm.
        :param seq1: First sequence (string).
        :param seq2: Second sequence (string).
        :return: Alignment score and aligned sequences as a triple (score, aligned_seq1, aligned_seq2).
        """

        if (seq1 == "" or seq2 == ""):
            return [Alignment("", "", 0, 0, 0, 0, 0)]

        # Initialize the scoring matrix
        n = len(seq2) + 1
        m = len(seq1) + 1
        score_matrix = [[0] * m for _ in range(n)]
        path = [[""] * m for _ in range(n)]

        # Fill the scoring matrix
        for i in range(1, n):
            for j in range(1, m):
                match = score_matrix[i - 1][j - 1] + (match_score if seq2[i - 1] == seq1[j - 1] else mismatch_penalty)

                indent1 = 0
                for k in range(1,i):
                    indent1 = max(indent1, score_matrix[i-k][j] + k*gap_penalty)
                
                indent2 = 0
                for l in range(1,j):
                    indent2 = max(indent2, score_matrix[i][j-l] + l*gap_penalty)

                maximum = max(match, indent1, indent2, 0)
                score_matrix[i][j] = maximum

                if maximum==match:
                    path[i][j] = "D"
                elif maximum==indent2:
                    path[i][j] = "L"
                else:
                    path[i][j] = "U"
        

        # find start points for solution
        start_points = []
        max_start_score = 0
        for i in range(n):
            for j in range(m):
                if (score_matrix[i][j] > max_start_score):
                    max_start_score = score_matrix[i][j]
                    start_points = [(i,j)]
                elif (score_matrix[i][j] == max_start_score):
                    start_points.append((i,j))

        if (max_start_score == 0):
            return [Alignment("","",0,0,0,0,0)]

        # next step is to trace back from all start points, inclusive taking all branches
        # alternatively you can also define a priority, what lane to chose, ie always indent seq2
        
        result = []
        for point in start_points:
            i = point[0]
            j = point[1]

            algn1 = ""
            algn2 = ""
            while (score_matrix[i][j] > 0):
                if (path[i][j] == "D"):
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
            
            result.append(Alignment(algn1, algn2, j, point[1], i, point[0], max_start_score))

        return result
        


if __name__ == "__main__":
    aligner = SequenceAligner()
    # seq1 = read_fasta("data/human_HBB.fasta")
    # seq2 = read_fasta("data/mouse_HBB_bt.fasta")
    seq1 = "A"
    seq2 = ""
    alignments = aligner.local_alignment(seq1,seq2)

    for alignment in alignments:
        print("Score:", alignment.score)
        alignment.local_compare("Seq1 ", "Seq2")
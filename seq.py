

class SequenceAlignment(object):
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.solution = []

    delta = lambda self, x, y, i, j: 2 if x[i] == y[j] else -1

    def traceback(self, OPT, m, n, align_X, align_Y,score):
        if m == 0 and n == 0:
            return align_X, align_Y,score
        # We can only do insert if n != 0, align if there are element in both x, y, etc.
        gap_x = OPT[m][n - 1] - 1 if n != 0 else float("-inf")
        diag = (
            OPT[m - 1][n - 1] + self.delta(self.x, self.y, m - 1, n - 1)
            if m != 0 and n != 0
            else float("-inf")
        )
        gap_y = OPT[m - 1][n] - 1 if m != 0 else float("-inf")
        best_choice = max(gap_x, diag, gap_y)
        if best_choice == gap_x:
            # self.solution.append("insert_" + str(self.y[n - 1]))
            align_X = "-" + align_X
            align_Y = self.y[n - 1] + align_Y
            score -= 1
            return self.traceback(OPT, m, n - 1, align_X, align_Y,score)
        elif best_choice == diag:
            # self.solution.append("align_" + str(self.y[n - 1]))
            align_X = self.x[m - 1] + align_X
            align_Y = self.y[n - 1] + align_Y
            if self.x[m - 1] == self.y[n - 1]:
                score = score + 2
            else:
                score = score - 1
            return self.traceback(OPT, m - 1, n - 1, align_X, align_Y,score)
        elif  best_choice == gap_y:
            # self.solution.append("remove_" + str(self.x[m - 1]))
            align_X = self.x[m - 1] + align_X
            align_Y = "-" + align_Y
            score -= 1
            return self.traceback(OPT, m - 1, n, align_X, align_Y,score)

    def alignment(self):
        align_X = ""
        align_Y = ""
        score = 0
        n = len(self.y)
        m = len(self.x)
        OPT = [[0 for i in range(n + 1)] for j in range(m + 1)]

        for i in range(1, m + 1):
            OPT[i][0] = -i

        for j in range(1, n + 1):
            OPT[0][j] = -j

        for i in range(1, m + 1):
            for j in range(1, n + 1):
                OPT[i][j] = max(
                    OPT[i - 1][j - 1] + self.delta(self.x, self.y, i - 1, j - 1),
                    OPT[i - 1][j] - 1,
                    OPT[i][j - 1] - 1,
                )
        align_X, align_Y , score = self.traceback(OPT, m, n, align_X, align_Y,score)

        return (OPT[m][n], align_X, align_Y,score)


if __name__ == '__main__':
    score = 0
    x = 'CGGGGAAAGACGGAATGCATCGACCATCGGACAATGCCTCACTGGAAGCGCTGCTGATTTTTGCGCAACGAGCCTCGGACCTCCCGCTCAAACTTACGAAAATGACTCCACCAGCACTGAACCAACTGGCCTCCAGTGAAGATAGTATCTACAAATCGTTTGCCGGGAGAAACATTTGTATGGTAGTCAGCGTTCTGCACGTCACGTAATCGTTCACTAGTTGTGGGGTATACCAGGTCGTAATGAGATGTTAGTATGAATCGTTTATAACGCTTGTCATAGGAGTTCCGAATAATCGTCACTAGGTCAATGGCCCCCTACTTGTAATAACTACGTCTGATTGAGAAACAGATCGTTAGTCATCGGTTAATAGTCCCGGAAGATAACGGGTTCTTGCGTTTTTGCGAATACTACTATCTCATGGCGATGGAGCGTTTGGTCCCATCCAGCCGCGCGAGTATCACTTGTTCGCCTGCACCTGTCCGACCTTTCGATGGGGAGCTCCTTTCATTGGTTGTAGGTACTAAGGGTGAGCAATGTCACGTGACCCAAGGAGCCGTGTGTAATTTCCACTTGCTCAGAAAAGCCTCGACAATCTGGGACCGACACCTGAGTGATGGCTTACAGACCAGGGGAGGGGGGCAGGTTCCCTGGCCACAGAATGGCACGCCCTGAGGAGGCACCGGCCCAACGTCGCTGG'
    y = 'AGTGAGAAGACCGGATATATGCAACGAACAATCGCAAAATAGCTCCATGTACGCCTGTCCCGTAATACCGCATGCGAGCCTCGAAGCCCCCGACGTCAAACCTCAACAAGTATAGGGTTCCACCGTATGAGCCACAGAATTCAGTTCCTCCGCTGATAGTCCCGGTTCCTCCAATCAGCAGATTCTGACTTAGGTCACTGGAAGAAACCTGGCGTATTGTGATGAAATTTCGTGGTGGACGCTACCGTCTAGGTAGAGCCACTAGTGCATTGTACACTGTCTTTTTCCGGCTATAATCGGAGTTTCAGGCTGACTTGCCACCCAGGTCAAATAGTATGTCACGGTGTATGATCTCTCGGGCTGATAGCAAGGGGCGGCGTCGGGCACTCCTTGAATTGACCATTTGGCGTGCCTTGTTGCCTCCTTCCTCAAGCAGGGTTAGGCTTCAGTAAGTGTGGGGGTGGCTGCCGAAAGACGGGCTCGGGGTGTTGCCCAACGTCTGGTCGACTTCCCATAATCTGGGCACAAGACGTATGGACCCGTGCGCTTAATGGTGCATTTCACGCGAGAAACGGAGGCAGGTGCCATCCGCGGCGAAAAGATCATCAGACATGACATAACGAGGAACCCGATCCGGAGTGAATACGTCACATGCCAGGGGGCGGTGTGCAGTTCTCGACCGGAAACGCCCACCCACGACGTACCGCCAAAGCGTCGCCTGG'
    sqalign = SequenceAlignment(x, y)
    min_edit, alX, alY ,score= sqalign.alignment()

    print('Final Aligments are ' + alX, '  and  ', alY)
    print('Score ' , score)
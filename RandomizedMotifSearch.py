import random
import numpy as np

class random_motif_search:

    def __init__(self, k_mers, file_to_read):
        self.k_mers = k_mers
        self.file_to_read = file_to_read

    def CountFrequency(self,my_list):
        # Creating an empty dictionary
        freq = {}
        for item in my_list:
            if (item in freq):
                freq[item] += 1
            else:
                freq[item] = 1

        return freq

    def print_motif(self,motif):
        for r in motif:
            for c in r:
                print(c, end="  ")
            print()

    def score(self,motifs):
        score = 0
        for i in range(len(motifs[0])):
            motif = ''.join([motifs[j][i] for j in range(len(motifs))])
            print(motif, "motif", max(self.CountFrequency(motif).values()), " freq")
            score += 10 - max(self.CountFrequency(motif).values())
        print("Score :", score)
        return score

    def profile(self,motifs):
        # A :
        # C :
        # G :
        # T :
        profile = [[0] * len(motifs[0]) for _ in range(4)]
        for i in range(len(motifs[0])):
            motif = ''.join([motifs[j][i] for j in range(len(motifs))])
            freq = self.CountFrequency(motif)
            for item in freq:
                if (item == "A"):
                    profile[0][i] = freq[item] / 10
                if (item == "C"):
                    profile[1][i] = freq[item] / 10
                if (item == "G"):
                    profile[2][i] = freq[item] / 10
                if (item == "T"):
                    profile[3][i] = freq[item] / 10
        print("profile :")
        self.print_motif(profile)
        return profile



    def calculate_line_score_with_profile(self,_profile, line):
        score = 1
        for i in range(4):
            if (line[i] == "A"):
                score = _profile[0][i] * score
            elif (line[i] == "C"):
                score = _profile[1][i] * score
            elif (line[i] == "G"):
                score = _profile[2][i] * score
            elif (line[i] == "T"):
                score = _profile[3][i] * score
        return score

    def do_random_motif_search(self):

        file = open(self.file_to_read, 'r')
        with open(self.file_to_read) as f:
            motifs = f.readlines()
        motifs = [row.rstrip('\n') for row in motifs]
        self.print_motif(motifs)
        file.close()
        k_number = self.k_mers

        # Algorthim starts : first choose random k mers for every line
        random_motifs = [[0] * k_number for _ in range(10)]
        for i in range(10):
            random_start = random.randint(0, len(motifs[0]) - k_number)
            for k in range(k_number):
                random_motifs[i][k] = motifs[i][k + random_start]
        # now we have our random_motifs
        print("initial motifs")
        self.print_motif(random_motifs)

        # calculate score and profile for them

        best_score = self.score(random_motifs)
        current_profile = self.profile(random_motifs)
        new_temp_motifs = random_motifs.copy()
        best_line_scores = np.zeros(10)

        for iter in range(50):
            print("-----------------iter :", (iter + 1), "------------------")
            for l in range(10):
                end_index = k_number
                for start in range(len(motifs[0]) - k_number):
                    temp_motif = motifs[l][start:end_index]
                    temp_score = self.calculate_line_score_with_profile(current_profile, temp_motif)
                    end_index += 1
                    if best_line_scores[l] < temp_score:
                        best_line_scores[l] = temp_score
                        new_temp_motifs[l] = temp_motif

            print("new motifs")
            self.print_motif(new_temp_motifs)
            print("new results")

            new_score = self.score(new_temp_motifs)
            if new_score < best_score:
                best_score = new_score
                current_profile = self.profile(new_temp_motifs)

            else:
                break

        print("final results")

        print("Final motifs:")
        self.print_motif(new_temp_motifs)
        print("score : ", best_score)
        self.print_motif(current_profile)



random_motif = random_motif_search(11,"input.txt")
random_motif.do_random_motif_search()

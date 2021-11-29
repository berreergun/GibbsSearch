import random

class gibbs_search:

    def __init__(self, k_mers, file_to_read):
        self.k_mers = k_mers
        self.file_to_read = file_to_read

    def CountFrequency(self, my_list):
        # Creating an empty dictionary
        freq = {}

        for item in my_list:
            if (item in freq):
                freq[item] += 1
            else:
                freq[item] = 1

        return freq

    def print_motif(self, motif):
        for r in motif:
            for c in r:
                print(c, end="  ")
            print()

    def score(self, motifs):
        score = 0
        for i in range(len(motifs[0])):
            motif = ''.join([motifs[j][i] for j in range(len(motifs))])
            print(motif, "motif", max(self.CountFrequency(motif).values()), " freq")
            score += 10 - max(self.CountFrequency(motif).values())
        print("Score :", score)
        return score

    def laplace_profile(self, profile, flag, k_mer):
        print("before laplacian profile :")
        self.print_motif(profile)
        if flag == 1:
            for x in range(len(profile)):
                for y in range(k_mer):
                    profile[x][y] = round((int(profile[x][y]) + 1) / 13, 5)

        else:
            for x in range(len(profile)):
                for y in range(10):
                    profile[x][y] = (int(profile[x][y])) / 10

        print("after laplacian profile :")
        self.print_motif(profile)
        return profile

    def profile(self, motifs, k_mer=10):
        # A :
        # C :
        # G :
        # T :
        profile = [[0] * len(motifs[0]) for _ in range(4)]
        flag = 0
        for i in range(len(motifs[0])):
            motif = ''.join([motifs[j][i] for j in range(len(motifs))])
            freq = self.CountFrequency(motif)

            for item in freq:
                if item != "A" or item != "T" or item != "C" or item != "G":
                    flag = 1

            for item in freq:
                if (item == "A"):
                    profile[0][i] = freq[item]
                if (item == "C"):
                    profile[1][i] = freq[item]
                if (item == "G"):
                    profile[2][i] = freq[item]
                if (item == "T"):
                    profile[3][i] = freq[item]
        print("profile :")
        self.print_motif(profile)
        return self.laplace_profile(profile, flag, k_mer)

    def calculate_line_score_with_profile(self, _profile, line):
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

    def print_con(self, current_profile, k_mer):
        # A :
        # C :
        # G :
        # T :
        pattern = ""
        for i in range(k_mer):
            A = current_profile[0][i]
            C = current_profile[1][i]
            G = current_profile[2][i]
            T = current_profile[3][i]

            if A >= C and A >= G and A >= T:
                pattern += 'A'
            elif C >= G and C >= T:
                pattern += 'C'
            elif G >= T:
                pattern += 'G'
            else:
                pattern += 'T'

            pattern = "".join(pattern)

        print("consensus string is ", pattern)

    def do_gibbs_search(self):
        ''' motifs = [['A', 'C', 'A', 'T', 'G', 'A', 'T', 'G', 'T', 'G'],
                  ['A', 'C', 'T', 'T', 'G', 'A', 'T', 'G', 'G', 'A'],
                  ['C', 'G', 'T', 'A', 'G', 'A', 'T', 'T', 'T', 'A'],
                  ['G', 'T', 'A', 'C', 'G', 'C', 'C', 'G', 'C', 'A'],
                  ['A', 'C', 'A', 'G', 'G', 'A', 'C', 'G', 'T', 'G'],
                  ['A', 'C', 'A', 'T', 'G', 'A', 'T', 'G', 'T', 'G'],
                  ['A', 'C', 'T', 'T', 'G', 'A', 'T', 'G', 'G', 'A'],
                  ['C', 'G', 'T', 'A', 'G', 'A', 'T', 'T', 'T', 'A'],
                  ['G', 'T', 'A', 'C', 'G', 'C', 'C', 'G', 'C', 'A'],
                  ['A', 'C', 'A', 'G', 'G', 'A', 'C', 'G', 'T', 'G']]'''

        file = open(self.file_to_read, 'r')
        with open(self.file_to_read) as f:
            motifs = f.readlines()
        motifs = [row.rstrip('\n') for row in motifs]
        self.print_motif(motifs)
        file.close()
        k_number = self.k_mers

        # Algorthm starts : first choose random k mers for every line
        random_motifs = [[0] * k_number for _ in range(10)]
        for i in range(10):
            random_start = random.randint(0, len(motifs[0]) - k_number)
            for k in range(k_number):
                random_motifs[i][k] = motifs[i][k + random_start]
        # now we have our random_motifs
        print("initial motifs")
        self.print_motif(random_motifs)
        # calculate initial score
        best_score = self.score(random_motifs)
        init_score = best_score
        iter = 0
        temp_num = 0
        while 1:
            iter += 1
            print("-----------------iter :", iter, "------------------")
            print("TEMP NUM: ", temp_num)
            # removed one motif
            x = random.randint(0, 9)
            print("removed line")
            print("x: ", x, " ", random_motifs[x])
            random_motifs.pop(x)
            print("removed motif")
            self.print_motif(random_motifs)

            current_profile = self.profile(random_motifs, k_number)
            for l in range(10):
                prob = []
                temp_motifs = []
                end_index = k_number
                for start in range(len(motifs[0]) - k_number):
                    temp_motif = motifs[x][start:end_index]
                    temp_score = self.calculate_line_score_with_profile(current_profile, temp_motif)
                    end_index += 1
                    temp_motifs.append(temp_motif)
                    prob.append(temp_score)

                indices = list(range(0, len(motifs[0]) - k_number))
                index = random.choices(indices, prob)[0]

            random_motifs.insert(x, temp_motifs[index])
            print("new motifs")
            self.print_motif(random_motifs)
            print("new results")
            new_score = self.score(random_motifs)
            if new_score < best_score:
                temp_num = 0
                best_score = new_score
                current_profile = self.profile(random_motifs, k_number)
            else:
                temp_num += 1

            if temp_num == 1000:
                print("final results")
                print("Final motifs:")
                self.print_motif(random_motifs)
                print("initial score : ", init_score)
                print("final score : ", best_score)
                self.print_motif(current_profile)
                self.print_con(current_profile, k_number)
                break

        #return best_score


random_motif = gibbs_search(10, "input.txt")
random_motif.do_gibbs_search()

'''
for i in range(11):
    best = random_motif.do_gibbs_search()
    print(best)
'''


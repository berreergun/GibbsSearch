import random
import numpy as np

def print_motif(motif):
    for r in motif:
        for c in r:
            print(c, end="  ")
        print()

def dna_write_to_file(dna):
    f = open("input.txt", "a")
    for x in dna:
        for y in x:
            f.write(y)
        f.write("\n")

def CountFrequency(my_list):
        # Creating an empty dictionary
        freq = {}

        for item in my_list:
            if (item in freq):
                freq[item] += 1
            else:
                freq[item] = 1

        return freq

def score(motifs):
        score = 0
        for i in range(len(motifs[0])):
            motif = ''.join([motifs[j][i] for j in range(len(motifs))])
            print(motif, "motif", max(CountFrequency(motif).values()), " freq")
            score += 10 - max(CountFrequency(motif).values())
        print("Score :", score)
        return score

def mutuation_generate(line):
    mutuated_line = line.copy()
    random_indexes = random.sample(range(1, 10), 4)
    for index in range(4):
        random_nuk = random.randint(0,3)
        if random_nuk == 0:
            if mutuated_line[random_indexes[index]] != "A":
                mutuated_line[random_indexes[index]] = "A"
            else: mutuated_line[random_indexes[index]] = "G"
        elif random_nuk == 1:
            if mutuated_line[random_indexes[index]] != "C":
                mutuated_line[random_indexes[index]] = "C"
            else:
                mutuated_line[random_indexes[index]] = "A"
        elif random_nuk == 2:
            if mutuated_line[random_indexes[index]] != "G":
                mutuated_line[random_indexes[index]] = "G"
            else:
                mutuated_line[random_indexes[index]] = "T"
        elif random_nuk == 3:
            if mutuated_line[random_indexes[index]] != "T":
                mutuated_line[random_indexes[index]] = "T"
            else:
                mutuated_line[random_indexes[index]] = "C"
    print("line after mutuation:  ", mutuated_line)
    return mutuated_line



dnas = [[0] * 500 for _ in range(10)]

for k in range(10):
    for l in range(500):
        nuk = random.randint(0, 3)
        if nuk == 0:
            dnas[k][l] = "A"
        elif nuk == 1:
            dnas[k][l] = "C"
        elif nuk == 2:
            dnas[k][l] = "G"
        elif nuk == 3:
            dnas[k][l] = "T"


our_motif = ["T","A","G","G","C","T","A","A","A","T"]
mutuatic_motifs = [[0] * 10 for _ in range(10)]


for line in range(10):
    mutation = random.randint(0, 490)
    mutated_line = mutuation_generate(our_motif)

    for start in range(10):
        mutuatic_motifs[line][start] = mutated_line[start]
        dnas[line][mutation+start] = mutated_line[start]


score(mutuatic_motifs)
dna_write_to_file(dnas)



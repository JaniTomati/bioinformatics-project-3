#!/usr/bin/python3
# -*- coding: utf-8 -*-

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Default values
cost = np.array([
    # A  C  G  T
    [0, 5, 2, 5],
    [5, 0, 5, 2],
    [2, 5, 0, 2],
    [5, 2, 5, 0]])
gap_cost = 5
look_up = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 0, "R": 0, "S": 0}


def read_in_sequences(file):
    """ Read in a fasta file containing sequences """
    sequences = []
    for seq_record in SeqIO.parse(file, "fasta"):
        sequences.append(seq_record.seq.upper())

    return sequences


def read_in_scoring_matrix(file):
    """ Read in the alphabet, scoring matrix and gap cost """
    alphabet = []
    scoring_matrix = []

    with open(file, "r") as f:
        lines = f.readlines()

    firstLine = True
    for line in lines:
        line_arr = line.rstrip().split()
        if firstLine == True:
            gapcost = int(line_arr[0])
            firstLine = False
        else:
            row = []
            for char in line_arr:
                if char == line_arr[0]: # add char to alphabet
                    alphabet.append(char)
                else:                   # add cost to scoring matrix
                    row.append(int(char))
            scoring_matrix.append(row)

    return gapcost, np.array(scoring_matrix), alphabet


def levenshtein_distance(seq1, seq2):
    """ Calculate the levenshtein distance between two sequences """
    n, m = len(seq1) + 1, len(seq2) + 1
    D = np.zeros((n, m), dtype=int)

    for i in range(n):
        for j in range(m):
            if i == 0 and j == 0:
                D[i, j] = 0
            elif i == 0:
                D[i, j] = j
            elif j == 0:
                D[i, j] = i
            else:
                indicator = 1
                if seq1[i - 1] == seq2[j - 1]:
                    indicator = 0
                D[i, j] = min(D[i - 1, j - 1] + indicator,
                              D[i, j - 1] + 1,
                              D[i - 1, j] + 1)

    return D[n - 1, m - 1]


def determine_center_sequence(sequences):
    """ Determine the center sequence such that the sum of the levenshtein distance is minimized """
    summed_distances = []
    for seq in sequences:
        sum = 0
        for seq_ in sequences:
            sum += levenshtein_distance(seq, seq_)
        summed_distances.append(sum)

    min = summed_distances[0]
    index = 0
    for i in range(1, len(sequences)):
        if min >= summed_distances[i]:
            min = summed_distances[i]
            index = i

    return index, sequences[index]


def optimal_alignment(seq1, seq2):
    """ Calculate the pairwise optimal alignment between two sequences using linear gapcost """
    m, n = len(seq1) + 1, len(seq2) + 1
    align_matrix = np.zeros((n, m))
    for i in range(n):
        for j in range(m):
            if j == 0 and i == 0:
                align_matrix[i, j] = 0
            elif j == 0:
                align_matrix[i, j] = align_matrix[i-1, j] + gap_cost
            elif i == 0:
                align_matrix[i, j] = align_matrix[i, j-1] + gap_cost

            else:
                val = min(
                    align_matrix[i-1, j-1] + cost[look_up[seq2[i-1]], look_up[seq1[j-1]]],
                    align_matrix[i-1, j] + gap_cost,
                    align_matrix[i, j-1] + gap_cost
                )
                align_matrix[i, j] = val

    return align_matrix


def traceback(opt, seq1, seq2):
    """ Traceback the optimal alignment of two sequences """
    s1 = []
    s2 = []
    i, j = len(seq2), len(seq1)
    while i > 0 and j > 0:
        if opt[i, j] == opt[i-1, j-1]+cost[look_up[seq1[j-1]], look_up[seq2[i-1]]]:
            s2.append(seq2[i-1])
            s1.append(seq1[j-1])
            i -= 1
            j -= 1

        elif opt[i-1, j] + gap_cost == opt[i,j]:
            s1.append("-")
            s2.append(seq2[i-1])
            i -= 1

        elif opt[i, j-1] + gap_cost == opt[i,j]:
            s2.append("-")
            s1.append(seq1[j-1])
            j -= 1

        elif i == 0:
            while j >= 0:
                s2.append("-")
                j -= 1
        elif j == 0:
            while i >= 0:
                s1.append("-")
                i -= 1

    s1 = "".join(s1[::-1])
    s2 = "".join(s2[::-1])

    return s1, s2


def extendExistingRow(current_row, new_gaps):
    """ Introduce new gaps into the existing rows of the matrix M """
    new_row = []
    for i in range(len(current_row)):
        if i in new_gaps:   # insert new gap
            new_row.append("-")
            new_row.append(current_row[i])
        else:
            new_row.append(current_row[i])

    return np.array(new_row)


def extendMSAMatrix(OPT, M):
    """ Extend the multiple alignment matrix M with a new optimal alignment OPT """
    extendedM = []

    i, j = 0, 0
    new_gaps = []
    new_row = [] # row that inserts the new sequence in the matrix M

    while i < len(M[0]):
        if M[0][i] == OPT[0][j]:
            new_row.append(OPT[1][j])
            j += 1
            i += 1
        elif M[0][i] == "-":
            new_row.append("-")
            i += 1
        elif OPT[0][j] == "-":
            new_row.append(OPT[1][j])
            new_gaps.append(i) # introduce new gap in position "i" of the matrix
            j += 1

    for row in M:
        extendedM.append(extendExistingRow(row, new_gaps))

    extendedM.append(np.array(new_row))
    return np.array(extendedM)


def main():
    sequences = ["ATTCT", "ACGT", "CTCGA", "ACGGT"]
    # sequences = read_in_sequences("testdata_short.fasta")
    center_index, center_sequence = determine_center_sequence(sequences)
    print("Determined center sequence:", center_sequence)

    for seq in sequences:
        opt = optimal_alignment(center_sequence, seq)
        alignment = traceback(opt, center_sequence, seq)
        print(alignment)

    # for i in range(len(sequences)):
    #     if i != center_index: # do not calculate optimal alignment between the center sequence and itself
    #         opt = optimal_alignment(sequences[center_index], sequences[i])
    #         alignment = traceback(opt, sequences[center_index], sequences[i])
    #
    #         if i != 0:
    #             M = extendMSAMatrix(alignment, M)
    #         else:   # set M to the first optimal alignment
    #             M = np.array([alignment[0], alignment[1]])
    #     print(i)
    # print(M)



if __name__ == '__main__':
    main()

#!/usr/bin/python3
# -*- coding: utf-8 -*-

import numpy as np

# Default values
cost = np.array([
    # A  C  G  T
    [10, 2, 5, 2],
    [2, 10, 2, 5],
    [5, 2, 10, 2],
    [2, 5, 2, 10]])
gap_cost = -5
look_up = {"A": 0, "C": 1, "G": 2, "T":3}


def levenshtein_distance(seq1, seq2):
    """ Calculate the levenshtein distance between two sequences """
    n = len(seq1) + 1
    m = len(seq2) + 1

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


def main():
    sequences = ["ATTCT", "ACGT", "CTCGA", "ACGGT"]
    center_index, center_sequence = determine_center_sequence(sequences)




if __name__ == '__main__':
    main()

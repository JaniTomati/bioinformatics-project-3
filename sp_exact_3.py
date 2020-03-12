#!/usr/bin/python3
# -*- coding: utf-8 -*-
# python3 sp_exact_3.py -s testdata_short.fasta -c scoring_matrix.txt

import sys
import argparse
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Default values
cost = np.array([
    # A  C  G  T
    [0, 5, 2, 5],
    [5, 0, 5, 2],
    [2, 5, 0, 5],
    [5, 2, 5, 0]])
gap_cost = 5
alphabet = ["A", "C", "G", "T"]
look_up = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 0, "R": 0, "S": 0}


def print_usage():
    """ Define a usage message """
    print("""
          Usage: sp_exact_3.py -s <sequence_file> -c <cost_file>\n
          where <sequence_file> contains a set of sequences over the alphabet
          {A,C,G,T,-} in FASTA format.
          And <cost_file> contains a cost matrix and a gap cost in a phylip-like
          format. If -c is not specified, the default values are used.
          """)


def parse_arguments():
    """ Parse the necessary arguments, otherwise use default """
    global cost, gap_cost, alphabet
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", help = "File that holds all sequences that are to be aligned")
    parser.add_argument("-c", help = "Txt file defining the cost matrix in a Phylip-like format")
    args = parser.parse_args()
    cost_file = args.c
    sequence_file = args.s

    if sequence_file is None:
        print_usage()
        sys.exit(1)

    if cost_file is None:
        print("Using default values:")
        print("Gapcost:", gap_cost)
        print("Subcost:\n", cost)
    else:
        gap_cost, cost, alphabet = read_in_scoring_matrix(cost_file)

    return sequence_file


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


def score_pair(a, b):
    """ Find the cost of the pair """
    if a == "-" and b == "-":
        return 0
    elif a == "-" and b != "-" or a != "-" and b == "-":
        return gap_cost
    elif a != "-" and b != "-":
        return cost[look_up[a], look_up[b]]


def SP(i, j, k):
    """ Calculate the sum of pairs score for three sequences """
    sum_of_pairs = 0
    sum_of_pairs += score_pair(i, j)
    sum_of_pairs += score_pair(i, k)
    sum_of_pairs += score_pair(j, k)
    return sum_of_pairs


def compute_msa(A, B, C):
    """ Iteratively compute the MSA between three sequences A, B, C """
    n, m, o = len(A) + 1, len(B) + 1, len(C) + 1
    T = np.zeros((n, m, o), dtype=int)

    for i in range(n):
        for j in range(m):
            for k in range(o):
                v0, v1, v2, v3, v4, v5, v6, v7 = sys.maxsize, sys.maxsize, sys.maxsize, sys.maxsize, sys.maxsize, sys.maxsize, sys.maxsize, sys.maxsize
                if i == 0 and j == 0 and k == 0:
                    v0 = 0
                if i > 0 and j > 0 and k > 0:
                    v1 = T[i-1][j-1][k-1] + SP(A[i - 1], B[j - 1], C[k - 1])
                if i > 0 and j > 0 and k >= 0:
                    v2 = T[i-1][j-1][k] + SP(A[i - 1], B[j - 1], "-")
                if i > 0 and j >= 0 and k > 0:
                    v3 = T[i-1][j][k-1] + SP(A[i - 1], "-", C[k - 1])
                if i >= 0 and j > 0 and k > 0:
                    v4 = T[i][j-1][k-1] + SP("-", B[j - 1], C[k - 1])
                if i > 0 and j >= 0 and k >= 0:
                    v5 = T[i-1][j][k] + SP(A[i - 1], "-", "-")
                if i >= 0 and j > 0 and k >= 0:
                    v6 = T[i][j-1][k] + SP("-", B[j - 1], "-")
                if i >= 0 and j >= 0 and k > 0:
                    v7 = T[i][j][k-1] + SP("-", "-", C[k - 1])
                T[i][j][k] = min(v0, v1, v2, v3, v4, v5, v6, v7)

    return T


def traceback(T, seq1, seq2, dim):
    """ Traceback the optimal alignment of two sequences """
    s1 = ""
    s2 = ""
    i, j = len(seq1), len(seq2)
    while i >= 0 and j >= 0:
        if opt[dim][i, j] == opt[dim][i-1, j-1] + cost[look_up[seq1[i-1]], look_up[seq2[j-1]]]:
            s1 += seq1[i-1]
            s2 += seq2[j-1]
            i -= 1
            j -= 1
        elif opt[dim][i,j] == opt[dim][i-1, j] + gap_cost:
            s1 += seq1[i-1]
            s2 += "-"
            i -= 1
        elif opt[dim][i,j] == opt[dim][i, j-1] + gap_cost:
            s1 += "-"
            s2 += seq2[j-1]
            j -= 1
        elif i == 0:
            while j >= 0:
                if j > 0:
                    s1 += ("-")
                j -= 1
        elif j == 0:
            while i >= 0:
                if i > 0:
                    s2 += ("-")
                i -= 1

    return s1[::-1], s2[::-1]


def backtracking(T, A, B, C):
    """ Create a traceback of the multiple sequence alignment """
    s1, s2, s3 = "", "", ""
    i, j, k = len(A), len(B), len(C)
    while i >= 0 and j >= 0 and k >= 0:
        if T[i][j][k] == T[i - 1][j - 1][k - 1] + SP(A[i - 1], B[j - 1], C[k - 1]):
            s1 += A[i-1]
            s2 += B[j-1]
            s3 += C[k-1]
            i -= 1
            j -= 1
            k -= 1
        elif T[i][j][k] == T[i-1][j-1][k] + SP(A[i - 1], B[j - 1], "-"):
            s1 += A[i-1]
            s2 += B[i-1]
            s3 += "-"
            i -= 1
            j -= 1
        elif T[i][j][k] == T[i-1][j][k-1] + SP(A[i - 1], "-", C[k - 1]):
            s1 += A[i-1]
            s2 += "-"
            s3 += C[k-1]
            i -= 1
            k -= 1
        elif T[i][j][k] == T[i][j-1][k-1] + SP("-", B[j - 1], C[k - 1]):
            s1 += "-"
            s2 += B[j-1]
            s3 += C[k-1]
            j -= 1
            k -= 1
        elif T[i][j][k] == T[i-1][j][k] + SP(A[i - 1], "-", "-"):
            s1 += A[i-1]
            s2 += "-"
            s3 += "-"
            i -= 1
        elif T[i][j][k] == T[i][j-1][k] + SP("-", B[j - 1], "-"):
            s1 += "-"
            s2 += B[j-1]
            s3 += "-"
            j -= 1
        elif T[i][j][k] == T[i][j][k-1] + SP("-", "-", C[k - 1]):
            s1 += "-"
            s2 += "-"
            s3 += C[k - 1]
            k -= 1
        elif j == 0 and k == 0:
            while i >= 0:
                if i > 0:
                    s2 += ("-")
                    S3 += ("-")
                i -= 1
        elif i == 0 and k == 0:
            while j >= 0:
                if j > 0:
                    s1 += ("-")
                    S3 += ("-")
                j -= 1
        elif i == 0 and j == 0:
            while k >= 0:
                if k > 0:
                    s1 += ("-")
                    S2 += ("-")
                k -= 1

    # print(s1[::-1] + "\n" + s2[::-1] + "\n" + s3[::-1])
    return (s1[::-1], s2[::-1], s3[::-1])


def main():
    sequence_file = parse_arguments()
    print("Reading sequences from", sequence_file, "\n")

    # sequences = ["ATTCT", "ACGT", "CTCGA", "ACGGT"]
    sequences = read_in_sequences(sequence_file)
    if len(sequences) > 3:
        print("Only using the first three sequences from the file.")

    T = compute_msa(sequences[0], sequences[1], sequences[2])
    alignment = backtracking(T, sequences[0], sequences[1], sequences[2])

    print("\nMSA")
    print("----------------------------------------")
    print(alignment[0] + "\n" + alignment[1] + "\n" + alignment[2])


if __name__ == '__main__':
    main()

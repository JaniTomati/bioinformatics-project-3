#!/usr/bin/python3
# -*- coding: utf-8 -*-

import numpy as np

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


def extendM(OPT, M):
    """ Extend the multiple alignment matrix M with a new optimal alignment OPT """
    extendedM = []

    j = 0
    i = 0
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
    M = np.array([
        ["A", "-", "-", "C", "G", "T"],
        ["A", "T", "T", "C", "-", "T"],
        ["C", "T", "-", "C", "G", "A"]])

    OPT = ["ACG-T", "ACGGT"]

    M = extendM(OPT, M)
    print(M)


if __name__ == '__main__':
    main()

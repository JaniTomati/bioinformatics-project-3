#!/usr/bin/python3
# -*- coding: utf-8 -*-

import time
import argparse
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Example call
# python3 sp_approx.py --C scoring_matrix.txt --P Y Seqs.fasta

# Default values
cost = np.array([
		[0, 5, 2, 5, 5],  # A
		[5, 0, 5, 2, 5],  # C
		[2, 5, 0, 5, 5],  # G
		[5, 2, 5, 0, 5],  # T
		[5, 5, 5, 5, 0]])  #-'

look_up = {'a':0, 'c':1, 'g':2, 't':3, 'A':0, 'C':1, 'G':2, 'T':3, '-':4, 'N':0, 'R':0, 'S':0}

gap_cost = 5



# Load scoring matrix
def read_in_scoring_matrix():
	""" Read in the alphabet, scoring matrix and gap cost """
	alphabet_size = None
	scoring_matrix = []
	alphabet = []

	with open(file, "r") as f:
		lines = f.readlines()

	for line in lines:
		if line == lines[0]:
			alphabet_size = int(line.rstrip())
		else:
			current_line = line.rstrip().replace(" ", "")
			row = []
			for char in current_line:
				if char == current_line[0]:
					alphabet.append(char)
				else:
					row.append(int(char))
			scoring_matrix.append(row)

	return (alphabet_size, alphabet, np.array(scoring_matrix))



# Optimal linear alignment of 2 sequences
def align_linear(seq1, seq2):
	n=len(seq1)+1
	m=len(seq2)+1
	S= np.zeros((n,m))
	for j in range(0,m):
		for i in range(0,n):
			if i==0 and j==0:
				S[i,j]=0
			elif i==0:
				S[i,j]=S[i,j-1]+gap_cost
			elif j==0:
				S[i,j]=S[i-1,j]+gap_cost
			else:
				S[i,j]=min(S[i-1,j]+gap_cost,
						   S[i,j-1]+gap_cost,
						   S[i-1,j-1]+cost[look_up[seq1[i-1]], look_up[seq2[j-1]]])

	return S

# Traceback algorithm - pairwise linear alignment
def backtrack(matrix, seq1, seq2):
	i=len(seq1)
	j=len(seq2)
	align1 = ""
	align2 = ""
	while i>0 or j>0:
		if matrix[i,j] == matrix[i-1,j-1]+cost[look_up[seq1[i-1]], look_up[seq2[j-1]]]:
			align1 = str(seq1[i-1]) + align1
			align2 = str(seq2[j-1]) + align2
			i-=1
			j-=1
		elif matrix[i,j] == matrix[i-1,j]+gap_cost:
			align1 = str(seq1[i-1]) + align1
			align2 = "-" + align2
			i -=1
		else:
			align2 = str(seq2[j-1]) + align2
			align1 = "-" + align1
			j -=1

	return [align1, align2]

## Merging of alignments
def merge_alignment(seq1, seq2, alignment):
	num1 = []
	num2 = []
	a = len(alignment)
	alignment.append([])
	for i in seq1[0]:
		num1.append(str(i))
	for j in seq2[0]:
		num2.append(str(j))
	i = 0
	j = 0
	k = 0
	while k < len(alignment[0]) or i < len(num1):
		if k >= len(alignment[0]):
			alignment[-1].append(num2[j])
			for b in range(0,a):
				alignment[b].append("-")
			i += 1
			j += 1
			k += 1
		elif i >= len(num1):
			alignment[-1].append("-")
			k +=1
		elif num1[i] == alignment[0][k]:
			alignment[-1].append(num2[j])
			i += 1
			j += 1
			k += 1
		elif num1[i] == "-":
			alignment[-1].append(num2[j])
			for b in range(0,a):
				alignment[b].insert(k,"-")
			i +=1
			j +=1
			k +=1
		elif alignment[0][k] == "-":
			alignment[-1].append("-")
			k +=1
	return alignment

## Find center string
def find_center(seqs):
	n = len(seqs)
	S= np.zeros((n,n))
	for i in range(0,n):
		for j in range(0,n):
			a = align_linear(seqs[i],seqs[j])[-1,-1]
			S[i,j] = a
	sum_val = []
	sum_val.append(np.sum(S,axis=1))
	return np.argmin(sum_val)

## Combining the different functions
def combine(list_of_seqs):
	center = find_center(list_of_seqs)
	list_of_seqs[0], list_of_seqs[center] = list_of_seqs[center], list_of_seqs[0]
	alignment = [[]]
	stored_alignments = []
	for i in list_of_seqs[0]:
		alignment[0].append(i)
	for i in range(1,len(list_of_seqs)):
		matrix = align_linear(list_of_seqs[0],list_of_seqs[i])
		pair_align = backtrack(matrix,list_of_seqs[0],list_of_seqs[i])
		stored_alignments.append(pair_align)
		new = merge_alignment([pair_align[0]], [pair_align[1]],alignment)
		alignment = new
	return alignment, stored_alignments

## Find the score of the alignment
def compute_sp_score(alignment):
	score = 0
	for i in range(len(alignment)):
		for j in range(i+1, len(alignment)):
			for c in range(len(alignment[i])):
				score = score +cost[look_up[alignment[i][c]], look_up[alignment[j][c]]]
	return score


def measure_times():
	""" Use this method to measure times of the algorithms using sequences of different length """
	sequences = []
	sequence_length = [5, 10, 50, 100, 500, 1000, 2500, 5000, 7500, 10000] # edit sequence

	for seq_record in SeqIO.parse("sequences/rand_sequences.fasta", "fasta"):
		sequences.append(seq_record.seq.upper())

	times_opt = []
	times_back = []
	i = 0
	while i < len(sequences):
		# measure optimal alignment algorithm
		start_opt = time.time()
		opt = optimal_alignment(sequences[i], sequences[i + 1])
		end_opt = time.time()
		times_opt.append(end_opt - start_opt)

		# measure trace back algorithm
		start_back = time.time()
		aligned = traceback(opt, sequences[i], sequences[i + 1])
		end_back = time.time()
		times_back.append(end_back - start_back)
		i += 2

	return times_opt, times_back


def test_alignment_algorithm():
	""" Test the alignment algorithm by using the sequences from the project2_examples.txt """
	sequences = []

	for seq_record in SeqIO.parse("sequences/test_sequences.fasta", "fasta"):
		sequences.append(seq_record.seq)

	i = 0
	while i < len(sequences):
		opt = optimal_alignment(sequences[i], sequences[i + 1])
		print("\nScore of the optimal global optimal_alignment: ", opt[len(sequences[i + 1]), len(sequences[i])])

		aligned = traceback(opt, sequences[i], sequences[i + 1])
		print(aligned[0], "\n" + aligned[1], "\n")
		i += 2


def verify_MSA(M, alignments):
    """ Get the pairwise alignments from the MSA matrix by deleting gap columns.
        Compare the alignments to the optimal alignments that have been calculated before """
    s1 = ""
    s2 = ""
    count_alignments = 0
    verified = True
    for i in range(1, len(M)):
        for j in range(len(M[i])):
            if not (M[0][j] == "-" and M[i][j] == "-"): # M[0] is center sequence
                s1 += M[0][j]
                s2 += M[i][j]
        if alignments[count_alignments][0] != s1:
            verified = False
            print("\n" + alignments[count_alignments][0] + "\n" + s1 + "\n")
        if alignments[count_alignments][1] != s2:
            verified = False
            print(alignments[count_alignments][1] + "\n" + s2 + "\n")

        count_alignments += 1
        s1 = ""
        s2 = ""

    return verified


def main():
	global gap_cost, cost, file
	# Arg parse
	parser = argparse.ArgumentParser()
	parser.add_argument("Seqs", help = "Fasta file path for sequences")
	parser.add_argument("--P", help = "Output strings of aligned sequences: Write Y")
	parser.add_argument("--C", help = "txt file with Phylip-like format")
	parser.add_argument("--G", help = "Gap cost")
	args = parser.parse_args()
	seqs = args.Seqs
	if args.G != None:
		gap_cost = np.array([args.G])
		gap_cost = gap_cost.astype("int32")[0]

	print_statement = args.P

	file = args.C


	if print_statement is not None:
		print_statement = print_statement.lower() # lower case letters

	if file is not None:
		alphabet_size, alphabet, cost = read_in_scoring_matrix()
		look_up = dict()
		for i in range(alphabet_size):
			look_up[alphabet[i].upper()] = i

	### Running the Code ###
	seq_list = []
	count_records = 0
	for seq_record in SeqIO.parse(seqs, "fasta"):
		if seq_record.seq != "":
			seq_list.append((seq_record.seq.upper()))
		else:
			print("Info: Skipped sequence", count_records, "because it was empty.\n")
		count_records += 1

	start_opt = time.time()
	opt, alignments = combine(seq_list)
	end_opt = time.time()

	verify = verify_MSA(opt, alignments)
	if verify:
		print("MSA verified.")
	else:
		print("MSA not verified.")

	op1 = []
	for st in opt:
		coll = ""
		for base in st:
			coll = coll + base
		op1.append(coll)
	print("\nAlignment: ")
	for s in op1:
		print(s)
	print("Making the optimal alignment took", end_opt - start_opt, "seconds.\n")

	start_back = time.time()
	aligned = compute_sp_score(opt)
	end_back = time.time()

	if print_statement == "y":
		print("The score of the alignment: ", aligned)
		print("Calculating the score", end_back - start_back, "seconds.\n")

	### Writing the aligned sequences into a fasta file ###
	#record = SeqRecord(Seq(opt),id="Seq1")

	#result_file = "result_linear.fasta"
	#SeqIO.write([record], result_file, "fasta")
	#print("Wrote results in fasta format to " +  result_file + ".")

	# test alignment algorithm using examples
	#test_alignment_algorithm()

	# function for taking the performance measures
	# times_opt, times_back = measure_times()
	# print("Optimal alignment:", times_opt)
	# print("Traceback", times_back)


if __name__ == '__main__':
	main()

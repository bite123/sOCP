# _*_ coding: UTF-8 _*_
# Author: PENG Zhao
# Version: 2023-10
""" Instruction
This script provides functions extracting features from an extended ORF.
""" 

import sys
import itertools
import Util

NTCHAR = ["A", "C", "G", "T"]

################################################################################
# Class Ntseq
################################################################################

class Ntseq:
	""" Refer to Seq.py from Biopython.
	A string with biological methods for feature extraction.
	"""

	def __init__(self, data):
		""" Create a Ntseq object. 
		>>> my_seq = Ntseq("GAGATGGCGGCGGAGCGGGGAGCCCGGCGACTCCTCAGC")
		"""
		if not isinstance(data, str):
			raise TypeError(
				"The sequence data should be a string."
			)
		self._data = data

	def __str__(self):
		# Return the full sequence as a python string, use str(my_seq).
		return self._data

	def length(self):
		# Return the length of the sequence, use len(my_seq).
		return len(self._data)

	def count(self, sub, start=0, end=sys.maxsize):
		""" Return a non-overlapping count, like that of a python string.
		>>> print(Ntseq("AAAA").count("AA"))
		2
		"""
		return str(self).count(str(sub), start, end)

	def count_overlap(self, sub, start=0, end=sys.maxsize):
		""" Return an overlapping count.
		>>> print(Ntseq("AAAA").count_overlap("AA"))
		3
		"""
		sub_str = str(sub)
		self_str = str(self)
		overlap_count = 0
		while True:
			start = self_str.find(sub_str, start, end) + 1
			if start != 0:
				overlap_count += 1
			else:
				return overlap_count

	def frame_split(self, k=3):
		# Return a list of codons, or a list of hexamer (k=6).
		res_list = []
		seq = str(self)
		codon_number = int(self.length()/k)
		for i in range(codon_number):
			res_list.append(seq[i*k:i*k+k])
		return res_list

	def translate(self):
		# Return the tranlated amino acid sequence for ORF.
		aa_list = []
		codon_list = self.frame_split()
		for codon in codon_list:
			aa = Util.CODON_TABLE[codon]
			aa_list.append(aa)
		return "".join(aa_list)

	def k_mer(self, k):
		# Return k-mer of the sequence.
		res_dict = {}
		for nt in itertools.product(NTCHAR, repeat=k):
			nt = "".join(nt)
			count = self.count_overlap(nt)
			total = self.length() - k + 1
			ratio = round(count/total, 8)
			res_dict[nt] = ratio
		return res_dict

	def g_gap(self, g):
		# Return g-bigap, the circumstance of two nucleotides with g gap-nucleotides.
		res_dict = {}
		total = self.length() - g - 1
		for nt in itertools.product(NTCHAR, repeat=2):
			nt = "".join(nt)
			res_dict[nt] = 0
		seq = str(self)
		for i in range(total):
			tar_str = seq[i] + seq[i+g+1]
			res_dict[tar_str] += 1
		for k,v in res_dict.items():
			v = round(v/total, 8)
			res_dict[k] = v
		return res_dict

	def g_bigap(self, g):
		# Return g-bigap, the circumstance of two DI-nucleotides with g gap-nucleotides.
		res_dict = {}
		total = self.length() - g - 3
		for nt in itertools.product(NTCHAR, repeat=4):
			nt = "".join(nt)
			res_dict[nt] = 0
		seq = str(self)
		for i in range(total):
			tar_str = seq[i:i+2] + seq[i+g+2:i+g+4]
			res_dict[tar_str] += 1
		for k,v in res_dict.items():
			v = round(v/total, 8)
			res_dict[k] = v
		return res_dict

	def fickett(self):
		# Return a parameter for each type of base, originated from fickett score, refer to DeepCPP.
		res_dict = {}
		codon_list = self.frame_split()
		for nt in NTCHAR:
			nt1, nt2, nt3 = 0, 0, 0
			for codon in codon_list:
				if codon[0] == nt:
					nt1 += 1
				if codon[1] == nt:
					nt2 += 1
				if codon[2] == nt:
					nt3 += 1
			f = max(nt1, nt2, nt3)/(min(nt1, nt2, nt3) + 1)
			res_dict[nt] = round(f, 8)
		return res_dict

	# The following features are dependent on collected data: codon bias, hexamer, and nucleotide bias.
	def codon_data(self):
		# Return codon numbers for data collection.
		res_dict = {}
		codon_list = self.frame_split()
		for codon in codon_list:
			if codon not in res_dict:
				res_dict[codon] = 1
			else:
				res_dict[codon] += 1
		return res_dict

	def codon_bias(self, codon_dict):
		# Return codon bias score.
		score = 0
		codon_list = self.frame_split()
		for codon in codon_list:
			score += float(codon_dict[codon])
		return round(score/len(codon_list), 8)

	def hexamer_data(self):
		# Return hexamer numbers for data collection.
		res_dict = {}
		hexa_list = []
		codon_list = self.frame_split(3)
		for i in range(len(codon_list) - 1):
			hexa = codon_list[i] + codon_list[i+1]
			hexa_list.append(hexa)
		for hexa in hexa_list:
			if hexa not in res_dict:
				res_dict[hexa] = 1
			else:
				res_dict[hexa] += 1
		return res_dict

	def hexamer_score(self, hexa_dict):
		# Return hexamer score.
		score = 0
		hexa_list = []
		codon_list = self.frame_split(3)
		for i in range(len(codon_list) - 1):
			hexa = codon_list[i] + codon_list[i+1]
			hexa_list.append(hexa)
		for hexa in hexa_list:
			score += float(hexa_dict[hexa])
		return round(score/len(hexa_list), 8)

	def nucleotide_bias_data(self):
		"""Return nucleotide bias for data collection and score calculation.
		Actually bases at -3, -2, -1, 4, 5, 6, Refer to DeepCPP.
		ATTENTIION!!!
		Unlike other methods in class Ntseq,
		this method requires -3, -2, and -1 site flanking ORF! 
		Which means the sequence treated in this method is extended ORF!
		This method also accepts base N!
		ATTENTIION!!!
		"""
		seq = str(self)
		res_list = [seq[0], seq[1], seq[2], seq[6], seq[7], seq[8]]
		return res_list

	def nucleotide_bias(self, ntbias_dict):
		"""Return nucleotide bias.
		ATTENTION!!!
		Unlike other methods in class Ntseq,
		this method requires -3, -2, and -1 site flanking ORF! 
		Which means the sequence treated in this method is extended ORF!
		This method also accepts base N and other non-ACGT Nucleotide!
		ATTENTION!!!
		"""
		score = []
		n_list = self.nucleotide_bias_data()
		for i in range(6):
			if n_list[i] not in NTCHAR:
				score.append(0.0)
				continue
			key = str(i) + "_" + n_list[i]
			score.append(round(float(ntbias_dict[key]), 8))
		return score


################################################################################
# Processing Function
################################################################################

# This function accepts an extended ORF and a list of feature IDs, 
# and extracts features from the ORF according to the ID list.
# The three dictionaries are required for data dependent features.
# The feature IDs are described in Data/Feature.list.
def feature_extraction(orf, codon_dict, hexa_dict, ntbias_dict, feature_list=[]):
	ntseq = Ntseq(orf[3:])
	ext_ntseq = Ntseq(orf)
	if not feature_list:
		feature_list = list(range(0, 1170))
	# initialization
	nt_fick_list = []
	ntbias_feature = []
	k1mer = []
	k2mer = []
	k3mer = []
	k4mer = []
	g1gap = []
	g2gap = []
	g3gap = []
	g1big = []
	g2big = []
	g3big = []
	# parsing features
	final_feature = []
	for i in feature_list:
		if i == 1:
			final_feature.append(ntseq.length())
		elif i in list(range(2,6)):
			if not nt_fick_list:
				nt_fick = ntseq.fickett()
				nt_fick_list = [nt_fick[k] for k in sorted(nt_fick.keys())]
			new_i = i - 2
			final_feature.append(nt_fick_list[new_i])
		elif i == 6:
			hexa_feature = ntseq.hexamer_score(hexa_dict)
			final_feature.append(hexa_feature)
		elif i == 7:
			codon_feature = ntseq.codon_bias(codon_dict)
			final_feature.append(codon_feature)
		elif i in list(range(8,14)):
			if not ntbias_feature:
				ntbias_feature = ext_ntseq.nucleotide_bias(ntbias_dict)
			new_i = i - 8
			final_feature.append(ntbias_feature[new_i])
		elif i in list(range(14,18)):
			if not k1mer:
				nt_mer = ntseq.k_mer(1)
				k1mer = [nt_mer[k] for k in sorted(nt_mer.keys())]
			new_i = i - 14
			final_feature.append(k1mer[new_i])
		elif i in list(range(18,34)):
			if not k2mer:
				nt_mer = ntseq.k_mer(2)
				k2mer = [nt_mer[k] for k in sorted(nt_mer.keys())]
			new_i = i - 18
			final_feature.append(k2mer[new_i])
		elif i in list(range(34,98)):
			if not k3mer:
				nt_mer = ntseq.k_mer(3)
				k3mer = [nt_mer[k] for k in sorted(nt_mer.keys())]
			new_i = i - 34
			final_feature.append(k3mer[new_i])		
		elif i in list(range(98,354)):
			if not k4mer:
				nt_mer = ntseq.k_mer(4)
				k4mer = [nt_mer[k] for k in sorted(nt_mer.keys())]
			new_i = i - 98
			final_feature.append(k4mer[new_i])
		elif i in list(range(354,370)):
			if not g1gap:
				nt_gap = ntseq.g_gap(1)
				g1gap = [nt_gap[k] for k in sorted(nt_gap.keys())]
			new_i = i - 354
			final_feature.append(g1gap[new_i])
		elif i in list(range(370,386)):
			if not g2gap:
				nt_gap = ntseq.g_gap(2)
				g2gap = [nt_gap[k] for k in sorted(nt_gap.keys())]
			new_i = i - 370
			final_feature.append(g2gap[new_i])
		elif i in list(range(386,402)):
			if not g3gap:
				nt_gap = ntseq.g_gap(3)
				g3gap = [nt_gap[k] for k in sorted(nt_gap.keys())]
			new_i = i - 386
			final_feature.append(g3gap[new_i])
		elif i in list(range(402,658)):
			if not g1big:
				nt_bigap = ntseq.g_bigap(1)
				g1big = [nt_bigap[k] for k in sorted(nt_bigap.keys())]
			new_i = i - 402
			final_feature.append(g1big[new_i])
		elif i in list(range(658,914)):
			if not g2big:
				nt_bigap = ntseq.g_bigap(2)
				g2big = [nt_bigap[k] for k in sorted(nt_bigap.keys())]
			new_i = i - 658
			final_feature.append(g2big[new_i])
		elif i in list(range(914,1170)):
			if not g3big:
				nt_bigap = ntseq.g_bigap(3)
				g3big = [nt_bigap[k] for k in sorted(nt_bigap.keys())]
			new_i = i - 914
			final_feature.append(g3big[new_i])
		else:
			raise TypeError(
				"A wrong number in feature list: " + str(i)
			)
	return final_feature


# Originated from feature_extraction().
# This function accepts a fasta file of extended ORFs and a list of feature IDs, 
# and generate a file saving features of ORFs according to the ID list.
def feature_extraction_batch(input_file, output_file, codon_dict, hexa_dict, ntbias_dict, feature_list=[]):
	seq_dict = Util.fasta_to_dict(input_file, parse=True, sep="")
	with open(output_file, "w") as f_out:
		for k,v in seq_dict.items():
			orf_id = k
			orf_feature = feature_extraction(v, codon_dict, hexa_dict, ntbias_dict, feature_list)
			new_line = orf_id + "\t" + "\t".join([str(i) for i in orf_feature]) + "\n"
			f_out.write(new_line)



# _*_ coding: UTF-8 _*_
# Author: PENG Zhao
# Version: 2023-10
""" Instruction
This script provides functions processing data dependent features,
i.e. codon bias, hexamer, and nucleotide bias.
""" 

import argparse
import math
from FeatureExtraction import Ntseq
import Util

HEXA_LIST = Util.load_hexamer_list()


# This function accepts a fasta file, and generates data files.
# The data files are used to generate data dictionaries for data dependent features. 
def data_collection(input_file, output_prefix, parse=True, sep=""):
	seq_dict = Util.fasta_to_dict(input_file, parse=parse, sep=sep)
	codon_file = output_prefix + ".codon_data.txt"
	hexa_file = output_prefix + ".hexa_data.txt"
	ntbias_file = output_prefix + ".ntbias_data.txt"

	# generating codon bias file
	dict_list = []
	for value in seq_dict.values():
		ntseq = Ntseq(value[3:]) # for ORF
		nt_dict = ntseq.codon_data()
		dict_list.append(nt_dict)
	sum_dict = Util.dict_sum(dict_list)
	Util.dict_to_tsv(sum_dict, codon_file)

	# generating hexamer file
	dict_list = []
	for value in seq_dict.values():
		ntseq = Ntseq(value[3:]) # for ORF
		nt_dict = ntseq.hexamer_data()
		dict_list.append(nt_dict)
	sum_dict = Util.dict_sum(dict_list)
	Util.dict_to_tsv(sum_dict, hexa_file)

	# generating nucleotide bias file
	fi_dict = {}
	for i in range(6):
		fi_dict[i] = {"A":0, "T":0, "G":0, "C":0, "N":0}
	for value in seq_dict.values():
		ntseq = Ntseq(value[0:]) # for extended ORF
		ntbias = ntseq.nucleotide_bias_data()
		for i in range(len(ntbias)):
			nt = ntbias[i]
			fi_dict[i][nt] += 1
	Util.dict_to_tsv(fi_dict, ntbias_file)


################################################################################
# The following four functions are applied in data_dict_building()
################################################################################
def dict_merged(dict_a, dict_b):
	res_dict = {}
	for key in sorted(dict_a.keys()):
		value_a = dict_a[key]
		value_b = dict_b[key]
		ratio = round(math.log(value_a/value_b), 8)
		res_dict[key] = [str(value_a), str(value_b), str(ratio)]
	return res_dict

def hexa_dict(hexa_file):
	with open(hexa_file) as f_in:
		count_dict = {}
		total = 0
		for line in f_in:
			elements = line.strip().split("\t")
			count_dict[elements[0]] = int(elements[1])
			total += int(elements[1])
		for hexa in HEXA_LIST:
			if hexa not in count_dict:
				count_dict[hexa] = 1
				total += 1
		res_dict = {}
		for key,value in count_dict.items():
			res_dict[key] = round(value/total, 8)
		return res_dict

def codon_dict(codon_file):
	with open(codon_file) as f_in:
		count_dict = {}
		for line in f_in:
			elements = line.strip().split("\t")
			count_dict[elements[0]] = int(elements[1])
		aa_dict = {}
		for key,value in Util.CODON_TABLE.items():
			if key not in count_dict:
				count_dict[key] = 1
			if value in aa_dict:
				aa_dict[value] += count_dict[key]
			else:
				aa_dict[value] = count_dict[key]
		res_dict = {}
		for key,value in count_dict.items():
			aa = Util.CODON_TABLE[key]
			ratio = round(value/aa_dict[aa], 8)
			res_dict[key] = ratio
		return res_dict

def ntbias_dict(ntbias_file):
	with open(ntbias_file) as f_in:
		site_dict = {}
		for line in f_in:
			elements = line.strip().split("\t")
			site_dict[elements[0]] = eval(elements[1])
		res_dict = {}
		for k,v in site_dict.items():
			total = 0
			for nt,count in v.items():
				if nt != "N":
					total += int(count)
			for nt,count in v.items():
				if nt != "N":
					new_key = k + "_" + nt
					ratio = round(count/total, 8)
				res_dict[new_key] = ratio
		return res_dict
################################################################################
# The above four functions are applied in data_dict_building()
################################################################################


# This function accepts data files from both positive and negative dataset,
# and generate files saving data dictionaries used for extracting data dependent features
def data_dict_building(pos_prefix, neg_prefix, output_prefix):
	pos_codon_file = pos_prefix + ".codon_data.txt"
	pos_hexa_file = pos_prefix + ".hexa_data.txt"
	pos_ntbias_file = pos_prefix + ".ntbias_data.txt"
	neg_codon_file = neg_prefix + ".codon_data.txt"
	neg_hexa_file = neg_prefix + ".hexa_data.txt"
	neg_ntbias_file = neg_prefix + ".ntbias_data.txt"
	codon_out = output_prefix + ".codonbias.tsv"
	hexa_out = output_prefix + ".hexamer.tsv"
	ntbias_out = output_prefix + ".ntbias.tsv"

	hexa_a = hexa_dict(pos_hexa_file)
	hexa_b = hexa_dict(neg_hexa_file)
	fi_hexa = dict_merged(hexa_a, hexa_b)
	Util.dict_to_tsv(fi_hexa, hexa_out, sep="\t")
	 
	codon_a = codon_dict(pos_codon_file)
	codon_b = codon_dict(neg_codon_file)
	fi_codon = dict_merged(codon_a, codon_b)
	Util.dict_to_tsv(fi_codon, codon_out, sep="\t")

	ntbias_a = ntbias_dict(pos_ntbias_file)
	ntbias_b = ntbias_dict(neg_ntbias_file)
	fi_ntbias = dict_merged(ntbias_a, ntbias_b)
	Util.dict_to_tsv(fi_ntbias, ntbias_out, sep="\t")


# This function is applied to load dictionaries building in data_dict_building().
def data_dict_read(input_file):
	res_dict = {}
	with open(input_file) as f_in:
		for line in f_in:
			elements = line.strip().split("\t")
			res_dict[elements[0]] = elements[3]
	return res_dict


################################################################################
# Main function
################################################################################
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Preparing information for data-based-features")
	parser.add_argument("-p", "--positive",
						help="Input positive dataset fasta", required=True)
	parser.add_argument("-n", "--negative",
						help="Input negative dataset fasta", required=True)
	parser.add_argument("-o", "--output",
						help="Output prefix", required=True)
	ARGS = parser.parse_args()

	data_collection(ARGS.positive, "tmp_pos", parse=True, sep="")
	data_collection(ARGS.negative, "tmp_neg", parse=True, sep="")
	data_dict_building("tmp_pos", "tmp_neg", ARGS.output)
# _*_ coding: UTF-8 _*_
# Author: PENG Zhao
# Version: 2023-10
""" Instruction
This script processes positive and negative dataset,
to get traning, testing, cross validation datasets.
"""

import argparse
import random
import Util
from FeatureExtraction import Ntseq

random.seed(1234)

################################################################################
# Sampling function
################################################################################

# This function accepts a dictionary, and randomly selects a number of keys.
def dict_select(seq_dict, sample_count):
	sample_dict = {}
	sample_list = random.sample(sorted(seq_dict), sample_count)
	for k,v in seq_dict.items():
		if k in sample_list:
			sample_dict[k] = v
	return sample_dict

# This function accepts a dictionary, and randomly divides it at ratio.
def dict_divide(seq_dict, ratio):
	sample_dict = {}
	nonsp_dict = {}
	sample_count = int(len(seq_dict) * ratio)
	sample_list = random.sample(sorted(seq_dict), sample_count)
	for k,v in seq_dict.items():
		if k in sample_list:
			sample_dict[k] = v
		else:
			nonsp_dict[k] = v
	return [sample_dict, nonsp_dict]

################################################################################
# Filter sequences in dictionary
################################################################################

# This function accepts a directory with sequences as values,
# and remove those key-values with sequences:
# 1) containing non-ACGT in coding region
# 2) containing an early stop codon
# 3) containing no stop codon
# 3) less than 11 codons (including stop codon)
def dict_filter(seq_dict):
	new_dict = {}
	del_dict = {}
	for key,value in seq_dict.items():
		ntseq = Ntseq(value[3:])
		if len(value)%3 != 0 or value[-3:] not in ["TAG", "TAA", "TGA"]:
			del_key = ">Our_of_frame#" + key
			del_dict[del_key] = value
		elif ntseq.count("A") + ntseq.count("C") + ntseq.count("G") + ntseq.count("T") != ntseq.length():
			del_key = ">N_include#" + key
			del_dict[del_key] = value
		else:
			aa = ntseq.translate()
			test_seq = aa[:-1]
			if "*" in test_seq:
				del_key = ">Stop_codon#" + key
				del_dict[del_key] = value
			elif len(aa) < 11:
				del_key = ">Too_short#" + key
				del_dict[del_key] = value
			else:
				new_dict[key] = value
	return [new_dict, del_dict]


################################################################################
# Balance positive and negative dataset
################################################################################

# This function accepts two dictionaries, and randomly selects samples from them,
# to maintain a ratio between them.
# If the ratio is set as 0, this step will be skipped.
def dict_balance(dict_a, dict_b, ratio=1):
	num_a = len(dict_a)
	num_b = len(dict_b)
	new_dict_a = dict_a
	new_dict_b = dict_b
	if ratio:
		if num_a > int(num_b*ratio):
			fin_a = int(num_b*ratio)
			new_dict_a = dict_select(dict_a, fin_a)
		else:
			fin_b = int(num_a/ratio)
			new_dict_b = dict_select(dict_b, fin_b)
	return [new_dict_a, new_dict_b]


################################################################################
# Main function
################################################################################

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Preprocessing positive and negative datasets")
	parser.add_argument("-p", "--pos",
						help="Input positive fasta", required=True)
	parser.add_argument("-n", "--neg",
						help="Input negative fasta", required=True)
	parser.add_argument("-o", "--output",
						help="Output prefix", required=True)
	parser.add_argument("-b", "--balance", default=0,
						help="Balacance ratio of positive to negative", required=False)
	parser.add_argument("-r", "--ratio", default=2,
						help="Ratio of training to testing", required=False)
	parser.add_argument("-f", "--fold", default=10,
						help="Folds of cross validation", required=False)
	ARGS = parser.parse_args()

	raw_pos_dict = Util.fasta_to_dict(ARGS.pos, parse=True, sep="")
	pos_dict = dict_filter(raw_pos_dict)[0]
	raw_neg_dict = Util.fasta_to_dict(ARGS.neg, parse=True, sep="")
	neg_dict = dict_filter(raw_neg_dict)[0]
	del_dict = dict_filter(raw_pos_dict)[1] | dict_filter(raw_neg_dict)[1]

	if ARGS.balance:
		pos_dict, neg_dict = dict_balance(pos_dict, neg_dict, ratio=float(ARGS.balance))

	training_ratio = float(ARGS.ratio)/(float(ARGS.ratio)+1)
	pos_training, pos_testing = dict_divide(pos_dict, training_ratio)
	neg_training, neg_testing = dict_divide(neg_dict, training_ratio)

	cv_fold = int(ARGS.fold)
	if cv_fold > 99:
		raise TypeError(
			"CV fold shoulb be no more than 99."
		)	
	# positive cv
	pos_cv_list = []
	pos_cv_dict = pos_training
	for i in range(cv_fold-1):
		div_ratio = 1/(cv_fold-i)
		sampling = dict_divide(pos_cv_dict, div_ratio)
		pos_cv_list.append(sampling[0])
		pos_cv_dict = sampling[1]
	pos_cv_list.append(pos_cv_dict)
	# negative cv  
	neg_cv_list = []
	neg_cv_dict = neg_training
	for i in range(cv_fold-1):
		div_ratio = 1/(cv_fold-i)
		sampling = dict_divide(neg_cv_dict, div_ratio)
		neg_cv_list.append(sampling[0])
		neg_cv_dict = sampling[1]
	neg_cv_list.append(neg_cv_dict)

	# Generating output files
	output_prefix = ARGS.output
	output_filter = output_prefix + ".filtering.fasta"
	Util.dict_to_fasta(del_dict, output_filter)
	output_pos_training = output_prefix + ".pos.training.fasta"
	Util.dict_to_fasta(pos_training, output_pos_training)
	output_pos_testing = output_prefix + ".pos.testing.fasta"
	Util.dict_to_fasta(pos_testing, output_pos_testing)
	output_neg_training = output_prefix + ".neg.training.fasta"
	Util.dict_to_fasta(neg_training, output_neg_training)
	output_neg_testing = output_prefix + ".neg.testing.fasta"
	Util.dict_to_fasta(neg_testing, output_neg_testing)

	with open(output_prefix+".cv_file_list.tsv", "w") as f_out:
		for i in range(cv_fold):
			cv_suffix = str(i).rjust(2,"0")
			cv_pos_training = output_prefix + ".cv" + str(i).rjust(2,"0") + ".pos.training.fasta"
			cv_pos_testing = output_prefix + ".cv" + str(i).rjust(2,"0") + ".pos.testing.fasta"
			cv_neg_training = output_prefix + ".cv" + str(i).rjust(2,"0") + ".neg.training.fasta"
			cv_neg_testing = output_prefix + ".cv" + str(i).rjust(2,"0") + ".neg.testing.fasta"
			new_line = "\t".join([cv_pos_training, cv_pos_testing, cv_neg_training, cv_neg_testing]) + "\n"
			f_out.write(new_line)
			for j in range(cv_fold):
				if j != i:
					Util.dict_to_fasta(pos_cv_list[j], cv_pos_training, mode="a")
					Util.dict_to_fasta(neg_cv_list[j], cv_neg_training, mode="a")
				else:
					Util.dict_to_fasta(pos_cv_list[j], cv_pos_testing)
					Util.dict_to_fasta(neg_cv_list[j], cv_neg_testing)

# _*_ coding: UTF-8 _*_
# Author: PENG Zhao
# Version: 2023-10
""" Instruction
This script provides functions applied in feature ablation.
"""

import argparse
from sklearn import *
import pandas as pd
import pymrmr
import FeatureExtraction
import ScikitModel
from FeatureDataProc import data_dict_read

def load_feature_cord():
	FEATURE_CORD = {}
	with open("Data/FeatureID.list") as f_in:
		for line in f_in:
			elements = line.strip().split("\t")
			FEATURE_CORD[elements[2]] = elements[0:2]
	return FEATURE_CORD



if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Processing feature ablation")
	parser.add_argument("-c", "--cv_file_list",
						help="File recording fasta files for cross validation", required=True)
	parser.add_argument("-o", "--output",
						help="Output prefix", required=True)
	parser.add_argument("-d", "--dict",
						help="Directory containing codonbias, hexamer and ntbias files", required=True)
	parser.add_argument("-f", "--feature_list",
						help="File recording feature subsets in feature ablation", required=True)
	parser.add_argument("-m", "--model",
						help="Model code for sklearn", required=True)
	ARGS = parser.parse_args()


	# Parsing input files
	codon_in = ARGS.dict + "/codonbias.tsv"
	codon_dict = data_dict_read(codon_in)
	hexa_in = ARGS.dict + "/hexamer.tsv"
	hexa_dict = data_dict_read(hexa_in)
	ntbias_in = ARGS.dict + "/ntbias.tsv"
	ntbias_dict = data_dict_read(ntbias_in)

	model = ARGS.model
	FEATURE_CORD = load_feature_cord()

	cv_file_list = []
	with open(ARGS.cv_file_list) as f_in:
		for line in f_in:
			elements = line.strip().split("\t")
			cv_file_list.append(elements)

	ablation_list = []
	with open(ARGS.feature_list) as f_in:
		for line in f_in:
			ablation_list.append(line.strip())


	# Ablation
	for i in range(len(ablation_list)):
		# parsing feature list
		key = "ablation" + str(i)
		feature_list = []
		for j in range(len(ablation_list)):
			if j != i:
				feature_location = FEATURE_CORD[ablation_list[j]]
				for k in range(feature_location[0], feature_location[1]):
					feature_list.append(k)

		# cross validation
		acc_dict[key] = []
		mcc_dict[key] = []
		for input_list in cv_file_list:
			FeatureExtraction.feature_extraction_batch(input_list[0], "tmp_pos_train_file", codon_dict, hexa_dict, ntbias_dict, feature_list)
			FeatureExtraction.feature_extraction_batch(input_list[1], "tmp_pos_test_file", codon_dict, hexa_dict, ntbias_dict, feature_list)
			FeatureExtraction.feature_extraction_batch(input_list[2], "tmp_neg_train_file", codon_dict, hexa_dict, ntbias_dict, feature_list)
			FeatureExtraction.feature_extraction_batch(input_list[3], "tmp_neg_test_file", codon_dict, hexa_dict, ntbias_dict, feature_list)

			X_train, Y_train = ScikitModel.load_feature("tmp_pos_train_file", "tmp_neg_train_file")
			X_test, Y_test = ScikitModel.load_feature("tmp_pos_test_file", "tmp_neg_test_file")
			clf = eval(model)
			acc, mcc = ScikitModel.model_train(X_train, Y_train, X_test, Y_test, clf)
			acc_dict[key].append(str(acc))
			mcc_dict[key].append(str(mcc))

	# generating output files
	with open(output_file+".ablation_acc.tsv", "w") as f_out:
		for key in sorted(acc_dict.keys()):
			new_line = key + "\t" + "\t".join(value) + "\n"
			f_out.write(new_line)


	with open(output_file+".ablation_mcc.tsv", "w") as f_out:
		for key in sorted(mcc_dict.keys()):
			new_line = key + "\t" + "\t".join(value) + "\n"
			f_out.write(new_line)




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
import Util
import ScikitModel
from FeatureDataProc import data_dict_read
from FeatureExtraction import feature_extraction_batch


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Processing feature ablation")
	parser.add_argument("-c", "--cv_file_list",
						help="File recording fasta files for cross validation", required=True)
	parser.add_argument("-o", "--output",
						help="Output prefix", required=True)
	parser.add_argument("-d", "--dict",
						help="Directory containing codonbias, hexamer and ntbias files", required=True)
	parser.add_argument("-f", "--feature_list",
						help="File recording features in feature ablation", required=True)
	parser.add_argument("-m", "--model",
						help="File recording the model", required=True)
	ARGS = parser.parse_args()


	# Parsing input files
	cv_file_list = []
	with open(ARGS.cv_file_list) as f_in:
		for line in f_in:
			elements = line.strip().split("\t")
			cv_file_list.append(elements)

	codon_in = ARGS.dict + "/codonbias.tsv"
	codon_dict = data_dict_read(codon_in)
	hexa_in = ARGS.dict + "/hexamer.tsv"
	hexa_dict = data_dict_read(hexa_in)
	ntbias_in = ARGS.dict + "/ntbias.tsv"
	ntbias_dict = data_dict_read(ntbias_in)

	with open(ARGS.model) as f_in:
		line = f_in.readline()
		model = line.strip("\n")

	ablation_dict = Util.load_feature_list_to_dictionary(ARGS.feature_list)
	ablation_dict["None"] = [] # The initial feature list 


	# cross validation
	acc_dict = {}
	mcc_dict = {}
	for k,v in ablation_dict.items():
		# parsing feature list
		new_key = "ablation:" + k
		feature_list = []
		for j,w in ablation_dict.items():
			if j != k:
				feature_list.extend(w)

		# training model
		acc_dict[new_key] = []
		mcc_dict[new_key] = []
		for input_list in cv_file_list:
			feature_extraction_batch(input_list[0], "tmp_pos_train_file", codon_dict, hexa_dict, ntbias_dict, feature_list)
			feature_extraction_batch(input_list[1], "tmp_pos_test_file", codon_dict, hexa_dict, ntbias_dict, feature_list)
			feature_extraction_batch(input_list[2], "tmp_neg_train_file", codon_dict, hexa_dict, ntbias_dict, feature_list)
			feature_extraction_batch(input_list[3], "tmp_neg_test_file", codon_dict, hexa_dict, ntbias_dict, feature_list)

			X_train, Y_train = ScikitModel.load_feature("tmp_pos_train_file", "tmp_neg_train_file")
			X_test, Y_test = ScikitModel.load_feature("tmp_pos_test_file", "tmp_neg_test_file")
			clf = eval(model)
			acc, mcc = ScikitModel.model_train(X_train, Y_train, X_test, Y_test, clf)
			acc_dict[new_key].append(str(acc))
			mcc_dict[new_key].append(str(mcc))


	# generating output files
	with open(ARGS.output+".ablation_acc.tsv", "w") as f_out:
		for key in sorted(acc_dict.keys()):
			new_line = key + "\t" + "\t".join(acc_dict[key]) + "\n"
			f_out.write(new_line)
	with open(ARGS.output+".ablation_mcc.tsv", "w") as f_out:
		for key in sorted(mcc_dict.keys()):
			new_line = key + "\t" + "\t".join(mcc_dict[key]) + "\n"
			f_out.write(new_line)




# _*_ coding: UTF-8 _*_
# Author: PENG Zhao
# Version: 2023-10
""" Instruction
This script provides functions for model training and saving.
"""

import argparse
from sklearn import *
import os
from joblib import dump, load
import Util
import ScikitModel
from FeatureDataProc import data_dict_read
from FeatureExtraction import feature_extraction_batch

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Processing model training")
	parser.add_argument("-i", "--input_list",
						help="File recording fasta files for model training", required=True)
	parser.add_argument("-o", "--output",
						help="Output directory", required=True)
	parser.add_argument("-d", "--dict",
						help="Directory containing codon_dict, hexa_dict and ntbias_dict files", required=True)
	parser.add_argument("-f", "--feature_list",
						help="File recording feature IDs in IFS", required=True)
	parser.add_argument("-m", "--model",
						help="Model code for sklearn", required=True)
	ARGS = parser.parse_args()

	# Parsing input files
	input_list = []
	with open(ARGS.input_list) as f_in:
		for line in f_in:
			input_list.append(line.strip())

	codon_in = ARGS.dict + "/codonbias.tsv"
	codon_dict = data_dict_read(codon_in)
	hexa_in = ARGS.dict + "/hexamer.tsv"
	hexa_dict = data_dict_read(hexa_in)
	ntbias_in = ARGS.dict + "/ntbias.tsv"
	ntbias_dict = data_dict_read(ntbias_in)

	with open(ARGS.model) as f_in:
		line = f_in.readline()
		model = line.strip("\n")

	feature_list = Util.load_feature_list(ARGS.feature_list)

	score_file = ARGS.output + ".test_score.list"


	# model training
	feature_extraction_batch(input_list[0], "tmp_pos_train_file", codon_dict, hexa_dict, ntbias_dict, feature_list)
	feature_extraction_batch(input_list[1], "tmp_pos_test_file", codon_dict, hexa_dict, ntbias_dict, feature_list)
	feature_extraction_batch(input_list[2], "tmp_neg_train_file", codon_dict, hexa_dict, ntbias_dict, feature_list)
	feature_extraction_batch(input_list[3], "tmp_neg_test_file", codon_dict, hexa_dict, ntbias_dict, feature_list)

	X_train, Y_train = ScikitModel.load_feature("tmp_pos_train_file", "tmp_neg_train_file")
	X_test, Y_test = ScikitModel.load_feature("tmp_pos_test_file", "tmp_neg_test_file")

	clf = eval(model)
	acc, mcc = ScikitModel.model_train(X_train, Y_train, X_test, Y_test, clf)
	dump(clf, "pretrained.model")


	# generating the output directory
	os.system("mkdir %s" % (ARGS.output))
	os.system("mv pretrained.model %s" % (ARGS.output))
	os.system("cp %s %s/codonbias.tsv" % (codon_in, ARGS.output))
	os.system("cp %s %s/hexamer.tsv" % (hexa_in, ARGS.output))
	os.system("cp %s %s/ntbias.tsv" % (ntbias_in, ARGS.output))
	os.system("cp %s %s/feature.list" % (ARGS.feature_list, ARGS.output))
	print("Training complete!\nACC: %f\nMCC: %f" % (acc, mcc))


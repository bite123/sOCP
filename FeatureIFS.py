# _*_ coding: UTF-8 _*_
# Author: PENG Zhao
# Version: 2023-10
""" Instruction
This script provides functions applied in mRMR-IFS for feature selection.
"""

import argparse
from sklearn import *
import pandas as pd
import pymrmr
import Util
import ScikitModel
from FeatureDataProc import data_dict_read
from FeatureExtraction import feature_extraction_batch

# This function accepts A list, 
# and returns B list containing incremental subsets containing the first n elements of A.
def feature_ifs(ranking_list):
	subset_list = []
	for i in range(len(ranking_list)):
		subset = ranking_list[:i+1]
		subset_list.append(subset)
	return subset_list

# This function accepts two lists:
# List A consists of several inner-layer list.
# List B consists of several numbers.
# It returns the filtered A, with inner-layer list selected according to List B as coordinates.
# To apply in this script, Number in List B = List coordinate + 1
def double_layer_list_selection(ori_list, number_list):
	new_list = []
	for l in ori_list:
		new_l = []
		for i in range(len(l)):
			if i+1 in number_list:
				new_l.append(l[i])
		new_list.append(new_l)
	return new_list

# This function accepts a double-layer list with isometric inner-layer list,
# and transform it to a new double-layer list by a matrix transposition.
def double_layer_list_transposition(ori_list):
	new_list = []
	for i in range(len(ori_list[0])):
		new_l = []
		for l in ori_list:
			new_l.append(l[i])
		new_list.append(new_l)
	return new_list


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Processing feature selection through mRMR-IFS")
	parser.add_argument("-c", "--cv_file_list",
						help="File recording fasta files for cross validation", required=True)
	parser.add_argument("-o", "--output",
						help="Output prefix", required=True)
	parser.add_argument("-d", "--dict",
						help="Directory containing codon_dict, hexa_dict and ntbias_dict files", required=True)
	parser.add_argument("-f", "--feature_list",
						help="File recording feature IDs in IFS", required=True)
	parser.add_argument("-m", "--model",
						help="Model code for sklearn", required=True)
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

	feature_list = Util.load_feature_list(ARGS.feature_list)
	feature_map = {} # record the feature ID
	for i in range(len(feature_list)):
		order = i + 1
		feature_map[str(order)] = str(feature_list[i])

	# cross validation
	mr_total = []
	acc_total = []
	mcc_total = []
	for input_list in cv_file_list:
		# loading files
		feature_extraction_batch(input_list[0], "tmp_pos_train_file", codon_dict, hexa_dict, ntbias_dict, feature_list)
		feature_extraction_batch(input_list[1], "tmp_pos_test_file", codon_dict, hexa_dict, ntbias_dict, feature_list)
		feature_extraction_batch(input_list[2], "tmp_neg_train_file", codon_dict, hexa_dict, ntbias_dict, feature_list)
		feature_extraction_batch(input_list[3], "tmp_neg_test_file", codon_dict, hexa_dict, ntbias_dict, feature_list)

		X_train, Y_train = ScikitModel.load_feature("tmp_pos_train_file", "tmp_neg_train_file")
		X_test, Y_test = ScikitModel.load_feature("tmp_pos_test_file", "tmp_neg_test_file")


		# Step 1, mRMR
		pos_training_in = pd.read_csv("tmp_pos_train_file", header=None, sep="\t")
		pos_training_in[[0]] = 1
		neg_training_in = pd.read_csv("tmp_neg_train_file", header=None, sep="\t")
		neg_training_in[[0]] = 0
		df = pd.concat([pos_training_in, neg_training_in])

		columns_id = []
		for i in df.columns:
			columns_id.append(str(i)) # column id should be string
		df.columns = columns_id
		feature_number = df.shape[1] - 1

		mr = pymrmr.mRMR(df, "MIQ", feature_number)


		# Step 2, IFS
		ranking_list = [int(i) for i in mr]
		subset_list = feature_ifs(ranking_list)
		acc_eval = []
		mcc_eval = []

		for subset in subset_list:
			new_X_train = double_layer_list_selection(X_train, subset)
			new_X_test = double_layer_list_selection(X_test, subset)
			clf = eval(model)
			acc, mcc = ScikitModel.model_train(new_X_train, Y_train, new_X_test, Y_test, clf)
			acc_eval.append(str(acc))
			mcc_eval.append(str(mcc))

		mr_total.append(mr)
		acc_total.append(acc_eval)
		mcc_total.append(mcc_eval)


	# generating output files
	mr_total = double_layer_list_transposition(mr_total)
	acc_total = double_layer_list_transposition(acc_total)
	mcc_total = double_layer_list_transposition(mcc_total)

	with open(ARGS.output+".mRMR_rank.tsv", "w") as f_out:
		for i in range(len(mr_total)):
			new_list = [feature_map[k] for k in mr_total[i]]
			new_line = "\t".join(new_list) + "\n"
			f_out.write(new_line)	
	with open(ARGS.output+".mRMR_acc.tsv", "w") as f_out:
		for i in range(len(acc_total)):
			new_line = "\t".join(acc_total[i]) + "\n"
			f_out.write(new_line)
	with open(ARGS.output+".mRMR_mcc.tsv", "w") as f_out:
		for i in range(len(mcc_total)):
			new_line = "\t".join(mcc_total[i]) + "\n"
			f_out.write(new_line)
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
import FeatureExtraction
import ScikitModel
from FeatureDataProc import data_dict_read
from FeatureAblation import load_feature_cord

def feature_ifs(ranking_list):
	subset_list = []
	for i in range(len(ranking_list)):
		subset = ranking_list[:i+1]
		subset_list.append(subset)
	return subset_list

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

	feature_list = []
	with open(ARGS.feature_list) as f_in:
		for line in f_in:
			feature_id = line.strip()
			feature_location = FEATURE_CORD[feature_id]
			for i in range(feature_location[0], feature_location[1]):
				feature_list.append(i)

	# cross validation
	mr_total = []
	acc_total = []
	mcc_total = []
	for input_list in cv_file_list:
		FeatureExtraction.feature_extraction_batch(input_list[0], "tmp_pos_train_file", codon_dict, hexa_dict, ntbias_dict, feature_list)
		FeatureExtraction.feature_extraction_batch(input_list[1], "tmp_pos_test_file", codon_dict, hexa_dict, ntbias_dict, feature_list)
		FeatureExtraction.feature_extraction_batch(input_list[2], "tmp_neg_train_file", codon_dict, hexa_dict, ntbias_dict, feature_list)
		FeatureExtraction.feature_extraction_batch(input_list[3], "tmp_neg_test_file", codon_dict, hexa_dict, ntbias_dict, feature_list)

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
			X_train = 从原X_train中取subset
			X_test = 
			clf = eval(model)
			clf = clf.fit(X_train, Y_train)
			pred = clf.predict(X_test)
			acc, mcc = model_train(X_train, Y_train, X_test, Y_test, clf)
			acc_eval.append(str(acc))
			mcc_eval.append(str(mcc))

		mr_total.append(mr)
		acc_total.append(acc_eval)
		mcc_total.append(mcc_eval)

	# generating output files
	with open(output_file+".mRMR_mr.tsv", "w") as f_out:
		for i in range(len(mr_total)):
			new_line = "\t".join(mr_total[i]) + "\n"
			f_out.write(new_line)	

	with open(output_file+".mRMR_acc.tsv", "w") as f_out:
		for i in range(len(acc_total)):
			new_line = "\t".join(acc_total[i]) + "\n"
			f_out.write(new_line)

	with open(output_file+".mRMR_mcc.tsv", "w") as f_out:
		for i in range(len(mcc_total)):
			new_line = "\t".join(mcc_total[i]) + "\n"
			f_out.write(new_line)
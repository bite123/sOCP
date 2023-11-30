# _*_ coding: UTF-8 _*_
# Author: PENG Zhao
# Version: 2023-10
""" Instruction
This script provides functions related to scikit-learn package,
and a one-step way to compare common algorithms.
"""

import argparse
import pandas as pd
import sklearn
from sklearn import *
from sklearn.linear_model import LogisticRegression
from sklearn import svm
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import accuracy_score, matthews_corrcoef
import Util
from FeatureDataProc import data_dict_read
from FeatureExtraction import feature_extraction_batch

# This function accepts feature files from positive and negative dataset,
# and returns two data lists: feature and label.
def load_feature(input_pos, input_neg):
	pos_feature = pd.read_csv(input_pos, header=None, sep="\t")
	neg_feature = pd.read_csv(input_neg, header=None, sep="\t")
	dim = pos_feature.columns.size
	X_pos = pos_feature.values[:,1:dim].astype(float).tolist()
	X_neg = neg_feature.values[:,1:dim].astype(float).tolist()
	X_total = X_pos + X_neg
	Y_total = [1.0] * len(X_pos) + [0.0] * len(X_neg)
	return [X_total, Y_total]

# This function accepts X data list (features) and Y data list (labels), from training and testing dataset.
# Then fitting a scikit-learn model (clf) to the data, and evaluating the model with accuracy and MCC.
def model_train(X_train, Y_train, X_test, Y_test, clf):
	clf = clf.fit(X_train, Y_train)
	pred = clf.predict(X_test)
	acc = accuracy_score(pred, Y_test)
	mcc = matthews_corrcoef(pred, Y_test)
	return [acc, mcc]


# Main function
''' The format of Cross-validation feature file list is a tsv file as below:
cv/cv01_pos_train    cv/cv01_pos_test    cv/cv01_neg_train    cv/cv01_neg_test
cv/cv02_pos_train    cv/cv02_pos_test    cv/cv02_neg_train    cv/cv02_neg_test
... ...
cv/cv10_pos_train    cv/cv10_pos_test    cv/cv10_neg_train    cv/cv10_neg_test
'''
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Comparing scikit-learn models")
	parser.add_argument("-c", "--cv_file_list",
						help="File recording fasta files for cross validation", required=True)    
	parser.add_argument("-o", "--output",
						help="Output prefix", required=True)
	parser.add_argument("-d", "--dict",
						help="Directory containing codonbias, hexamer and ntbias files", required=True)
	parser.add_argument("-m", "--model_list",
						help="File recording candidata scikit-learn models", required=True)
	parser.add_argument("-f", "--feature_list",
						help="File recording feature list used in model training", required=True)
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

	model_list = []
	with open(ARGS.model_list) as f_in:
		for line in f_in:
			model_list.append(line.strip("\n"))

	feature_list = Util.load_feature_list(ARGS.feature_list)

	# cross validation
	acc_dict = {}
	mcc_dict = {}
	for model in model_list:
		acc_dict[model] = []
		mcc_dict[model] = []
		for input_list in cv_file_list:
			feature_extraction_batch(input_list[0], "tmp_pos_train_file", codon_dict, hexa_dict, ntbias_dict, feature_list)
			feature_extraction_batch(input_list[1], "tmp_pos_test_file", codon_dict, hexa_dict, ntbias_dict, feature_list)
			feature_extraction_batch(input_list[2], "tmp_neg_train_file", codon_dict, hexa_dict, ntbias_dict, feature_list)
			feature_extraction_batch(input_list[3], "tmp_neg_test_file", codon_dict, hexa_dict, ntbias_dict, feature_list)		

			X_train, Y_train = load_feature("tmp_pos_train_file", "tmp_neg_train_file")
			X_test, Y_test = load_feature("tmp_pos_test_file", "tmp_neg_test_file")
			clf = eval(model)
			acc, mcc = model_train(X_train, Y_train, X_test, Y_test, clf)
			acc_dict[model].append(str(acc))
			mcc_dict[model].append(str(mcc))

	# generating output files
	with open(ARGS.output+".model_acc.tsv", "w") as f_out:
		for key,value in acc_dict.items():
			new_line = key + "\t" + "\t".join(value) + "\n"
			f_out.write(new_line)
	with open(ARGS.output+".model_mcc.tsv", "w") as f_out:
		for key,value in mcc_dict.items():
			new_line = key + "\t" + "\t".join(value) + "\n"
			f_out.write(new_line)
	
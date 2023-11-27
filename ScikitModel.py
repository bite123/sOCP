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
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Comparing scikit-learn models")
	parser.add_argument("-c", "--cross_validation",
						help="Cross-validation feature file list", required=True)    
	parser.add_argument("-o", "--output",
						help="Output prefix", required=True)
	parser.add_argument("-d", "--dict",
						help="Directory containing codonbias, hexamer and ntbias files", required=True)
	parser.add_argument("-m", "--model_list",
						help="A file recording candidata scikit-learn models", required=True)
	parser.add_argument("-f", "--feature_list",
						help="A file recording feature list used in model training", required=True)
	ARGS = parser.parse_args()

	''' The format of Cross-validation feature file list is a tsv file as below:
	cv/cv01_pos_train    cv/cv01_pos_test    cv/cv01_neg_train    cv/cv01_neg_test
	cv/cv02_pos_train    cv/cv02_pos_test    cv/cv02_neg_train    cv/cv02_neg_test
	... ...
	cv/cv10_pos_train    cv/cv10_pos_test    cv/cv10_neg_train    cv/cv10_neg_test
	'''
	cv_feature_list = []
	with open(ARGS.cross_validation) as f_in:
		for line in f_in:
			elements = line.strip().split("\t")
			cv_feature_list.append(elements)

	model_list = []
	with open(ARGS.model_list) as f_in:
		for line in f_in:
			model_list.append(line.strip("\n"))

	acc_dict = {}
	mcc_dict = {}
	for model in model_list:
		acc_dict[model] = []
		mcc_dict[model] = []
		for i in range(len(cv_feature_list)):
			X_train, Y_train = load_feature(cv_feature_list[i][0], cv_feature_list[i][2])
			X_test, Y_test = load_feature(cv_feature_list[i][1], cv_feature_list[i][3])
			clf = eval(model)
			acc, mcc = model_train(X_train, Y_train, X_test, Y_test, clf)
			acc_dict[model].append(str(acc))
			mcc_dict[model].append(str(mcc))

	with open(ARGS.output+".model_acc.tsv", "w") as f_out:
		for key,value in acc_dict.items():
			new_line = key + "\t" + "\t".join(value) + "\n"
			f_out.write(new_line)
	with open(ARGS.output+".model_mcc.tsv", "w") as f_out:
		for key,value in mcc_dict.items():
			new_line = key + "\t" + "\t".join(value) + "\n"
			f_out.write(new_line)
	
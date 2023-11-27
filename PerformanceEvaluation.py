# _*_ coding: UTF-8 _*_
# Author: PENG Zhao
# Version: 2023-10
""" Instruction
This script provides functions evaluating model performance.
"""

import argparse

# This function is used in case dividing by zero.
def div(a, b):
	if not b:
		c = "NaN"
	else:
		c = round(a/b, 8)
	return c

################################################################################
# Preprocessing the prediction file and the actual label file
################################################################################

# The prediction file is generated in PredictingCodingPotential.py,
# and the actual label file is generated in dataset processing.
def combine_pred_act(output_file, prediction_file, actual_file):
	prediction_dict = Util.tsv_to_dict(prediction_file)
	actual_dict = Util.tsv_to_dict(actual_file)
	with open(output_file, "w") as f_out:
		for key,value in prediction_dict.items():
			if key in actual_dict:
				act_value = actual_dict[key]
			new_line = "\t".join([key, value, act_value])
			f_out.write(new_line)

################################################################################
# Evaluating function
################################################################################

# This function accepts a preprocessed file from combine_pred_act(),
# and returns a list including: 
# [TPR, TNR, FPR, FNR, PRE, ACC, F1S, HM, MCC]
# The preprocessed file contains:
# Column 1, the ORF ID
# Column 2, the predict value, a score ranged from 0 to 1
# Column 3, the actual value, 0/1
def eval_metric(input_file, cutoff=0.5):
	with open(input_file) as f_in:
		tp, tn, fp, fn = 0, 0, 0, 0
		for line in f_in:
			elements = line.strip().split("\t")
			pred = float(elements[1])
			act = int(float(elements[2]))
			if act:
				if pred < cutoff:
					fn += 1
				else:
					tp += 1
			else:
				if pred < cutoff:
					tn += 1
				else:
					fp += 1
	p = tp + fn
	n = tn + fp
	tpr = div(tp, p)
	tnr = div(tn, n)
	fpr = div(fp, n)
	fnr = div(fn, p)
	pre = div(tp, tp + fp)
	acc = div(tp + tn, p + n)
	if pre != "NaN" and tpr != "NaN":
		f1s = div(2 * pre * tpr, pre + tpr)
	else:
		f1s = "NaN"
	if tpr != "NaN" and tnr != "NaN":
		hm = div(2 * tpr * tnr, tpr + tnr)
	else:
		hm = "NaN"
	mcc_deno = ((tp + fn) * (tp + fp) * (tn + fp) * (tn + fn))**(1/2)
	mcc = div(tp * tn - fp * fn, mcc_deno)
	return [tpr, tnr, fpr, fnr, pre, acc, f1s, hm, mcc]

# This function accepts a preprocessed file from combine_pred_act(),
# returns a list used in roc_curve() and pr_curve().
def parse_for_roc(input_file):
	with open(input_file) as f_in:
		score_dict = {}
		label_dict = {}
		i, j = 0, 0
		for line in f_in:
			elements = line.strip().split("\t")
			pred = float(elements[1])
			act = int(float(elements[2]))
			score_dict[i] = pred
			label_dict[i] = act
			i += 1
			j += act
		p = j
		n = i - j
	return [score_dict, label_dict, p, n]

################################################################################
# Generating a file recording metrics for cutoffs
################################################################################

# This function accepts a preprocessed file from combine_pred_act(),
# returns evaluation metrics for cutoffs,
def cutoff_metric(input_file, output_file):
	# getting total predicted score list
	total_score = []
	with open(input_file) as f_in:
		for line in f_in:
			elements = line.strip().split("\t")
			pred = float(elements[1])
			total_score.append(pred)
	total_score.extend([0.0, 1.0])
	total_score = list(set(total_score))
	# output
	with open(output_file, "w") as f_out:
		for score in total_score:
			metric_list = eval_metric(input_file, cutoff = score)
			metric_list = [str(i) for i in metric_list]
			new_line = str(score) + "\t" + "\t".join(metric_list) + "\n"
			f_out.write(new_line)
	return eval_metric(input_file, cutoff=0.5)

################################################################################
# AUC ROC/PR and curve plotting
################################################################################

# This function accepts a preprocessed file from combine_pred_act(),
# returns AUC ROC and generates a file for ROC plotting. 
def roc_curve(input_file, output_file):
	score_dict, label_dict, p, n = parse_for_roc(input_file)
	tp, fp, tpr, fpr, auc = 0, 0, 0, 0, 0
	plot_list = []
	score_tuplelist = sorted(score_dict.items(), key=lambda x: x[1], reverse=True)
	for score_tuple in score_tuplelist:
		score = score_tuple[1]
		label = label_dict[score_tuple[0]]
		if label:
			tp += 1
			tpr = div(tp, p)
			adding_area = 0
		else:
			fp += 1
			new_fpr = div(fp, n)
			adding_area = round((new_fpr - fpr) * tpr, 8)
			fpr = new_fpr
		plot_list.append([fpr, tpr])
		auc += adding_area
	with open(output_file, "w") as f_out:
		for plot in plot_list:
			line = str(plot[0]) + "\t" + str(plot[1]) + "\n"
			f_out.write(line)
	return auc

# This function accepts a preprocessed file from combine_pred_act(),
# returns AUC PR and generates a file for PR plotting. 
def pr_curve(input_file, output_file):
	score_dict, label_dict, p, n = parse_for_roc(input_file)
	tp, fp, tpr, pre, auc = 0, 0, 0, 1, 0
	plot_list = []
	score_tuplelist = sorted(score_dict.items(), key=lambda x: x[1], reverse=True)
	for score_tuple in score_tuplelist:
		score = score_tuple[1]
		label = label_dict[score_tuple[0]]
		if label:
			tp += 1
		else:
			fp += 1
		new_tpr = div(tp, p)
		new_pre = div(tp, (tp + fp))
		adding_area = round((new_tpr - tpr) * (new_pre + pre)/2, 8)
		tpr, pre = new_tpr, new_pre
		plot_list.append([tpr, pre])
		auc += adding_area
	with open(output_file, "w") as f_out:
		for plot in plot_list:
			line = str(plot[0]) + "\t" + str(plot[1]) + "\n"
			f_out.write(line)
	return auc


################################################################################
# Main function
################################################################################
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Evaluating preformance of model")
	parser.add_argument("-i", "--input",
						help="Input file recording prediction scores", required=True)
	parser.add_argument("-a", "--actual",
						help="File recording actual labels", required=True)
	parser.add_argument("-o", "--output",
						help="Output prefix", required=True)
	ARGS = parser.parse_args()

	prediction_file = ARGS.input
	act_file = ARGS.actual
	output_prefix = ARGS.output
	output_preproc = output_prefix + ".preprocessed.tsv"
	output_metric = output_prefix + ".metric.tsv"
	output_roc = output_prefix + ".roc.tsv"
	output_pr = output_prefix + ".pr.tsv"

	combine_pred_act(output_preproc, prediction_file, actual_file)

	metric_id = [TPR, TNR, FPR, FNR, PRE, ACC, F1S, HM, MCC]
	metric_list = cutoff_metric(output_preproc, output_metric)
	print("Metrics (cutoff = 0.5):")
	for i in range(len(metric_id)):
		print(metric_id[i] + "\t" + str(metric_list[i]))

	auc_roc = roc_curve(output_preproc, output_roc)
	print("AUC ROC:")
	print(auc_roc)

	auc_pr = pr_curve(output_preproc, output_pr)
	print("AUC PR:")
	print(auc_pr)


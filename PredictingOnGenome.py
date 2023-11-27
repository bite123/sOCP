# _*_ coding: UTF-8 _*_
# Author: PENG Zhao
# Version: 2023-10
""" Instruction
This script accepts a genome fasta,
and predicts potential coding smORFs using a pre-trained model.
""" 

import argparse
import string
from joblib import dump, load
from FeatureDataProc import data_dict_read
from FeatureExtraction import Ntseq, feature_extraction
import Util

################################################################################
# Processing function
################################################################################
def load_cutoff(input_file):
	# This function reads the cutoff configuration file.
	start_codon_cutoff = {}
	with open(input_file) as f_in:
		for line in f_in:
			elements = line.strip().split("\t")
			start_codon = elements[0]
			cutoff = float(elements[1])
			if cutoff:
				start_codon_cutoff[start_codon] = cutoff
	return start_codon_cutoff

# This function accepts a nucleotide sequence, and returns its reverse complement.
# Upper case will be transferred.
# The undefined nucleotide (N, R, Y, etc.) will all be transfered to N.
def reverse_complement(seq):
	seq = seq[::-1].upper()
	map_dict = {}
	for i in string.ascii_uppercase:
		map_dict[i] = "N"
	map_dict.update({"A":"T", "T":"A", "G":"C", "C":"G"})
	new_seq = "".join([map_dict[i] for i in seq])
	return new_seq

################################################################################
# Main function
################################################################################
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Predicting smORFs from genome")
	parser.add_argument("-i", "--input",
						help="Input fasta file", required=True)
	parser.add_argument("-o", "--output",
						help="Output prefix", required=True)
	parser.add_argument("-m", "--model",
						help="Directory containing model files", required=True)
	parser.add_argument("-c", "--cutoff",
						help="Cutoff configuration file", required=True)
	parser.add_argument("-l", "--min_len",
						help="Minimum number of codons (excluding stop codon)", required=False)
	parser.add_argument("-L", "--max_len",
						help="Maximum number of codons (excluding stop codon)", required=False)    
	ARGS = parser.parse_args()

	################################################################################
	# Processing input, output, and cutoff file
	################################################################################
	seq_dict = Util.fasta_to_dict(ARGS.input, parse=True, sep="")
	start_codon_cutoff = load_cutoff(ARGS.cutoff)

	output_fasta_file = ARGS.output + ".fasta"
	output_upsite_file = ARGS.output + ".upsite.list"
	output_bed_file = ARGS.output + ".bed"

	################################################################################
	# Processing length setting
	################################################################################
	if not ARGS.min_len:
		min_len = 10
	else:
		min_len = int(ARGS.min_len) - 1
	if not ARGS.max_len:
		max_len = 100
	else:
		max_len = int(ARGS.max_len)

	################################################################################
	# Processing model files
	################################################################################
	codon_in = ARGS.model + "/codonbias.tsv"
	codon_dict = data_dict_read(codon_in)
	hexa_in = ARGS.model + "/hexamer.tsv"
	hexa_dict = data_dict_read(hexa_in)
	ntbias_in = ARGS.model + "/ntbias.tsv"
	ntbias_dict = data_dict_read(ntbias_in)

	f_list_in = ARGS.model + "/feature_id.list"
	feature_list = Util.load_feature_list(f_list_in)

	model_in = ARGS.model + "/pretrained.model"
	clf = load(model_in)

	################################################################################
	# Predicting and generating output files
	################################################################################
	with open(output_fasta_file, "w") as fa_out, open(output_upsite_file, "w") as ef_out, open(output_bed_file, "w") as bed_out:
		for key,value in seq_dict.items():
			seq_id = key
			seq = value.upper()
			# processing the sense strand
			for i in range(len(seq)):
				if seq[i:i+3] in start_codon_cutoff:
					start_position = i
					start_codon = seq[i:i+3]
					for j in range(i+3, i+3*max_len+3, 3):
						if set(seq[j:j+3])|{'A', 'C', 'G', "T"} != {'A', 'C', 'G', "T"}:
							break
						if seq[j:j+3] in ["TAA","TGA","TAG"]:
							len_orf = j - i + 3
							if len_orf < 3*min_len + 3:
								break
							end_position = j + 3
							orf = seq[i:j+3]
							if i < 3:
								extended_orf = (3-i)*"N" + seq[0:j+3]
							else:
								extended_orf = seq[i-3:j+3]
							# predicting coding potential
							orf_feature = feature_extraction(extended_orf, codon_dict, hexa_dict, ntbias_dict, feature_list)
							score = clf.predict_proba([orf_feature])[0][1]
							if score >= start_codon_cutoff[start_codon]:
								score = str(round(1000*score))
								orf_id = seq_id + "#" + str(start_position) + "#" + "+" + "#" + start_codon + "#" + score
								bed = [seq_id, str(start_position), str(end_position), orf_id, "+", score]
								fa_out.write(">" + orf_id + "\n" + orf + "\n")
								ef_out.write(orf_id + "\t" + extended_orf[0:3] + "\n")						
								bed_out.write("\t".join(bed_l) + "\n")
							break
			# processing the anti-sense strand
			rc_seq = reverse_complement(seq)
			for i in range(len(rc_seq)): 
				if rc_seq[i:i+3] in start_codon_cutoff:
					end_position = len(rc_seq) - i
					start_codon = rc_seq[i:i+3]
					for j in range(i+3, i+3*max_len+3, 3):
						if set(rc_seq[j:j+3])|{'A', 'C', 'G', "T"} != {'A', 'C', 'G', "T"}:
							break
						if rc_seq[j:j+3] in ["TAA","TGA","TAG"]:
							len_orf = j - i + 3
							if len_orf < 3*min_len + 3:
								break
							start_position = len(rc_seq) - (j + 3)
							orf = rc_seq[i:j+3]
							if i < 3:
								extended_orf = (3-i)*"N" + rc_seq[0:j+3]
							else:
								extended_orf = rc_seq[i-3:j+3]
							# predicting coding potential
							orf_feature = feature_extraction(extended_orf, codon_dict, hexa_dict, ntbias_dict, feature_list)
							score = clf.predict_proba([orf_feature])[0][1]
							if score >= start_codon_cutoff[start_codon]:
								score = str(round(1000*score))
								orf_id = seq_id + "#" + str(end_position) + "#" + "-" + "#" + start_codon + "#" + score
								bed = [seq_id, str(start_position), str(end_position), orf_id, "-", score]
								fa_out.write(">" + orf_id + "\n" + orf + "\n")
								ef_out.write(orf_id + "\t" + extended_orf[0:3] + "\n")
								bed_out.write("\t".join(bed) + "\n")
							break


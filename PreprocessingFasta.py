# _*_ coding: UTF-8 _*_
# Author: PENG Zhao
# Version: 2024-10
""" Instruction
This script provides functions preprocessing the input fasta corresonding file,
to acquire extended ORFs or the upstreaming 3-nt for sOCP models.
"""

import argparse
import re
import Util
from Bio import SeqIO
from Bio.Seq import Seq

################################################################################
# ncRNA
################################################################################

# This function accepts a SHORT nucleotide sequence, 
# and returns the extended ORF with the maximum length within the sequence.
# UPPER CASE IS NEEDED!
# SHORT indicates the sequence is commonly an RNA, not a chromosome.
# For extended ORF, if there is an incomplete upstreaming 3-nt, using N instead.
# If no ORF >= 11 codons is found, return "".
# Also returning length, start position and end position.
def extract_orf_from_short_seq(seq):
	ful_len = len(seq)
	orf_dict = {}
	for i in range(ful_len):
		if seq[i:i+3] == "ATG":
			start_position = i + 1
			for j in range(i+3, ful_len-2, 3):
				if seq[j:j+3] in ["TAA", "TGA", "TAG"]:
					end_position = j + 3
					len_orf = j - i + 3
					orf_dict[start_position] = [end_position, len_orf]
					break
	f_start = ful_len
	max_len = 32
	for key,value in orf_dict.items():
		if 303 >= value[1] > max_len:
			max_len = value[1]
			f_start = key
		elif value[1] == max_len:
			if key < f_start:
				f_start = key
	if f_start == ful_len:
		return ["", 0, 0, 0]
	else:
		f_end = orf_dict[f_start][0]
		seq = "NNN" + seq
		max_seq = seq[f_start-1:f_end+3]
		max_len = orf_dict[f_start][1]
		return [max_seq, f_start, f_end, max_len]

# This function extracts the longest ORF from an input ncRNA fasta.
def preproc_ncRNA(input_file, output_file):
	seq_dict = Util.fasta_to_dict(input_file, parse=True, sep="")
	with open(output_file, "w") as f_out:
		for k,v in seq_dict.items():
			orf = extract_orf_from_short_seq(v)[0]
			if orf:
				f_out.write(k)
				f_out.write(orf + "\n")


################################################################################
# CDS and mRNA
################################################################################

# This function extracts the upsteaming 3-nt, based on a CDS file and its corresponding mRNA file.
# The two fasta files can be retrieved from NCBI RefSeq, and will have the following formats:
# 1) The CDS and mRNA are in one-to-one correspondence in order.
# 2) There is location information in the header line, e.g.:
# >lcl|NM_001419531.1_prot_NP_001406460.1_1 ... [location=256..4701] ...

def preproc_CDS_mRNA(input_cds, input_mrna, output_file):
	# parse CDS to get location
	pat = re.compile(r'location=\D*(\d+)')
	with open(input_cds) as f_in:
		pos_dict = {}
		headline_dict = {}
		sequence_dict = {}
		i = 0
		for line in f_in:
			if line.startswith(">"):
				i += 1
				headline_dict[i] = line
				sequence_dict[i] = []
				start_position = int(re.search(pat, line).group(1))
				pos_dict[i] = start_position
			else:
				sequence_dict[i].append(line)

	# parse mRNA to get the upstreaming 3-nt
	with open(input_mrna) as f_in:
		mrna_dict = {}
		up_dict = {}
		i = 0
		for line in f_in:
			if line.startswith(">"):
				i += 1
				mrna_dict[i] = []
			else:
				mrna_dict[i].append(line.strip())
		for j in range(len(mrna_dict)):
			j += 1
			mrna = "".join(mrna_dict[j])
			start_position = pos_dict[j]
			if start_position <= 3:
				up_site = "N" * (3-start_position+1) + mrna[0:start_position-1]
			else:
				up_site = mrna[start_position-4:start_position-1]
			up_dict[j] = up_site

	# generate output files
	with open(input_cds) as f_in, open(output_file,"w") as f_out:
		for i in range(len(up_dict)):
			i += 1
			headline = headline_dict[i]
			f_out.write(headline)
			seqline = up_dict[i]
			for seq in sequence_dict[i]:
				seqline += seq.strip()
			seqline += "\n"
			f_out.write(seqline)


################################################################################
# CDS, genome and bed
################################################################################

# This function accepts <class 'Bio.Seq.Seq'>, strand, start site, and end site,
# returns a list including the start codon and the upstreaming 3-nt 
def seq_upsite(seq, strand, start, end):
	if strand == "+":
		start_codon = str(seq[start:start+3]).upper()
		up_codon = str(seq[start-3:start]).upper()
	elif strand == "-":
		start_codon = str(seq[end-3:end].reverse_complement()).upper()
		up_codon = str(seq[end:end+3].reverse_complement()).upper()
	return [start_codon, up_codon]

# This function extracts the upsteaming 3-nt, based on a CDS file and its corresponding bed.
# The genome fasta file is required as well.
# If the chromosome names in the first column of bed file do not match
# the header line of the genome chromosomes, a tsv file as follow is required:
# Column 1: Chromosome ID in genome fasta; Column 2: Chromosome ID in bed

def preproc_CDS_bed(input_cds, input_genome, input_bed, output_file, chr_id_map_file=""):
	# parse chr_id_map_file
	if chr_id_map_file:
		chr_map = {}
		with open(chr_id_map_file) as f_in:
			for line in f_in:
				elements = line.strip().split("\t")
				chr_map[elements[0]] = elements[1]

	# read fasta
	seq_dict = {}
	chr_list = list(SeqIO.parse(input_genome, "fasta"))
	for i in range(len(chr_list)):
		if chr_map and chr_list[i].id in chr_map:
			seq_id = chr_map[chr_list[i].id]
		else:
			seq_id = chr_list[i].id
		seq_dict[seq_id] = chr_list[i].seq

	# parse bed file
	up_dict = {}
	start_dict = {}
	with open(input_bed) as f_in:
		for line in f_in:
			elements = line.strip().split("\t")
			s_id = elements[3]
			s_chr = elements[0]
			s_strand = elements[5]
			s_start = int(elements[1])
			s_end = int(elements[2])
			seq = seq_dict[s_chr]
			start_codon, up_codon = seq_upsite(seq, s_strand, s_start, s_end)
			up_dict[s_id] = up_codon
			start_dict[s_id] = start_codon

	# generate output file
	with open(input_cds) as f_in:
		headline_dict = {}
		sequence_dict = {}
		i = 0
		for line in f_in:
			if line.startswith(">"):
				i += 1
				headline_dict[i] = line
				sequence_dict[i] = []
			else:
				sequence_dict[i].append(line)

	with open(output_file, "w") as f_out:
		for i in range(len(headline_dict)):
			i += 1
			f_out.write(headline_dict[i])
			cds_id = headline_dict[i].lstrip(">").strip()
			seqline = ""
			for seq in sequence_dict[i]:
				seqline += seq.strip()
			if cds_id in start_dict and seqline[0:3] == start_dict[cds_id]:
				seqline = up_dict[cds_id] + seqline + "\n"
			else:
				seqline = "NNN" + seqline + "\n"
			f_out.write(seqline)


################################################################################
# Filtering extended ORFs
################################################################################
# This function accepts extended ORFs, and removed those:
# 1) not multiples of 3
# 2) have non-ACTG nucleotides within
# 3) have early stop codons 
# 4) have an out-of-range length (not within 11-101 codons, including the stop codon)
def check_extended_orf(extended_orf, min_l=10, max_l=100):
	min_length = 3*min_l + 6
	max_length = 3*max_l + 6
	length = len(extended_orf)
	if length%3 == 0 and length <= max_length and length >= min_length:
		if_removed = 0
		coding_region = extended_orf[3:length-3]
		for i in range(len(coding_region)//3):
			codon = coding_region[3*i:3*i+3].upper()
			if codon not in Util.CODON_TABLE or Util.CODON_TABLE[codon] == "*":
				if_removed = 1
				break
	else:
		if_removed = 1
	return if_removed

def filtering_fasta(input_file, output_file):
	with open(input_file) as f_in:
		headline_dict = {}
		sequence_dict = {}
		i = 0
		for line in f_in:
			if line.startswith(">"):
				i += 1
				headline_dict[i] = line
				sequence_dict[i] = []
			else:
				sequence_dict[i].append(line)

	with open(output_file, "w") as f_out:
		for i in range(len(headline_dict)):
			i += 1
			seqline = ""
			for seq in sequence_dict[i]:
				seqline += seq.strip()
			if check_extended_orf(seqline) == 0:
				f_out.write(headline_dict[i])
				f_out.write(seqline + "\n")


################################################################################
# Main function
################################################################################
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Preprocessing to get extended ORFs")
	parser.add_argument("-t", "--type",
						help="Type of preprocessing: ncRNA, mRNA, genome", required=True)
	parser.add_argument("-i", "--input",
						help="Input ORF fasta", required=True)
	parser.add_argument("-o", "--output",
						help="Output extended ORF fasta", required=True)
	parser.add_argument("-m", "--mrna",
						help="Parameter for type mRNA: input mRNA fasta", required=False)
	parser.add_argument("-g", "--genome",
						help="Parameter for type genome: input genome fasta", required=False)
	parser.add_argument("-b", "--bed",
						help="Parameter for type genome: input bed", required=False)
	parser.add_argument("-l", "--list", default="",
						help="Parameter for type genome: input chromosome ID mapping list", required=False)
	ARGS = parser.parse_args()

	input_file = ARGS.input
	output_file = ARGS.output
	proc_type = ARGS.type
	if proc_type == "ncRNA":
		preproc_ncRNA(input_file, output_file)
	elif proc_type == "mRNA":
		preproc_CDS_mRNA(input_file, ARGS.mrna, output_file)
	elif proc_type == "genome":
		preproc_CDS_bed(input_file, ARGS.genome, ARGS.bed, output_file, ARGS.list)
	elif proc_type == "filter":
		filtering_fasta(input_file, output_file)
	else:
		raise TypeError(
			"Type (-t|--type) should be ncRNA, mRNA or genome."
		)



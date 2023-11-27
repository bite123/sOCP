# _*_ coding: UTF-8 _*_
# Author: PENG Zhao
# Version: 2023-10
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
	max_len = 33
	for key,value in orf_dict.items():
		if value[1] >= max_len:
			max_len = value[1]
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
# 3) The sequence typeset is 70 nucleotides per line. 

def preproc_CDS_mRNA(input_cds, input_mrna, output_file):
	# parse CDS to get location
	pat = re.compile(r'location=(join\()?(\d+)')
	with open(input_cds) as f_in:
		pos_dict = {}
		i = 0
		for line in f_in:
			if line.startswith(">"):
				i += 1
				start_position = int(re.search(pat, line).group(2))
				pos_dict[i] = start_position

	# parse mRNA to get the upstreaming 3-nt
	with open(input_mrna) as f_in:
		up_dict = {}
		i = 0
		for line in f_in:
			if line.startswith(">"):
				i += 1
				start_position = pos_dict[i]
				if start_position <= 3:
					line = f_in.readline()
					up_site = "N" * (3-start_position+1) + line[0:start_position-1]
				else:
					jump_line = start_position//70
					residual_site = start_position%70
					if residual_site <= 3:
						jump_line -= 1
						residual_site += 70
						for j in range(jump_line):
							line = f_in.readline()
						first_line = f_in.readline().strip()
						second_line = f_in.readline()
						total_line = first_line + second_line
						up_site = total_line[residual_site-4:residual_site-1]
					else:
						for j in range(jump_line):
							line = f_in.readline()
						total_line = f_in.readline()
						up_site = total_line[residual_site-4:residual_site-1]					
				up_dict[i] = up_site

	# generate output files
	with open(input_cds) as f_in, open(output_file,"w") as f_out:
		i = 0
		for line in f_in:
			if line.startswith(">"):
				i += 1
				f_out.write(line)
				line = f_in.readline()
				up_site = up_dict[i]
				new_line = up_site + line
				f_out.write(new_line)
			else:
				f_out.write(line)


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
	with open(input_cds) as f_in, open(output_file, "w") as f_out:
		for line in f_in:
			if line.startswith(">"):
				cds_id = line.lstrip(">").strip()
				f_out.write(line)
				line = f_in.readline()
				if cds_id in start_dict and line[0:3] == start_dict[cds_id]:
					new_line = up_dict[cds_id] + line
				else:
					new_line = "NNN" + line
				f_out.write(new_line)
			else:
				f_out.write(line)


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
	else:
		raise TypeError(
			"Type (-t|--type) should be ncRNA, mRNA or genome."
		)



# _*_ coding: UTF-8 _*_
# Author: PENG Zhao
# Version: 2023-10
""" Instruction
This script provides functions utilized in other scripts.
"""

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y",                        
    "TGT": "C", "TGC": "C",             "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "TAG": "*", "TGA": "*", "TAA": "*"
}

# This function accepts a fasta file, and returns a dictionary, with record line as key, and sequences as value.
# If set the paramater parse False, the sequence will remain as list of lines.
# If set the paramater sep, the record line will be separated, and the first part will be used as the key. 
# Lower case in the sequences will be transferred to upper by default.
def fasta_to_dict(fasta, parse=True, sep=""):
    seq_dict = {}
    with open(fasta) as f_in:
        for line in f_in:
            if line.startswith("#"):
                continue
            elif line.startswith(">"):
                record = line
                if sep:
                    record = line.strip().lstrip(">").split(sep)[0]
                seq_dict[record] = []
            else:
                new_line = line.upper()
                seq_dict[record].append(new_line)
        if parse:
            new_dict = {}
            for k,v in seq_dict.items():
                new_v = ""
                for l in v:
                    new_v += l.strip()
                new_dict[k] = new_v
            seq_dict = new_dict
        return seq_dict

# This function accepts a PARSED but NO-SEP result form fasta_to_dict(),
# and transforms it back into fasta.
def dict_to_fasta(input_dict, output_file, mode="w"):
    with open(output_file, mode) as f_out:
        for k,v in input_dict.items():
            f_out.write(k)
            f_out.write(v + "\n")

# This function accepts a list of dictionaries with integer values, and returns a merge one.
# For a key in multiple dictionaries, the values will be summed. 
def dict_sum(dict_list):
    res_dict = {}
    for d in dict_list:
        for k,v in d.items():
            if k not in res_dict:
                res_dict[k] = v
            else:
                res_dict[k] += v
    return res_dict

# This function accepts a dictionary, and saves it in a tsv file. 
# Parameter sep is used if the value is a list.
def dict_to_tsv(input_dict, output_file, sep=""):
    with open(output_file, "w") as f_out:
        for k,v in input_dict.items():
            if sep:
                new_v = []
                for i in v:
                    new_v.append(str(i))
                new_v = sep.join(new_v)
            else:
                new_v = str(v)
            line = str(k) + "\t" + new_v + "\n"
            f_out.write(line)

# This function accepts a tsv file, and generates a dictionary.
# Parameter sep is set True if there are two or more columns as value.
def tsv_to_dict(input_file, sep=False):
    output_dict = {}
    with open(input_file) as f_in:
        for line in f_in:
            elements = line.strip("\n").split("\t")
            key = elements[0]
            if sep:
                value = elements[1:]
            else:
                value = elements[1]
            output_dict[key] = value
    return output_dict

# This fuction load hexamer list in Data directory
def load_hexamer_list():
    hexa_list_file = "Data/Hexamer.list"
    hexa_list = []
    with open(hexa_list_file) as f_in:
        for line in f_in:
            hexa_list.append(line.strip())
    return hexa_list

# These two fuctions load a feature list. The total feature list is provided in Data directory.
def load_feature_list(input_file):
    feature_list = []
    with open(input_file) as f_in:
        for line in f_in:
            elements = line.strip().split("\t")
            start = int(elements[0])
            end = int(elements[1])
            feature_list.extend(range(start,end))
    return feature_list

def load_feature_list_to_dictionary(input_file):
    feature_dict = {}
    with open(input_file) as f_in:
        for line in f_in:
            elements = line.strip().split("\t")
            f_id = elements[2]
            start = int(elements[0])
            end = int(elements[1])
            feature_dict[f_id] = list(range(start,end))
    return feature_dict


# _*_ coding: UTF-8 _*_
# Author: PENG Zhao
# Version: 2023-10
""" Instruction
This script predicts coding potential of smORFs using a pre-trained model.
""" 
    
import argparse
from joblib import dump, load
from FeatureDataProc import data_dict_read
from FeatureExtraction import Ntseq, feature_extraction
import Util

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Predicting coding potential")
    parser.add_argument("-i", "--input",
                        help="Input fasta file", required=True)
    parser.add_argument("-o", "--output",
                        help="Output prefix", required=True)
    parser.add_argument("-m", "--model",
                        help="Directory containing model files", required=True)
    ARGS = parser.parse_args()

    ################################################################################
    # Processing input and output file
    ################################################################################
    seq_dict = Util.fasta_to_dict(ARGS.input, parse=True, sep="")
    output_feature_file = ARGS.output + ".feature.tsv"
    output_score_file = ARGS.output + ".score.tsv"

    ################################################################################
    # Processing model files
    ################################################################################
    codon_in = ARGS.model + "/codonbias.tsv"
    codon_dict = data_dict_read(codon_in)
    hexa_in = ARGS.model + "/hexamer.tsv"
    hexa_dict = data_dict_read(hexa_in)
    ntbias_in = ARGS.model + "/ntbias.tsv"
    ntbias_dict = data_dict_read(ntbias_in)

    f_list_in = ARGS.model + "/feature.list"
    feature_list = Util.load_feature_list(f_list_in)

    model_in = ARGS.model + "/pretrained.model"
    clf = load(model_in)

    ################################################################################
    # Predicting
    ################################################################################
    feature_res = {}
    score_res = {}
    for key,value in seq_dict.items():
        orf_id = key.lstrip(">").strip()
        orf_feature = feature_extraction(value, codon_dict, hexa_dict, ntbias_dict, feature_list)
        feature_res[orf_id] = orf_feature
        # predict
        orf_feature = [orf_feature]
        score = clf.predict_proba(orf_feature)[0][1]
        score_res[orf_id] = score

    ################################################################################
    # Generating files
    ################################################################################    
    Util.dict_to_tsv(feature_res, output_feature_file, sep="\t")
    Util.dict_to_tsv(score_res, output_score_file, sep=False)




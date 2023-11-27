# sOCP

> Predicting coding potential of smORFs

# Usage

The sOCP software consists of 3 functions:

1. Predicting coding potential of smORFs in human, using a pre-trained model.
2. Providing a pipeline to train a model for predicting coding potential of smORFs.
3. Providing a pipeline to predict smORFs from the genome using a prepared model。

## 1 Getting started

### 1.1 Download 

```
git clone https://github.com/bite123/sOCP.git
cd sOCP
```

### 1.2 Hardware and Software Requirements

The memory requirement is XXX

sOCP requires python version 3.10 or higher, and packages XXXX. 

## 2 Predicting coding potential of smORFs in human (Eg01)

The pre-trained model for human is provided in **sOCP/HumanModel**, including:

- **codonbias.tsv**
- **feature.list**
- **hexamer.tsv**
- **ntbias.tsv**
- **pretrained.model**

> As seen, running sOCP model requires not only the model itself, but also a list of features and three TSV files related to codonbias, hexamer, and ntbias. As the directory of model could be any word, these files are universally named. 

Taking the FASTA file **Examples/Eg01.fasta** as an example. The command lines are as below:

```
python PredictingCodingPotential.py \
-i ./Examples/Eg01.fasta \
-o Eg01 \
-m ./HumanModel
```

- **-i** The path to the input FASTA file.
- **-o** The prefix for output files.
- **-m** The path to the model directory.

The output file **Eg01.feature.tsv** and **Eg01.score.tsv** will be generated, recording feature information and predicting scores for input sequences separately.

## 3 Preprocessing sequences for predicting

The input sequence for sOCP includes the complete smORF and the upstreaming three nucleotides, named extended ORF in sOCP. 

However, common ORF sequences lack the upstreaming three nucleotides. sOCP provides scripts to extract extended ORFs from mRNA and ncRNA. 

There are three modes of preprocessing methods.

### 3.1 Extracting extended ORFs from ncRNA (Eg02)

```
Corresponding description in PreprocessingFasta.py
# This function accepts a SHORT nucleotide sequence, 
# and returns the extended ORF with the maximum length within the sequence.
# UPPER CASE IS NEEDED!
# SHORT indicates the sequence is commonly an RNA, not a chromosome.
# For extended ORF, if there is an incomplete upstreaming 3-nt, using N instead.
# If no ORF >= 11 codons is found, return "".
# Also returning length, start position and end position.
```

The command lines are as below:

```
python PreprocessingFasta.py \
-t ncRNA \
-i ./Examples/Eg02.ncrna.fasta \
-o Eg02.ext_orf.fasta
```

- **-t** The type of the input sequences.
- **-i** The path to the input FASTA file.
- **-o** The path to the output FASTA file.

The output file **Eg02.ext_orf.fasta** will be generated. Noted that of 10 input sequences, only in 7 sequences an ORF is found.

### 3.2 Extracting extended ORFs based on mRNA (Eg03)

```
Corresponding description in PreprocessingFasta.py
# This function extracts the upsteaming 3-nt, based on a CDS file and its corresponding mRNA file.
# The two fasta files can be retrieved from NCBI RefSeq, and will have the following formats:
# 1) The CDS and mRNA are in one-to-one correspondence in order.
# 2) There is location information in the header line, e.g.:
# >lcl|NM_001419531.1_prot_NP_001406460.1_1 ... [location=256..4701] ...
# 3) The sequence typeset is 70 nucleotides per line. 
```

The command lines are as below:

```
python PreprocessingFasta.py \
-t mRNA \
-i ./Examples/Eg03.orf.fasta \
-m ./Examples/Eg03.mrna.fasta \
-o Eg03.ext_orf.fasta
```

- **-t** The type of the input sequences.
- **-i** The path to the input ORF FASTA file.
- **-m** The path to the mRNA FASTA file. 
- **-o** The path to the output FASTA file.

The output file **Eg03.ext_orf.fasta** will be generated. 

### 3.3 Extracting extended ORFs based on genome and bed annotation (Eg04)

```
Corresponding description in PreprocessingFasta.py
# This function extracts the upsteaming 3-nt, based on a CDS file and its corresponding bed.
# The genome fasta file is required as well.
# If the chromosome names in the first column of bed file do not match
# the header line of the genome chromosomes, a tsv file as follow is required:
# Column 1: Chromosome ID in genome fasta; Column 2: Chromosome ID in bed 
```

The command lines are as below:

```
python PreprocessingFasta.py \
-t genome \
-i ./Examples/Eg04.orf.fasta \
-g ./Examples/Eg04.genome.fasta \
-b ./Examples/Eg04.bed \
-l ./Examples/Eg04.chr_mapping.list \
-o Eg04.ext_orf.fasta
```

- **-t** The type of the input sequences.
- **-i** The path to the input ORF FASTA file.
- **-g** The path to the genome FASTA file. 
- **-b** The path to the annotation BED file. 
- **-l** The path to the file mapping the chromosome IDs between BED and genome FASTA. 
- **-o** The path to the output FASTA file.

The output file **Eg04.ext_orf.fasta** will be generated. 

## 4 Preparing dataset (Eg05)

In **Section 4** to **Section** **8**, there will be a series of methods showing how to train a model for predicting coding potential of smORFs.

First of all, the positive dataset (extended ORFs from mRNA) and the negative dataset (extended ORFs from mRNA) are required. You can retrieve these FASTA files from NCBI and Ensemble database, and preprocess them to get extended ORFs, as shown in Section 3. Remove redundant sequences between the positive dataset and negative dataset using the CD-HIT software is also recommended.

In this section, the positive dataset and the negative dataset will be divided for training and testing separately. The training part will be divided further for cross validation.

The command lines are as below:

```
python DatasetProc.py \
-p ./Examples/Eg05.pos.fasta \
-n ./Examples/Eg05.neg.fasta \
-b 1 \
-r 2 \
-f 10 \
-o Eg05
```

- **-p** The path to the positive dataset. to the positive dataset (extended ORFs from mRNA).
- **-n** The path to the positive dataset. to the negative dataset (extended ORFs from ncRNA).
- **-b** The ratio of the sequence numbers between the positive and negative. 
- **-r** The ratio of the sequence numbers between training and testing.  
- **-f** Folds of cross validation.
- **-o** The prefix for output files.

The output file includes:

- **Eg05.filtering.fasta**: Sequences filtering from the input FASTA, caused by out of frame, containing non-ACGT nucleotides, containing no stop codon, and less than too short (less than 11 codons).
- **Eg05.pos.training.fasta**, **Eg05.pos.testing.fasta**: Divided positive dataset.
- **Eg05.neg.training.fasta**, **Eg05.neg.testing.fasta**: Divided negative dataset.
- **Eg05.cv_file_list.tsv**: List of cross validation files.
- **Eg05.cv01** to **Eg05.cv09** files: 10-fold cross validation, four files in each (positive training, positive testing, negative training, negative testing).

## 5 Preparing for data-dependent-features (Eg06)

Codon bias, hexamer score and nucleotide bias in sOCP are called data-dependent-features. To extract these three features, three TSV files integrating positive and negative information should be prepared first.

The total training data (not single subset of cross validation folds) from **Section 4 Eg05**will be used. The command lines are as below:

```
python FeatureDataProc.py \
-p ./Examples/Eg05_output/Eg05.pos.training.fasta \
-n ./Examples/Eg05_output/Eg05.neg.training.fasta \
-o Eg06
```

- **-p** The path to the positive dataset.
- **-n** The path to the negative dataset.
- **-o** The prefix for output files.

The output file includes **Eg06.codonbias.tsv**, **Eg06.hexamer.tsv**, and **Eg06.ntbias.tsv**.

***Note*:** 

- *Some sOCP commands will generate "tmp_" files simultaneously. These files are useless henceforth and could be deleted.*
- *The count of missing codon and hexamer will be set as 1 to avoid  divided-by-zero. However, this will introduce many results of "count 1" when the dataset is small, such as in **Eg06.hexamer.tsv**.*

## 6 Model selection (Eg07)

After acquiring the cross validation datasets, sOCP provides a method to select a model from the sci-kit learn package (sklearn) with a good performance.

First, a list of candidate models should be prepared. The sklearn codes of some models are listed in **Data/Models.list**.  Any code available in sklearn could be used as well.

In this example we set **Eg07.model.list** as:

```
svm.SVC(kernel="rbf")
ensemble.RandomForestClassifier(max_features="sqrt")
ensemble.GradientBoostingClassifier(learning_rate=0.8)
```

A list of input FASTA files is also required, recording paths to the cross validation FASTA files generated in **Section 4 Eg05**.

The TSV files generated in **Section 5 Eg06** is also requred. They are moved to **Eg07_dict** and renamed to the universal names (**codonbias.tsv**, **hexamer.tsv** and **ntbias.tsv**).

At last, we need an sOCP feature list.

### sOCP feature tips

> In Data directory **FeatureMapping.list** and **DetailedFeatureMapping.list** show what format an sOCP feature list should be:
>
> - The first and second column indicate the start and end location a & b, that is, a list as `[a, ..., b-1]`. In other words, the length of the feature list equals b-a. Actually, the list could acquired by the python command `list(range(a,b))`.  
> - The third column is the ID column. It will only be used in a few steps. And the ID could be set as any word as you like.
>
> sOCP provides a complete feature set as a 1169-d vector (the location 0 of the list is ORF ID and not considered as a feature). However, a complicate feature set is not recommended if your training dataset is small.

In this example we set **Eg07.feature.list** as:

```
0       1       orf_id
1       2       orf_length
2       6       fickett
6       7       hexamer
7       8       codon_bias
8       14      nt_bias
14      18      1_mer
18      34      2_mer
```

The command lines are as below:

```
python ScikitModel.py \
-c ./Examples/Eg07.cv_file_list.tsv \
-d ./Examples/Eg07_dict \
-m ./Examples/Eg07.model.list \
-f ./Examples/Eg07.feature.list \
-o Eg07
```

- **-c** The path to the list of cross validation files.
- **-d** The directory containing codonbias, hexamer and ntbias files.
- **-m** The path to the list of models.
- **-f** The path to the list of features.
- **-o** The prefix for output files.



## 7 Feature selection

### 7.1 Feature ablation

### 7.2 mRMR-IFS

## 8 Training model and performance evaluation

### 8.1 Model training and application

### 8.2 Performance evaluation

## 9 Predicting smORFs from the genome

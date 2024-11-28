# sOCP

> Predicting coding potential of smORFs

The sOCP software consists of 3 functions:

1. Predicting coding potential of smORFs in human, using a pre-trained model.
2. Providing a pipeline to train a model for predicting coding potential of smORFs.
3. Providing a method to predict smORFs from the genome using a prepared model。

**Updated history:**

- 2024 Nov, code optimization in PreprocessingFasta.py as below: 
  - Fix a bug in extraction of extended ORFs from ncRNA.
  - Extraction of extended ORFs from mRNA will no longer be limited to input FASTA files with 70 words one line.
  - Output sequences will be merged into one line.
  - The regular expression pattern in extraction of extended ORFs from mRNA is changed to accept more forms.
  - The "filter" mode is added.

# Basic function

## 1 Getting started

### 1.1 Download 

```
git clone https://github.com/bite123/sOCP.git
cd sOCP
```

### 1.2 Dependencies

Running sOCP requires python version 3.10 or higher, and packages Bio, joblib, pandas, pymrmr, and sklearn.

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

The output files **Eg01.feature.tsv** and **Eg01.score.tsv** will be generated, recording feature information and predicting scores for input sequences separately.

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
# Only ORFs with 11 - 101 codons (including the stop codon) will be retained
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

***Note*:** 

- *This function is generally to process mRNA sequences from NCBI RefSeq or other database, in which both mRNA and ORF within the mRNA are provided.*
- *If you have the mRNA sequence but have no idea where the ORF starts and ends, you can only use the function in **3.1** instead, which require merely the RNA sequence as input. However, it is not recommended.*

```
Corresponding description in PreprocessingFasta.py
# This function extracts the upsteaming 3-nt, based on a CDS file and its corresponding mRNA file.
# The two fasta files can be retrieved from NCBI RefSeq, and will have the following formats:
# 1) The CDS and mRNA are in one-to-one correspondence in order.
# 2) There is location information in the header line, e.g.:
# >lcl|NM_001419531.1_prot_NP_001406460.1_1 ... [location=256..4701] ... 
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

## 3.4 Filtering extended ORFs

```
Corresponding description in PreprocessingFasta.py
# This function accepts extended ORFs, and removed those:
# 1) not multiples of 3
# 2) have non-ACTG nucleotides within
# 3) have early stop codons 
# 4) have an out-of-range length (not within 11-101 codons, including the stop codon)
```

Note that extraction in **3.2** and **3.3** only generate extended ORFs according to the location and will not consider if the ORF sequence is suitable for training.  

If you want to filter extended ORFs according to the above criteria, use the "filter" mode as below:

```
python PreprocessingFasta.py \
-t filter \
-i input.fasta \
-o output.fasta
```

# Training a model to predict coding potential of smORFs

In **Section 4** to **Section** **8**, there will be a series of methods showing how to train a model for predicting coding potential of smORFs.

## 4 Preparing dataset (Eg05)

First of all, the positive dataset (extended ORFs from mRNA) and the negative dataset (extended ORFs from ncRNA) are required. You can retrieve these FASTA files from NCBI and Ensemble database, and preprocess them to get extended ORFs, as shown in **Section 3**. 

**Note:**

- You can filter the dataset according to ORF length manually or using the "filter" mode as **3.4** described.

- Remove redundant sequences between the positive dataset and negative dataset, or within each dataset using the CD-HIT software is also recommended.

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

The output files include:

- **Eg05.filtering.fasta**: Sequences filtering from the input FASTA, caused by out of frame, containing non-ACGT nucleotides, containing no stop codon, and less than too short (less than 11 codons).
- **Eg05.pos.training.fasta**, **Eg05.pos.testing.fasta**: Divided positive dataset.
- **Eg05.neg.training.fasta**, **Eg05.neg.testing.fasta**: Divided negative dataset.
- **Eg05.cv_file_list.tsv**: List of cross validation files.
- **Eg05.cv01** to **Eg05.cv09** files: 10-fold cross validation, four files in each (positive training, positive testing, negative training, negative testing).

## 5 Preparing for data-dependent-features (Eg06)

Codon bias, hexamer score and nucleotide bias in sOCP are called data-dependent-features. To extract these three features from an extended ORF, three TSV files integrating positive and negative information should be prepared first.

The total training data (not single subset of cross validation folds) from **Section 4 Eg05** will be used. The command lines are as below:

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

First, a list of candidate models should be prepared. The sklearn codes of some models are listed in **Data/Models.list**.  Any model code available in sklearn could be used here as well.

In this example we set **Eg07.model.list** as:

```
svm.SVC(kernel="rbf")
ensemble.RandomForestClassifier(max_features="sqrt")
ensemble.GradientBoostingClassifier(learning_rate=0.8)
```

***Note*:** 

- *For any other running step of sOCP except Model selection, the model has been determined, so only the first line of model.list will be considered.*

Second, we need a list of input FASTA files recording paths to the cross validation FASTA files, with the format as:

```
Corresponding description in ScikitModel.py
''' The format of Cross-validation feature file list is a tsv file as below:
cv/cv01_pos_train    cv/cv01_pos_test    cv/cv01_neg_train    cv/cv01_neg_test
cv/cv02_pos_train    cv/cv02_pos_test    cv/cv02_neg_train    cv/cv02_neg_test
... ...
cv/cv10_pos_train    cv/cv10_pos_test    cv/cv10_neg_train    cv/cv10_neg_test
'''
```

Here we used the cross validation files generated in **Section 4 Eg05**.

Third, the TSV files generated in **Section 5 Eg06** are also requred. They are moved to **Eg07_dict** and renamed to the universal names (**codonbias.tsv**, **hexamer.tsv** and **ntbias.tsv**).

At last, we need an sOCP feature list.

### Important: features in sOCP

> In Data directory **FeatureMapping.list** and **DetailedFeatureMapping.list** show what format an sOCP feature list should be:
>
> - The first and second column of a feature indicate the start and end location a & b, that is, locations including `[a, ..., b-1]`. In other words, the length of the feature equals b-a. Actually, the locations could acquired by the python function `list(range(a,b))`.  
> - The third column is the ID of the feature. It will only be used in a few steps. And the ID could be set as any word as you like.
> - You could create  a new feature line like `1	33	base`，`6	8	hexamer_codonbias`, only to be aware not exceeding 1170.
>
> sOCP provides a complete feature set as a 1169-d vector (the location 0 of the list is ORF ID and not considered as a feature). However, a complicate feature set is not recommended if your training dataset is small.
>
> Initially the **base** feature subset `1	33	base` (including ORF ID, ORF length, Fickett score, Hexamer score, Codon bias, nucleotide bias, 1-mer, and 2-mer) is a good choise for a training dataset with thousands of sequences. Other k-mer, g-gap, and g-bigap could be added to feature list if you have a huge dataset or the base feature subset performed not well.  

In this example we set **Eg07.feature.list** as below. Only a 13-d subset excluding `1	2	orf_lenth`  is used, because the example dataset is too small and will not  tell difference between models if using ORF length:

```
2       6       fickett
6       7       hexamer
7       8       codon_bias
8       14      nt_bias
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

The output files include **Eg07.model_acc.tsv** and **Eg07.model_mcc.tsv**, recording the evaluation scores ACC and MCC for the candidate models (each column represents one fold of cross validation). You could determine which model will be used in the following steps according to the scores.

## 7 Feature selection

The initial feature list could be optimized by removing features useless for model training, that is, feature selection. There are two methods of  feature selection in sOCP: feature ablation and mRMR-IFS.

### 7.1 Feature ablation (Eg08)

In feature ablation, from the initial feature list, one feature line is excluded at one time to form a new subset. Then all new subsets is evaluated by cross validation, to show if any feature line in the initial feature list could removed.

For example, if an ablation experiment is applied to **Eg07.feature.list**:

```
2       6       fickett
6       7       hexamer
7       8       codon_bias
8       14      nt_bias
```

One feature line is excluded at one time, resulting in four subsets, each with three features.

We used the same cross validation files and TSV-containing directory as in **Eg07**. However, given the model has been determined, **Eg08.model.list** is only one line as:

```
ensemble.RandomForestClassifier(max_features="sqrt")
```

The command lines are as below:

```
python FeatureAblation.py \
-c ./Examples/Eg07.cv_file_list.tsv \
-d ./Examples/Eg07_dict \
-m ./Examples/Eg08.model.list \
-f ./Examples/Eg07.feature.list \
-o Eg08
```

- **-c** The path to the list of cross validation files.
- **-d** The directory containing codonbias, hexamer and ntbias files.
- **-m** The path to the file recording the model.
- **-f** The path to the list of features.
- **-o** The prefix for output files.

The output files include **Eg08.ablation_acc.tsv** and **Eg08.ablation_mcc.tsv**, recording the evaluation scores ACC and MCC for the ablation subsets (each column represents one fold of cross validation). You could determine if one of the features should be removed from the initial features.

The feature ablation could be repeated until satisfied. For example, we determine to remove codon bias from the initial feature list based on **Eg08** output files, then a further ablation could be processed on the remaining feature list:

```
2       6       fickett
6       7       hexamer
8       14      nt_bias
```

### 7.2 mRMR-IFS (Eg09)

The mRMR-IFS method ranks all features, incrementally generates the feature subsets consisting of the top-ranked *N* features one-by-one, and tells how these feature subsets perform in model training. 

Note that here the *n*-dim feature is divided into scalars, e.g., nt_bias is divided into 6 scalars in this step. See **Data/DetailedFeatureMapping.list**. 

Therefore mRMR-IFS is commonly applied for a small feature list, and a complicate feature list is better to be shortened by feature ablation first.

To show the difference between the two feature selection methods, we applied the same input files as in **Eg08**:

```
python FeatureIFS.py \
-c ./Examples/Eg07.cv_file_list.tsv \
-d ./Examples/Eg07_dict \
-m ./Examples/Eg08.model.list \
-f ./Examples/Eg07.feature.list \
-o Eg09
```

The output files **Eg09.ablation_rank.tsv** indicates the mRMR rank order (each column for one fold of cross validation), and the incremental feature subsets have generated based on the order. For example, the first column `6,10,12,3,2,11,5,9,4,7,8,13` indicates that, in the first fold of cross validation, the incremental subsets will be `[6], [6,10], ... [6,10,12,3,2,11,5,9,4,7,8,13]`, twelve subsets in total. The numbers  are feature IDs according to the input feature list. See details in **Data/DetailedFeatureMapping.list**. 

The output files **Eg09.ablation_acc.tsv** and **Eg09.ablation_mcc.tsv** records the evaluation scores ACC and MCC, each column for one fold of cross validation, each line for an incremental subset, according to **Eg09.ablation_rank.tsv**. You could determine if a feature subset is chosen out of the initial features.

## 8 Training model (Eg10)

Once model and feature list have been determined by cross validation, at last the model will be trained and saved for further prediction.

The data-dependent-feature directory, model, and feature list are described above (using the initial feature dataset). For the training and testing files, the FASTA files generated in **Eg05** will be used.

```
Examples/Eg05_output/Eg05.pos.training.fasta
Examples/Eg05_output/Eg05.pos.testing.fasta
Examples/Eg05_output/Eg05.neg.training.fasta
Examples/Eg05_output/Eg05.neg.testing.fasta
```

The command lines are as below.:

```
python ModelTraining.py \
-i ./Examples/Eg10.model_training.list \
-d ./Examples/Eg07_dict \
-m ./Examples/Eg08.model.list \
-f ./Examples/Eg07.feature.list \
-o Eg10_model
```

- **-i** The path to the list of four FASTA files: positive training, positive testing, negative training, and negative testing.
- **-d** The directory containing codonbias, hexamer and ntbias files.
- **-m** The path to the file recording the model.
- **-f** The path to the list of features.
- **-o** The output directory.

The output directory **Eg10_model** will be generated. As described in **Section 2**, there are five universally named files for a trained sOCP model. The trained model could now be applied for predicting using  **PredictingCodingPotential.py**, as seen in **Eg11**.

## 9 Performance evaluation (Eg11)

As ACC and MCC has been printed once model training is complete, sOCP provides functions to explore more performance scores of model.

First, **Eg05.pos.testing.fasta** and **Eg05.neg.testing.fasta** are concatenated to generate **Eg11.fasta**.

Then, use the model trained in **Eg10** to get prediction result as in Eg01 :

```
python PredictingCodingPotential.py \
-i ./Examples/Eg11.fasta \
-o Eg11 \
-m ./Examples/Eg10_model
```

The output files **Eg11.feature.tsv** and **Eg11.score.tsv** will be generated.

For performance evaluation, a label file **Eg11.label.tsv** is required as well. There are two columns in **Eg11.label.tsv**, and the first column indicates the ORF ID, the same with **Eg11.score.tsv**.  The second column of **Eg11.label.tsv** is 1 or 0, indicates whether the ORF is from positive or negative dataset.

> Here is an example command showing how to generate the label file from the score file. Number 34 indicates in **Eg11.score.tsv**, the first 34 ORFs are positive.

>```
>awk 'BEGIN{FS="\t";OFS="\t"}NR<=34{print $1,1}NR>34{print $1,0}' Eg11.score.tsv > Eg11.label.tsv 
>```

At last, the performance evaluation is processed as:

```
python PerformanceEvaluation.py \
-i ./Examples/Eg11.score.tsv \
-a ./Examples/Eg11.label.tsv \
-o Eg11
```

The evaluation metrics will be printed to screen at once. There will be these output files for details:

- **Eg11.preprocessed.tsv**, a pre-preprocessed file integrating score and label files together.
- **Eg11.metric.tsv**. Since the printed metrics are calculated for the default cutoff 0.5, this file records metrics (TPR, TNR, FPR, FNR, PRE, ACC, F1S, HM, and MCC) according to other cutoffs.
- **Eg11.roc.tsv**, providing data for ROC curve plotting.
- **Eg11.pr.tsv**, providing data for PR curve plotting.

Of course, **PredictingCodingPotential.py** could help you to acquire evaluation metrics from any other source, as long as the score (ranged from 0 to 1) and label files are provided.

# Predicting smORFs from genome

In **Section 10** and **Section** **11**, we will show how to use a trained model to predict smORFs from genome.

## 10 Set cutoff for prediction

Before prediction, a cutoff list is required. See **Data/CutoffConf.list**:

```
ATG	0.5
CTG	0.75
GTG	1
ACG	1
TTG	0
ATT	0
ATC	0
ATA	0
AAG	0
AGG	0
```

The first column indicates the start codon of smORFs, and the second column indicates the score cutoffs. Only smORFs with scores higher than its corresponding cutoff will be reserved during prediction. Cutoff 0 means smORFs with that start codon will not be considered.

Non-ATG start codons vary largely in efficiency and frequency. You could adjust their cutoff to your interest.  Or set all cutoffs of non-ATG start codons as 0, if you only concern about ATG-started smORFs.

## 11 Predicting smORFs from the genome (Eg12)

As an example, we will use the model trained in **Eg10**. And a ~700k sequence from human genome will be used as the input genome, split into ChrA and ChrB.

The command lines are as below:

```
python PredictingOnGenome.py \
-i ./Examples/Eg12.genome.fasta \
-m ./Examples/Eg10_model \
-c ./Examples/Eg12.cutoff.list \
-o Eg12
```

- **-i** The path to the input genome FASTA file.
- **-m** The path to the model directory.
- **-c** The path to the list of cutoffs.
- **-o** The prefix for output files.
- **-l** The minimum number of codons (excluding stop codon).  This parameter is optional and set as 10 by default.
- **-L** The maximum number of codons (excluding stop codon). This parameter is optional and set as 100 by default.

The predicted smORF ID ChrA#1213#+#CTG#760 indicates the smORF starts at Position 1213 on Chromosome ChrA, Strand + (forward), with CTG as its start codon, and is scored 0.760 by sOCP. 

The output files includes:

- **Eg12.socp.fasta**, recording the predicted smORFs.
- **Eg12.socp.bed**, recording the basic information of the smORFs.
- **Eg12.socp.upsite.tsv**, recording the the upstreaming three nucleotides of smORFs, as the nucleotides have been used in prediction.

It usually takes a long time to predict smORFs from genome. While multi-processing is not available in sOCP, you might split the genome FASTA into small files, (e. g. one chromosome in one file) and run sOCP In parallel.

# Contact

For more details, please see the paper:

**sOCP: a framework predicting smORF coding potential based on TIS and in-frame features, and effectively applicated in the human genome**

If you have any problem, please contact Peng Zhao (bite123jz@163.com).

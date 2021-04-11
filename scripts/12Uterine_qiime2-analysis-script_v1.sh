#!/bin/bash

## Script for Qiime2 analysis of Uterus Microbiome Project
## V1: 4.5.2020


#Preparations
#FASTQ files were loaded on a AWS x2large instance with 150 gb storage
#Manifest file was created by linux find command written to csv and ms excel


## 1-Join F- and R file and denoise using DADA2
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path qiime2/uterus_manifestfile.txt \
  --output-path qiime2/uterus-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2
	
  
   qiime dada2 denoise-paired \
  --i-demultiplexed-seqs qiime2/uterus-demux.qza \
  --p-trim-left-f 13 \
  --p-trim-left-r 13 \
  --p-trunc-len-f 250 \
  --p-trunc-len-r 240 \
  --p-max-ee-f 3 \
  --p-max-ee-r 4 \
  --p-n-threads 0\
  --o-table qiime2/table.qza \
  --o-representative-sequences qiime2/rep-seqs.qza \
  --o-denoising-stats qiime2/denoising-stats.qza
  
qiime feature-table summarize \
  --i-table qiime2/table.qza \
  --o-visualization qiime2/table.qzv \
  --m-sample-metadata-file uterine_metadata.txt


### Train classifier### NAIVE Bayesian classifier for Taxonomic assignment, based on Greengeenes 13.8 99% OTUs
echo "/n Start extracting reads" 

qiime feature-classifier extract-reads \
  --i-sequences 99_otus.qza \
  --p-f-primer AGAGTTTGATYMTGGCTCAG \
  --p-r-primer ATTACCGCGGCTGCTGG \
  --p-min-length 300 \
  --p-max-length 600 \
  --o-reads ref-seqs.qza

echo "/n End extracting reads"

echo "/n Train classifier - START"
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ref-seqs.qza \
  --i-reference-taxonomy ref-taxonomy.qza \
  --o-classifier classifier.qza
  
echo "/n Train classifier END"
  
### Taxonomic classification
qiime feature-classifier classify-sklearn \
  --i-classifier training-feature-classifier/classifier.qza \
  --i-reads qiime2/rep-seqs.qza \
  --o-classification qiime2/taxonomy.qza

qiime feature-table filter-samples \
  --i-table qiime2/table.qza \
  --m-metadata-file uterine_metadata.txt \
  --p-where "[Group] IN ('RIF', 'RPL', 'CTRL')" \
  --o-filtered-table qiime2/table_filt.qza
    
 qiime tools export qiime2/table_filt.qza --output-dir qiime2
 
 qiime tools export qiime2/taxonomy.qza --output-dir qiime2
  
  
  
  ## Make Tree
 qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences qiime2/rep-seqs.qza \
  --o-alignment qiime2/aligned-rep-seqs.qza \
  --o-masked-alignment qiime2/masked-aligned-rep-seqs.qza \
  --o-tree qiime2/unrooted-tree.qza \
  --o-rooted-tree qiime2/rooted-tree.qza
  
qiime tools export \
  unrooted-tree.qza \
  --output-dir qiime2

  
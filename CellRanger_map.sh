#!/bin/bash

## Cell Ranger
 
# mapping mouse samples using cellranger
 
module load cellranger/7.1.0
 
# running counting on SRR16409869-SRR16409878 sample
cellranger count --id=SRR16409869_78 \
--fastqs=/home/projects/misc/evoec/data/mouse_ec/SRR16409869_78 \
--transcriptome=/home/databases/genomes/processed/Mus_mus_38.primary/indexes/genome/cellranger/refdata-gex-mm10-2020-A
 
 
# mapping human ec  samples using cellranger
 
module load cellranger/7.1.0
 
# running counting on SRR16928876 sample
cellranger count --id=SRR16928876 \
--fastqs=/home/projects/misc/evoec/data/humanec/SRR16928876 \
--transcriptome=/home/databases/genomes/processed/Hom_sap_38.primary/indexes/genome/cellranger/refdata-gex-GRCh38-2020-A
 
 
 
# running CellRanger on 10X snRNA-seq pig samples using refseq reference files
 
module load cellranger/7.1.0
 
# filtering gtf file
cellranger mkgtf \
/home/projects/misc/evoec/development/sbb_mapping/cellranger/pig/ref_genome/genomic.gtf \
/home/projects/misc/evoec/development/sbb_mapping/cellranger/pig/ref_genome/filt_genome.gtf \
--attribute=gene_biotype:protein_coding
 
# building the index using cellranger filtered gtf
cellranger mkref \
  --genome=Sus_genome \
  --fasta=/home/projects/misc/evoec/development/sbb_mapping/cellranger/pig/ref_genome/GCF_000003025.6_Sscrofa11.1_genomic.fna \
  --genes=/home/projects/misc/evoec/development/sbb_mapping/cellranger/pig/ref_genome/filt_genome.gtf
 
 
 
# mapping pig  ec  samples using cellranger
 
module load cellranger/7.1.0
 
# running counting on pig adult  sample
cellranger count --id=Pig_adult2 \
--fastqs=/home/projects/misc/evoec/data/pigdata \
--transcriptome=/home/projects/misc/evoec/development/sbb_mapping/cellranger/pig/indexing/Sus_genome \
--chemistry=SC3Pv2
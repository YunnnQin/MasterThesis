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


## STARsolo

# Mapping human samples with STARsolo

# Map to  human reference genome GRCh38 (Ensembl release 98)
# Create genome index using STAR
# Download human genome reference 	Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz	
# Define the path to human genome
GenomeDir=/work/Reference/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa

# Download human gene annotation GTF  Homo_sapiens.GRCh38.109.gtf.gz
# Define the path to gene annotation
GTFfile=/work/Reference/Homo_sapiens.GRCh38.109.gtf

# Define the path to output directory
OutputDir=/work/Reference/ref_human

# Run STAR to get the genome index
STAR --runMode genomeGenerate \
     --genomeDir $OutputDir \
     --genomeFastaFiles $GenomeDir \
     --sjdbGTFfile $GTFfile \
     --runThreadN 16

# Map Sample HSB231_EC, HSB237_EC, HSB628_EC to reference 
# Define the path to the reference genome index
GENOME_INDEX=/work/mapping_files/Reference/ref_human


# Define the path to whitelist
# Barcode length is 28!
WHITELIST_DIR=/work/mapping_files/whitelist/3M-february-2018.txt

# Define the paths to the input file folders
INPUT_DIR_R1=/work/human_dataset/R1.fastq
INPUT_DIR_R2=/work/human_dataset/R2.fastq


# Iterate over the input files and run STAR on each pair
for R1_file in $INPUT_DIR_R1/*.fastq.gz
do
    # Define the input file pair
    R2_file=$INPUT_DIR_R2/$(basename "$R1_file" | sed 's/_1/_2/')

    # Define the output filename
    output=$(basename "$R1_file" _R1.fastq)

    # Run STAR on the input file pair
   STAR --genomeDir $GENOME_INDEX \
     --readFilesIn $R2_file $R1_file \
     --soloType CB_UMI_Simple \
     --soloCBwhitelist $WHITELIST_DIR \
     --soloUMIlen 12 \
     --soloFeatures Gene GeneFull SJ Velocyto \
     --soloOutFileNames human \
     --outSAMtype BAM SortedByCoordinate \
     --readFilesCommand gunzip -c

    # Print a message indicating that the mapping is complete
    echo "Mapping of $R1_file and $R2_file is complete"
done



## PIG DATASET

# Create Genome index using STAR

# Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa.gz 
# Genome reference
# Define the path to genome reference
GenomeDir=/work/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa
# Sus_scrofa.Sscrofa11.1.109.gtf.gz
# Gene annotation
# Define the path to gene annotation
GTFfile=/work/Sus_scrofa.Sscrofa11.1.109.gtf

# Define the path to output directory
OutputDir=/work/Reference/ref_pig

# Run STAR to get the genome index
STAR --runMode genomeGenerate \
     --genomeDir $OutputDir \
     --genomeFastaFiles $GenomeDir \
     --sjdbGTFfile $GTFfile \
     --runThreadN 16


# Merge 4 replicates for Adult_nc sample into one fastq file

# Merge R1
cat SRR9705091_1.fastq.gz SRR9705092_1.fastq.gz SRR9705093_1.fastq.gz SRR9705094_1.fastq.gz > merged_1.fastq.gz
# Merge R2
cat SRR9705091_2.fastq.gz SRR9705092_2.fastq.gz SRR9705093_2.fastq.gz SRR9705094_2.fastq.gz > merged_2.fastq.gz


# Mapping Sample Adult_nc to reference using STARsolo
# Define the path to the reference genome index
GENOME_INDEX=/work/mapping_files/Reference/ref_pig
# Define the path to the output directory
OUTPUT_DIR=/work/mapped_output
# Define the path to whitelist
# Barcode length = 26 
WHITELIST_DIR=/work/mapping_files/whitelist/737K-august-2016.txt 

# Define the paths to the input file folders
# Using merged fastq files from 4 replicates!
INPUT_DIR_R1=/work/pig_dataset/R1.fastq/merged_1.fastq.gz
INPUT_DIR_R2=/work/pig_dataset/R2.fastq/merged_2.fastq.gz

  # Run STAR on the input file pair
   STAR --genomeDir $GENOME_INDEX \
     --readFilesIn $INPUT_DIR_R2 $INPUT_DIR_R1 \
     --soloType CB_UMI_Simple \
     --soloCBwhitelist $WHITELIST_DIR \
     --soloFeatures Gene GeneFull SJ Velocyto \
     --soloOutFileNames output \
     --outSAMtype BAM SortedByCoordinate \
     --readFilesCommand gunzip -c

# Map to  mouse reference genome GRCm39 (mm10) Refseq (the latest version) NCBI
# Create genome index using STAR
# Download mouse genome reference 
# Define the path to human genome
GenomeDir=/work/Reference/HGCF_000001635.27_GRCm39_genomic.fna

# Download human gene annotation GTF 
# Define the path to gene annotation
GTFfile=/work/Reference/genomic.gtf

# Define the path to output directory
OutputDir=/work/Reference/ref_mouse

# Run STAR to get the genome index
STAR --runMode genomeGenerate \
     --genomeDir $OutputDir \
     --genomeFastaFiles $GenomeDir \
     --sjdbGTFfile $GTFfile \
     --runThreadN 16


# Mapping 9 mouse ENT Samples to reference using STARsolo

# Define input and output directories
INPUT_DIR=/work/mouse_dataset
OUTPUT_DIR=/work/mouse_map
# Define the genome index and whitelist directory
GENOME_INDEX=/work/mapping_files/Reference/ref_mouse 
WHITELIST_DIR=/work/mapping_files/whitelist/737K-august-2016.txt 
# Define an associative array to match R1 and R2 files for each sample
declare -A SAMPLE_PAIRS=(
    ["69"]="78"
    ["70"]="80"
    ["71"]="81"
    ["72"]="82"
    ["73"]="83"
    ["74"]="84"
    ["75"]="85"
    ["76"]="86"
)

# Loop over the samples
for SAMPLE_ID in "${!SAMPLE_PAIRS[@]}"; do
    # Get the corresponding R2 sample ID
    SAMPLE_ID_R2=${SAMPLE_PAIRS[$SAMPLE_ID]}

    # Define the input file names
    INPUT_DIR_R1=${INPUT_DIR}/SRR164098${SAMPLE_ID}_1.fastq.gz
    INPUT_DIR_R2=${INPUT_DIR}/SRR164098${SAMPLE_ID_R2}_2.fastq.gz

    # Define the output file names
    OUTPUT_PREFIX=${OUTPUT_DIR}/${SAMPLE_ID}
    OUTPUT_FILE=${OUTPUT_PREFIX}.bam

    # Run the starsolo command in the background
    STAR --genomeDir $GENOME_INDEX \
         --readFilesIn $INPUT_DIR_R2 $INPUT_DIR_R1 \
         --soloType CB_UMI_Simple \
         --soloCBwhitelist $WHITELIST_DIR \
         --soloFeatures Gene GeneFull SJ Velocyto \
         --soloOutFileNames $OUTPUT_PREFIX \
         --outSAMtype BAM SortedByCoordinate \
         --readFilesCommand gunzip -c \
         --outTmpDir $OUTPUT_PREFIX \
         &

done

# Wait for all the background jobs to finish
wait



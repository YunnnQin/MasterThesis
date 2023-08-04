## Human dataset

# Download raw fastq using sra-toolkit
# Three samples: GSM5687879	HSB231_EC, GSM5687887	HSB237_EC, GSM5687889	HSB628_EC
# Prefetch SRR16928876 SRR16928884 SRR16928886
./sratoolkit.3.0.1-centos_linux64/bin/prefetch SRR16928876
./sratoolkit.3.0.1-centos_linux64/bin/prefetch SRR16928884
./sratoolkit.3.0.1-centos_linux64/bin/prefetch SRR16928886

# Download fastq.gz file using fasterq-dump
./sratoolkit.3.0.1-centos_linux64/bin/fastq-dump --gzip --split-files ./SRR16928876/SRR16928876.sra
./sratoolkit.3.0.1-centos_linux64/bin/fastq-dump --gzip --split-files ./SRR16928884/SRR16928884.sra
./sratoolkit.3.0.1-centos_linux64/bin/fastq-dump --gzip --split-files ./SRR16928886/SRR16928886.sra


## Mouse dataset
#(single-end data)

# Download raw fastq using sra-toolkit simulteously
# Sample name: GSM5629170	Mouse 10x ENT
# Run: 27 runs  SRR16409868 to SRR16409894 (9 replicates)
seq 80 87 | xargs -P 8 -I {} sh -c "./sratoolkit.3.0.1-centos_linux64/bin/fasterq-dump --split-files ./SRR164098{}/SRR164098{}.sra "
./sratoolkit.3.0.1-centos_linux64/bin/fastq-dump --gzip --split-files ./SRR16409878/SRR16409878.sra
seq 69 77 | xargs -P 9 -I {} sh -c "./sratoolkit.3.0.1-centos_linux64/bin/fasterq-dump --split-files ./SRR164098{}/SRR164098{}.sra "
# R1 and R2 pairs
#SRR16409869 -- SRR16409878
#SRR16409870 -- SRR16409880
#SRR16409871 -- SRR16409881
#SRR16409872 -- SRR16409882
#SRR16409873 -- SRR16409883
#SRR16409874 -- SRR16409884
#SRR16409875 -- SRR16409885
#SRR16409876 -- SRR16409886
#SRR16409877 -- SRR16409887


## pig dataset

# Adult sample
# Runs: SRR9705091 - SRR9705094
seq 91 94 | xargs -P 4 -I {} sh -c "./sratoolkit.3.0.1-centos_linux64/bin/fasterq-dump --split-files ./SRR97050{}/SRR97050{}.sra "

# E70 sample
# Runs: SRR9705123 - SRR9705126
seq 23 26 | xargs -P 4 -I {} sh -c "./sratoolkit.3.0.1-centos_linux64/bin/fasterq-dump --split-files ./SRR97051{}/SRR97051{}.sra "
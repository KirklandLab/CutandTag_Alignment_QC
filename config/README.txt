Run hg38 samples in snakemake pipeline
- config.csv 
    + change bowtie2 genome index
    + change effective genome size
- samples.csv
    + change samples.csv to ones provided below


sample,fastq1,fastq2,sampleType,Control
H3.3_A,resources/trimmed_50_K27ac_S11_R1_001.fastq.gz,resources/trimmed_50_K27ac_S11_R2_001.fastq.gz,control,NA
H3.3_B,resources/trimmed_50_K27me3_S9_R1_001.fastq.gz,resources/trimmed_50_K27me3_S9_R2_001.fastq.gz,control,NA

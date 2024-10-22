Run hg38 samples in snakemake pipeline
- config.csv 
    + change bowtie2 genome index
    + change effective genome size
    + change macs2 genome size
- samples.csv
    + change samples.csv to ones provided below


sample,fastq1,fastq2,sampleType,Control
trimmed_50_K27ac,resources/trimmed_50_K27ac_S11_R1_001.fastq.gz,resources/trimmed_50_K27ac_S11_R2_001.fastq.gz,control,NA
trimmed_50_K27me3,resources/trimmed_50_K27me3_S9_R1_001.fastq.gz,resources/trimmed_50_K27me3_S9_R2_001.fastq.gz,control,NA

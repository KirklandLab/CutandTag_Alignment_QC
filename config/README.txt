Run hg38 samples in snakemake pipeline
- config.csv 
    + change bowtie2 genome index
    + change bamCoverage effective genome size
    + change macs2 genome size
- samples.csv
    + change samples.csv to ones provided driectly below



sample,fastq1,fastq2,sampleType,Control
K27ac_50_trimmed,resources/trimmed_50_K27ac_S11_R1_001.fastq.gz,resources/trimmed_50_K27ac_S11_R2_001.fastq.gz,control,NA
K27me3_50_trimmed,resources/trimmed_50_K27me3_S9_R1_001.fastq.gz,resources/trimmed_50_K27me3_S9_R2_001.fastq.gz,control,NA




Sample naming recommendation for correct plot output
- "Histone" + "_" + "Replicate" + "Any other identifier"
- Examples:
    + K27ac_50
    + K27me3_5
    + H3K27me3_rep1
    + H3K4me3_rep2_set1
    + H3K27ac_rep3_control
    + H3K27ac_rep3_treatment

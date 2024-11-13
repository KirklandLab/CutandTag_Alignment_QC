# CutandTag_Analysis_Snakemake
![Cut&Tag](/images/Cut&Tag.png)
+ OpenAI. (2024). Cartoon image of scissors cutting DNA and protein and DNA playing tag. DALL-E. Retrieved from OpenAI.

# 1) Project Description
CutAndTag_Analysis_Pipeline is a Snakemake pipeline adapted from the protocol by Ye Zheng, Kami Ahmad, and Steven Henikoff (dx.doi.org/10.17504/protocols.io.bjk2kkye). This pipeline is designed to process Cut-and-Tag sequencing data to facilitating the analysis of chromatin accessibility and DNA-protein interactions. It includes steps for quality control, read alignment, BAM to BigWig conversion, peak calling, and various visualizations, enabling high-resolution profiling of protein-DNA binding sites.

The pipeline provides automated quality checks, including FastQC and FastQ Screen reports, and performs read alignment with Bowtie2. Outputs include raw and scaled BigWig files, which allow users to visualize signal intensity across the genome. Additionally, peak calling is performed using MACS2, and fragment length and alignment summary plots are generated for in-depth data exploration. The workflow is automated with Snakemake, and dependencies are managed through module environments, ensuring reproducibility and flexibility.

A compact dataset is included within the repository for testing purposes, along with example scripts for analyzing publicly available Cut-and-Tag datasets. This pipeline extends the original protocol, offering a robust framework for both routine analysis and more complex studies.

# 2) Instructions to run on Slurm managed HPC
2A. Clone repository
```
git clone https://github.com/JK-Cobre-Help/CutandTag_Analysis_Snakemake.git
```
2B. Load modules
```
module purge
module load slurm python/3.10 pandas/2.2.3 numpy/1.22.3 matplotlib/3.7.1
```
2C. Modify Samples file
```
vim samples.csv
```
2D. Dry Run
```
snakemake -npr
```
2E. Run on HPC with config.yml options
```
sbatch --wrap="snakemake -j 999 --use-envmodules --latency-wait 60 --cluster-config config/cluster_config.yml --cluster 'sbatch -A {cluster.account} -p {cluster.partition} --cpus-per-task {cluster.cpus-per-task}  -t {cluster.time} --mem {cluster.mem} --output {cluster.output}'"
```

# 3) Explanation of samples.csv
Note. Make sure to check sample.csv before each run

The samples.csv file in the config folder has paths to the test fastq files. You must replace those paths with those for your own fastq files. The first column of each row is the sample name. This name will be used for all output files. Columns 2 and 3 are the paths to the paired fastq files. Column 4 is the sample type (either "treatment" or "control"). Column 5 is the name of the corresponding Control sample for each treated sample (use "NA" if the sample is a control).

| sample             | fastq1                        | fastq2                        | sampleType | Control   |
|--------------------|-------------------------------|-------------------------------|------------|-----------|
| K27ac_50_trimmed   | K27ac_50_trimmed_R1.fastq.gz  | K27ac_50_trimmed_R2.fastq.gz  | control    | NA        |
| K27me3_50_trimmed  | K27me3_50_trimmed_R1.fastq.gz | K27me3_50_trimmed_R1.fastq.gz | control    | NA        |


Sample naming recommendation for correct plot output
- "Histone" + "_" + "Replicate" + "Any other identifier"
- Examples:
    + K27ac_50
    + K27me3_5
    + K27ac_50_trimmed
    + H3K27me3_rep1
    + H3K4me3_rep2_set1
    + H3K27ac_rep3_control
    + H3K27ac_rep3_treatment

# 4) Citations
Zheng, Y., Ahmad, K., & Henikoff, S. (2019). CUT&Tag for efficient epigenomic profiling of small samples and single cells. Protocols.io, dx.doi.org/10.17504/protocols.io.bjk2kkye

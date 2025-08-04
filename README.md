![GitHub Release](https://img.shields.io/github/v/release/KirklandLab/CutandTag_Analysis_Snakemake)
![GitHub Release Date](https://img.shields.io/github/release-date/KirklandLab/CutandTag_Analysis_Snakemake)
![GitHub repo size](https://img.shields.io/github/repo-size/KirklandLab/CutandTag_Analysis_Snakemake)
![GitHub last commit](https://img.shields.io/github/last-commit/KirklandLab/CutandTag_Analysis_Snakemake)
![GitHub contributors](https://img.shields.io/github/contributors/KirklandLab/CutandTag_Analysis_Snakemake)
![GitHub Downloads (all assets, all releases)](https://img.shields.io/github/downloads/KirklandLab/CutandTag_Analysis_Snakemake/total)
![GitHub commits since latest release](https://img.shields.io/github/commits-since/KirklandLab/CutandTag_Analysis_Snakemake/latest)
[![DOI](https://zenodo.org/badge/873121124.svg)](https://doi.org/10.5281/zenodo.15232228)
[![License: Apache 2.0](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

# CutandTag_Analysis_Snakemake
![Cut&Tag](/images/Cut&Tag.png)
+ OpenAI. (2024). Cartoon image of scissors cutting DNA and protein and DNA playing tag. DALL-E. Retrieved from OpenAI.

# 1) Project Description
CutAndTag_Analysis_Pipeline is a Snakemake pipeline adapted from the protocol by Ye Zheng, Kami Ahmad, and Steven Henikoff (dx.doi.org/10.17504/protocols.io.bjk2kkye). This pipeline is designed to process Cut-and-Tag sequencing data to facilitating the analysis of chromatin accessibility and DNA-protein interactions. It includes steps for quality control, read alignment, BAM to BigWig conversion, peak calling, and various visualizations, enabling high-resolution profiling of protein-DNA binding sites.

The pipeline provides automated quality checks, including FastQC and FastQ Screen reports, and performs read alignment with Bowtie2. Outputs include raw and scaled BigWig files, which allow users to visualize signal intensity across the genome. Additionally, peak calling is performed using MACS2, and fragment length and alignment summary plots are generated for in-depth data exploration. The workflow is automated with Snakemake, and dependencies are managed through module environments, ensuring reproducibility and flexibility.

A compact dataset is included within the repository for testing purposes, along with example scripts for analyzing publicly available Cut-and-Tag datasets. This pipeline extends the original protocol, offering a robust framework for both routine analysis and more complex studies.

Downstream analysis can be performed in the [CutandTag_ReplicatePeak_Analysis](https://github.com/KirklandLab/CutandTag_ReplicatePeak_Analysis) snakemake pipeline. This pipeline starts with already aligned and filtered BAM files, focusing on the identification of reproducible peaks, the generation of consensus peak sets, and the visualization of overlaps and signal distributions across multiple samples or experimental conditions.

# 2) Explanation of `samples.csv`

**IMPORTANT:** Always update `samples.csv` with your sample, FASTQ file paths, histone, and replicate information before running.

This file should be placed in the `config/` directory and must contain the following columns:

| sample               | fastq1                                   | fastq2                                   | histone            | replicate    |
|----------------------|------------------------------------------|------------------------------------------|--------------------|--------------|
| H3K27ac_1_control    | path/to/H3K27ac_1_control_R1.fastq.gz    | path/to/H3K27ac_1_control_R2.fastq.gz    |H3H27ac_Control     |1             |
| H3K27ac_2_control    | path/to/H3K27ac_2_contorl_R1.fastq.gz    | path/to/H3K27ac_2_control_R1.fastq.gz    |H3K27ac_Control     |2             |
| H3K27me3_1_treatment | path/to/H3K27me3_1_treatment_R1.fastq.gz | path/to/H3K27me3_1_treatment_R1.fastq.gz |H3K27me3_Treatment  |1             |
| H3K27me3_2_treatment | path/to/H3K27me3_2_treatment_R1.fastq.gz | path/to/H3K27me3_2_treatment_R1.fastq.gz |H3K27me3_Treatment  |2             |

+ **sample**: Unique name used for all output files.
+ **fastq1 / fastq2**: Paths to R1 and R2 FASTQ files.
+ **histone**: Mark name (used in plots). Include control or treatment label in `histone`_`control` format
+ **replicate**: Include replicate number (e.g., `1`, `2`, `A`, `B`, etc.)

Use informative `sample` names that match your design. Examples:

+ `H3K27ac_rep1_control`
+ `H3K4me3_rep2_treatment`
+ `H3.3_A_trimmed`
+ `H2AX_1_H3K27ac`

Include both the histone and samlpe type in `histone` column. This determines how samples are grouped together. Following these practices will improve downstream plot labeling and reproducibility metrics.

# 3) Explanation of config.yml
Note. Make sure to check config.yml for the appropriate genome alignment

The config.yml is used to identify the file path of the bowtie2 genome index, specify effective genome size and genome for macs2. There is also information about specific modules and version numbers to maintain dependencies in the snakemake workflow. Running the mm10 genome does not require any modifications to the config.yml. When using the hg38 genome the following need to be modified with the information provided in the config.yml but commented out.

Run hg38 samples in snakemake pipeline
- config.yml 
    + change bowtie2 genome index file path
    + change bamCoverage effective genome size
    + change macs2 genome size


# 4) Examples of Output from Test Files in Repo

Below are some example plots generated by the pipeline.  

| 1. **Alignment Summary Plot**                                                         | 2. **Peak Summary Plot**                                                                   |
| :-----------------------------------------------------------------------------------: | :----------------------------------------------------------------------------------------: |
| <img src="/images/alignment_summary_plot.png" width="300">                            | <img src="/images/peak_summary_plot.png" width="300">                                      |
| *Summary of alignment stats*                                                          | *Summary of called peak stats*                                                             |

| 3. **Fragment Length Plot**                                                           | 4. **Fragment Count Correlation Plot**                                                     |
| :-----------------------------------------------------------------------------------: | :----------------------------------------------------------------------------------------: |
| <img src="/images/fragment_length_plot.png" width="300">                              | <img src="/images/fragCount_correlation_plot.png" width="300">                             |
| *Fragment Length summary plots*                                                       | *Correlation plot of fragment counts*                                                      |


# 5) Instructions to run on Slurm managed HPC
5A. Download verson controlled repository
```
wget https://github.com/KirklandLab/CutandTag_Analysis_Snakemake/releases/download/v1.0.4/CutandTag_Analysis_Snakemake-1.0.4.tar.gz
tar -xzf CutandTag_Analysis_Snakemake-1.0.4.tar.gz
rm CutandTag_Analysis_Snakemake-1.0.4.tar.gz
cd CutandTag_Analysis_Snakemake-1.0.4
```
5B. Load modules
```
module purge
module load slurm python/3.10 pandas/2.2.3 numpy/1.22.3 matplotlib/3.7.1
```
5C. Modify samples and config file
```
vim samples.csv
vim config.yml
```
5D. Dry Run
```
snakemake -npr
```
5E. Run on HPC with config.yml options
```
sbatch --wrap="snakemake -j 999 --resources mem_mb=200000 --use-envmodules --latency-wait 300 --cluster-config config/cluster_config.yml --cluster 'sbatch -A {cluster.account} -p {cluster.partition} --cpus-per-task {cluster.cpus-per-task}  -t {cluster.time} --mem {cluster.mem} --output {cluster.output}'"
```

# 6) Citations
Zheng, Y., Ahmad, K., & Henikoff, S. (2019). CUT&Tag for efficient epigenomic profiling of small samples and single cells. Protocols.io, dx.doi.org/10.17504/protocols.io.bjk2kkye

![GitHub Release](https://img.shields.io/github/v/release/KirklandLab/CutandTag_Alignment_QC)
![GitHub Release Date](https://img.shields.io/github/release-date/KirklandLab/CutandTag_Alignment_QC)
![GitHub repo size](https://img.shields.io/github/repo-size/KirklandLab/CutandTag_Alignment_QC)
![GitHub last commit](https://img.shields.io/github/last-commit/KirklandLab/CutandTag_Alignment_QC)
![GitHub contributors](https://img.shields.io/github/contributors/KirklandLab/CutandTag_Alignment_QC)
![GitHub Downloads (all assets, all releases)](https://img.shields.io/github/downloads/KirklandLab/CutandTag_Alignment_QC/total)
![GitHub commits since latest release](https://img.shields.io/github/commits-since/KirklandLab/CutandTag_Alignment_QC/latest)
[![DOI](https://zenodo.org/badge/873121124.svg)](https://doi.org/10.5281/zenodo.15232228)
[![License: Apache 2.0](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

# CutandTag_Alignment_QC
![Cut&Tag](/images/Cut&Tag.png)
*(Image generated with DALL-E. OpenAI, 2024: Cartoon of scissors cutting DNA and DNA tagging protein)*

## 1) Project Description  

**CutandTag_Alignment_QC** is a Snakemake workflow adapted from the protocol by Ye Zheng, Kami Ahmad, and Steven Henikoff ([dx.doi.org/10.17504/protocols.io.bjk2kkye](https://dx.doi.org/10.17504/protocols.io.bjk2kkye)). This pipeline is designed to process Cut-and-Tag sequencing data to facilitate the analysis of chromatin accessibility and DNA-protein interactions. It uses raw FASTQs as input and includes automated steps for read alignment, quality control, signal track generation (raw and scaled), peak calling, and various visualizations, enabling high-resolution profiling of protein-DNA binding sites. The pipeline performs automated quality checks using FastQC, MultiQC, and FastQ Screen reports while also aligning reads with Bowtie2. Outputs include raw and scaled BigWig files, which allow users to visualize signal intensity across the genome. Additionally, peak calling is performed using MACS2, and fragment length and alignment summary plots are generated for detailed data exploration. The workflow is automated with Snakemake and dependencies are managed through module environments, ensuring reproducibility and flexibility.  

Downstream analysis can be performed in the [CutandTag_ReplicatePeak_Analysis](https://github.com/KirklandLab/CutandTag_ReplicatePeak_Analysis) snakemake workflow. This pipeline starts with aligned BAM files, focusing on the identification of reproducible peaks, the generation of consensus peak sets, and the visualization of overlaps and signal distributions across multiple samples and experimental conditions.  

### **Key Features**

+ **Raw QC & Contamination Check**
  + Generates **FastQC** reports for raw FASTQ files
  + Aggregates reports into a single **MultiQC** summary
  + Performs contamination detection via **FastQ Screen**

+ **Read Alignment & Coverage**
  + Aligns reads using **Bowtie2**
  + Converts SAM to sorted BAM files and generates BAM indices
  + Produces **raw**, **CPM-normalized**, and **custom-scaled** BigWig tracks using DeepTools  
    + **CPM** scaling enables comparison across samples by normalizing to 1 million reads  
    + **Custom scale factor** normalization adjusts coverage more precisely based on total mapped reads, ideal for datasets with substantial read depth variation

+ **Peak Calling**
  + Performs narrow peak calling with **MACS2**
  + Uses a customizable q-value threshold and genome size

+ **Summary Plots**
  + Alignment stats (total reads, alignment rate)
  + Fragment length distribution
  + FRiP scores and peak count summaries
  + Fragment count correlation between replicates

+ **Modular & Reproducible**
  + Designed for HPC environments using module-based environments
  + Uses `samples.csv` and `config.yml` for full customization
  + Easily integrated with the downstream pipeline: [CutandTag_ReplicatePeak_Analysis](https://github.com/KirklandLab/CutandTag_ReplicatePeak_Analysis)

---

## 2) Intended Use Case

This pipeline is ideal for processing Cut-and-Tag datasets:

+ Performs **initial alignment and QC** of Cut-and-Tag sequencing data  
+ Creates BigWig tracks for genome browser visualization  
+ Summarizes sample-level metrics  

**Note**: This pipeline does **not** merge replicates or define consensus peaks. That is handled by the companion pipeline linked above.

---

## 3) Dependencies and Configuration

All parameters and module versions are specified in `config/config.yml`.

### **Key fields include:**

+ `bowtie2_genome`: path to the Bowtie2 index for the reference genome (e.g., mm10 or hg38)
+ `effective_genome_size`: used by DeepTools to calculate coverage normalization
+ `genome_size`: effective genome size string used by MACS2 (e.g., `mm`, `hs`)
+ `binSize`: bin size for BigWig coverage generation
+ `macs2_qvalue`: q-value threshold for MACS2 peak calling
+ `fastqc`, `fastq_screen`, `multiqc`, `bowtie2`, `samtools`, `deeptools`, `bedtools`, `macs2`, `R`, `bioconductor`: module names and versions for use on an HPC

### **Changing Genomes**

+ To switch from mm10 to hg38:
  + Update `bowtie2_genome` path to the new index
  + Change `effective_genome_size` to the appropriate value (e.g., `2913022398` for hg38)
  + Set `genome_size` to `hs` for human in MACS2

### **Tool Versions and Resources**

+ The `config/config.yml` file defines module versions to load with `--use-envmodules`
+ The `config/cluster_config.yml` file defines resource usage per rule (e.g., memory, time, cores)

---

## 4) Tools & Modules

This workflow uses the following tools through environment modules on an HPC system:

+ **FastQC**: for assessing raw FASTQ file quality  
+ **FastQ Screen**: for detecting sample contamination via database alignment  
+ **MultiQC**: for aggregating all QC outputs into a single HTML report  
+ **Bowtie2**: for aligning paired-end reads to a reference genome  
+ **Samtools**: for converting and indexing SAM/BAM files  
+ **DeepTools**: `bamCoverage`: generates BigWig tracks (raw, CPM, or scaled)  
+ **MACS2**: for narrow peak calling using BAMPE mode  
+ **Bedtools**: for fragment extraction, cleanup, and binning  
+ **R**: for generating summary plots (alignment stats, fragment length, correlation, peak summary)  
+ **Bioconductor**: for R-based QC visualizations  
+ **Python**: for calculating scale factors and performing workflow logic

**Note**: All tool paths, versions, and parameters are centrally managed in `config.yml`, making the pipeline reproducible, portable, and easy to maintain.

---

## 5) Example Data  

A compact dataset is included within the repository for testing purposes (provided by Jacob Kirkland), along with example scripts for analyzing publicly available Cut-and-Tag datasets. This pipeline extends the original protocol, offering a robust framework for both routine analysis and more complex studies.  

---

## 6) Explanation of `samples.csv`

**IMPORTANT:** Always update `config/samples.csv` with your sample, FASTQ file paths, histone, and replicate information before running.

The `samples.csv` must contain the following columns:

| sample               | fastq1                                   | fastq2                                   | histone            | replicate    |
|----------------------|------------------------------------------|------------------------------------------|--------------------|--------------|
| H3K27ac_1_control    | path/to/H3K27ac_1_control_R1.fastq.gz    | path/to/H3K27ac_1_control_R2.fastq.gz    |H3K27ac_Control     |1             |
| H3K27ac_2_control    | path/to/H3K27ac_2_control_R1.fastq.gz    | path/to/H3K27ac_2_control_R1.fastq.gz    |H3K27ac_Control     |2             |
| H3K27me3_1_treatment | path/to/H3K27me3_1_treatment_R1.fastq.gz | path/to/H3K27me3_1_treatment_R1.fastq.gz |H3K27me3_Treatment  |1             |
| H3K27me3_2_treatment | path/to/H3K27me3_2_treatment_R1.fastq.gz | path/to/H3K27me3_2_treatment_R1.fastq.gz |H3K27me3_Treatment  |2             |

+ **sample**: Unique name used for all output files.
+ **fastq1 / fastq2**: Paths to R1 and R2 FASTQ files.
+ **histone**: Mark name (used in plots). Include control or treatment label in `histone`_`control` format
+ **replicate**: Include replicate number (e.g., `1`, `2`, `A`, `B`, etc.)

Include both the histone and sample type in `histone` column. This determines how samples are grouped together. Following these practices will improve downstream plot labeling and reproducibility metrics.  

Use informative `sample` names that match your design and be sure that fastq file names do not include "." before the file extension.  

Example sample names:  
+ `H3K27ac_rep1_control`  
+ `H3K4me3_rep2_treatment`  
+ `H3-3_A_trimmed`  
+ `H2AX_1_H3K27ac`  

Example FASTQ names:
+ `H3K27ac_1_control_R1.fastq.gz`
+ `JKTS192_S63_R1_001.fastq.gz`

Example FASTQ names to **not** use:
+ `H3.3_A_S21_R1_001.fastq.gz`
+ `Sample.Rep.A_S09_R1.fastq.gz`

---

## 7) Output Overview

| Category                     | Output Folder/Format                                               |
|-----------------------------|---------------------------------------------------------------------|
| **QC**                      | `results/qc/fastqc/`, `results/qc/multiqc/`, `results/qc/fastq_screen/` |
| **Alignments**              | `results/alignment/sam/`, `results/alignment/bam/` (sorted + .bai) |
| **BigWig Tracks**           | `results/alignment/bigwig/` (raw, CPM, scaled)                     |
| **Peak Calls**              | `results/peakCalling/` `.narrowPeak` files                         |
| **Fragment Length Files**   | `results/alignment/sam/*_fragmentLen.txt`, summary plot            |
| **FRiP Score & Coverage**   | `results/alignment/bam/*_frip.txt`, `*_fragmentsCount.bin500.bed`  |
| **Plots**                   | `results/plots/` for alignment summary, FRiP, and peak stats       |

**Note**: BigWig Tracks Explanation  
  + **CPM**: normalizes coverage to 1 million reads using `--normalizeUsing CPM`  
  + **Custom scale factor**: calculated as `1 / (total_mapped_reads / 1,000,000)` and applied with `--scaleFactor`, allowing more precise normalization across varying sequencing depths  

---

## 8) Example Output Plots

Below are example plots generated by this pipeline.  

| 1. **Alignment Summary Plot**                                                         | 2. **Peak Summary Plot**                                                                   |
| :-----------------------------------------------------------------------------------: | :----------------------------------------------------------------------------------------: |
| <img src="/images/alignment_summary_plot.png" width="300">                            | <img src="/images/peak_summary_plot.png" width="300">                                      |
| *Summary of alignment stats*                                                          | *Summary of called peak stats*                                                             |

| 3. **Fragment Length Plot**                                                           | 4. **Fragment Count Correlation Plot**                                                     |
| :-----------------------------------------------------------------------------------: | :----------------------------------------------------------------------------------------: |
| <img src="/images/fragment_length_plot.png" width="300">                              | <img src="/images/fragCount_correlation_plot.png" width="300">                             |
| *Fragment Length summary plots*                                                       | *Correlation plot of fragment counts*                                                      |

---

## 9) Instructions to run on Slurm managed HPC  
9A. Download version controlled repository
```
wget https://github.com/KirklandLab/CutandTag_Alignment_QC/releases/download/v1.0.4/CutandTag_Alignment_QC-1.0.4.tar.gz
tar -xzf CutandTag_Alignment_QC-1.0.4.tar.gz
rm CutandTag_Alignment_QC-1.0.4.tar.gz
cd CutandTag_Alignment_QC-1.0.4
```
9B. Load modules
```
module purge
module load slurm python/3.10 pandas/2.2.3 numpy/1.22.3 matplotlib/3.7.1
```
9C. Modify samples and config file
```
vim samples.csv
vim config.yml
```
9D. Dry Run
```
snakemake -npr
```
9E. Run on HPC with config.yml options
```
sbatch --wrap="snakemake -j 999 --resources mem_mb=200000 --use-envmodules --latency-wait 300 --cluster-config config/cluster_config.yml --cluster 'sbatch -A {cluster.account} -p {cluster.partition} --cpus-per-task {cluster.cpus-per-task}  -t {cluster.time} --mem {cluster.mem} --output {cluster.output}'"
```

---

## 10) Citation  

If you use this workflow in your research, please cite both the original CUT&Tag protocol and this pipeline:  

+ Zheng, Y., Ahmad, K., & Henikoff, S. (2019). CUT&Tag for efficient epigenomic profiling of small samples and single cells. *Protocols.io*. https://dx.doi.org/10.17504/protocols.io.bjk2kkye  

+ **Boyd, K.A.** (2025). *CutandTag_Alignment_QC: A reproducible Snakemake workflow for alignment, quality control, and signal generation from Cut-and-Tag sequencing data*. *Zenodo*. https://doi.org/10.5281/zenodo.15232228  

[![DOI](https://zenodo.org/badge/873121124.svg)](https://doi.org/10.5281/zenodo.15232228)  

---

## 11) Authorship & Contributions  

**Kevin A. Boyd** – Designed and implemented the Snakemake workflow for a Slurm-managed HPC environment, modularized the pipeline structure, implemented all processing steps, integrated spike-in normalization support, designed quality control plots, and created the documentation.   

**Jacob Kirkland** – Principal Investigator; provided project direction, conceptual guidance, and experimental data for pipeline development.  

This work was developed under the guidance of Jacob Kirkland as part of a COBRE-funded collaborative effort. While the pipeline was built specifically for use within the Kirkland Lab, it is broadly applicable to Cut-and-Tag data analysis in other research settings.  

---

## 12) License

This project is licensed under the **Apache 2.0**. See the [LICENSE](LICENSE) file for details.  

[![License: Apache 2.0](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)


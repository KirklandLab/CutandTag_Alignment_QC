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

--- 

## 1) Project Description

**CutandTag_Alignment_QC** is a Snakemake workflow adapted from the protocol by Ye Zheng, Kami Ahmad, and Steven Henikoff ([dx.doi.org/10.17504/protocols.io.bjk2kkye](https://dx.doi.org/10.17504/protocols.io.bjk2kkye)). This pipeline is designed to process CUT&Tag sequencing data from raw paired end FASTQ files through raw FASTQ quality control, read alignment, optional duplicate capping, optional random downsampling, BigWig signal track generation, per sample peak calling, fragment length analysis, and QC visualization.

The workflow aligns reads with **Bowtie2**, creates sorted and indexed BAM files, optionally limits PCR duplicate burden using fragment coordinate duplicate capping, optionally downscales samples to a defined fragment depth, and outputs a final analysis BAM for each sample. Final analysis BAMs are then used to generate raw, CPM normalized, and target scaled BigWig files for genome browser visualization, call per sample peaks with MACS2, calculate raw and final fragment length distributions, compute FRiP scores, generate fragment count correlations, and produce alignment, duplicate/downsampling, and peak summary plots.

Automated raw FASTQ quality checks are done using **FastQC**, **MultiQC**, and **FastQ Screen**. These raw FASTQ QC steps are controlled by a config toggle, making it possible to skip them when rerunning the workflow after an initial QC has already been completed. The workflow is automated with Snakemake and dependencies are managed through environment modules, supporting reproducible execution on a Slurm managed HPC system.

Downstream analysis can be performed in the [CutandTag_ReplicatePeak_Analysis](https://github.com/KirklandLab/CutandTag_ReplicatePeak_Analysis) Snakemake workflow. This companion pipeline starts with aligned BAM files and focuses on identifying reproducible peaks, generating consensus peak sets, and visualizing overlaps and signal distributions across multiple samples and experimental conditions.

### **Key Features**
*Note: Optional features can be toggled on and off in the `config` file found in `config/config.yml`.  

+ **Optional: Raw QC & Contamination Check**
  + Generates **FastQC** reports for raw FASTQ files
  + Aggregates reports into a single **MultiQC** summary
  + Performs contamination detection on R1 via **FastQ Screen**
  + Can be toggled on/off using `use_fastq_qc`

+ **Optional: Adapter Trimming**
  + Can be toggled on/off using `use_trimming`
  + Uses **Cutadapt** to remove the configured adapter sequence from paired-end reads
  + Produces trimmed paired-end FASTQ files and a per-sample trimming report
  + When trimming is enabled, Bowtie2 aligns the trimmed FASTQs using local alignment mode
  + When trimming is disabled, Bowtie2 aligns the original FASTQs using end-to-end alignment mode

+ **Read Alignment & BAM Processing**
  + Aligns paired end reads using **Bowtie2**
  + Converts SAM to BAM, then sorts and indexes BAM files
  + Produces sorted BAM files for each sample
  + Defines a final analysis BAM based on duplicate capping and downsampling settings
  + Includes preconfigured reference settings for the mouse **mm10**/**mm39** and human **hg38** genome builds
  + Can be used with other organisms or genome assemblies by providing a compatible Bowtie2 index and updating the corresponding genome-size settings in the `config/config.yml`

+ **Optional: Duplicate Capping**
  + Can be toggled on/off using `use_duplicate_cap`
  + Caps duplicate fragments by exact genomic fragment coordinates
  + Retains up to `duplicate_cap_max` identical fragments at each fragment position
  + Helps reduce the influence of PCR duplicates without forcing complete duplicate removal
  + Produces per sample duplicate metrics and duplicate/downsampling summary plots

+ **Optional: Random Downsampling**
  + Can be toggled on/off using `use_downsampling`
  + Supports three target selection modes:
    + `manual`
    + `lowest`
    + `lowest_with_floor`
  + All downsampling methods:
    + Supports either a fixed seed for reproducibility or `"random"` to generate a new seed when downsampling is performed
    + Downsample only samples above the resolved target
    + Leaves samples below the target unchanged
    + Produces per sample downsampling metrics and final analysis fragment counts

+ **Read Alignment & Coverage**
  + Aligns reads using **Bowtie2**
  + Converts SAM to sorted BAM files and generates BAM indices
  + Produces **raw**, **CPM normalized**, and **target scaled** BigWig tracks using DeepTools
    + **Raw** BigWigs show direct coverage from the final analysis BAM
    + **CPM** BigWigs normalize coverage using `bamCoverage --normalizeUsing CPM`
    + **Target scaled** BigWigs scale signal using a resolved fragment depth target or empirical scaling reference
      + *When downsampling is **enabled**, the scaling target is the resolved downsampling target*
      + *When downsampling is **disabled**, no BAM downsampling target is used; the target scaled BigWig reference is resolved to the empirical lowest final analysis fragment count*

+ **Peak Calling**
  + Performs controll free per sample peak calling with **MACS2**
  + Uses paired end BAM mode
  + Uses a customizable q-value threshold and genome size setting
  + Uses `--keep-dup all`, so duplicate handling is determined by the upstream duplicate-capping configuration rather than by MACS2
  + The workflow does not currently assign a matched IgG or other control BAM to each sample during MACS2 peak calling

+ **Summary Plots**
  + Alignment statistics
  + Duplicate burden and downsampling behavior
  + Final analysis fragment length distributions
  + FRiP scores and peak count summaries
  + Fragment count correlation between samples

+ **Modular & Reproducible**
  + Designed for HPC environments using module based software environments
  + Uses `samples.csv` and `config.yml` for full customization
  + Uses Snakemake to track dependencies and rule outputs
  + Easily integrated with the downstream pipeline: [CutandTag_ReplicatePeak_Analysis](https://github.com/KirklandLab/CutandTag_ReplicatePeak_Analysis)

+ **Extra Script: Heatmap Generation**
  + A helper script is included in `scripts/make_heatplot.sh` for visualizing BigWig signals over defined BED regions
  + Instructions for customizing and running the script are provided as comments in the script
 
+ **Extra Script: FASTQ Combining**
  + A helper script is included in `scripts/combine_fastqs.sh` for combining multiple FASTQ files from the same biological sample
  + This is useful when additional sequencing is returned by the sequencing core, such as reads from multiple lanes, sequencing runs, or sequencing submissions that should be analyzed as one sample
  + Instructions for customizing and running the script are provided as comments in the script

---

## 2) Intended Use Case

This pipeline is designed for the initial processing and quality assessment of paired end CUT&Tag sequencing data. The main goal is to produce clean, well documented alignment files, browser tracks, and sample level metrics that can be used to evaluate data quality before downstream replicate aware analysis.

+ Performs **initial alignment and QC** of CUT&Tag sequencing data
+ Optionally runs raw FASTQ QC using **FastQC**, **MultiQC**, and **FastQ Screen**
+ Optionally trims paired-end adapter sequence using **Cutadapt**
+ Creates sorted BAM files and config dependent final analysis BAMs
+ Optionally caps duplicate fragments to reduce PCR duplicate burden
+ Optionally downsamples samples to a manual, lowest sample, or floor filtered target depth
+ Generates raw, CPM normalized, and target scaled BigWig tracks for genome browser visualization
+ Calls per sample MACS2 peaks for QC and sample level assessment
+ Calculates raw and final fragment length distributions
+ Calculates FRiP scores using per sample peak calls
+ Summarizes sample level alignment, duplication, downsampling, fragment length, FRiP, peak count, and correlation metrics
+ Provides final analysis BAMs and signal tracks for downstream workflows

**Note**: This pipeline is intended for alignment, preprocessing, signal generation, and sample level QC. It does **not** merge replicates, define consensus peaks, or perform replicate aware differential analysis. Those steps are handled by the companion pipeline linked above.

---

## 3) Dependencies and Configuration

All parameters and module versions are specified in `config/config.yml`.

### **Key fields include:**

+ `use_trimming`: whether to trim paired-end FASTQ files before alignment
+ `cutadapt_adapter_r1`: adapter sequence removed from R1
+ `cutadapt_adapter_r2`: adapter sequence removed from R2
+ `cutadapt_minimum_length`: minimum read length retained after trimming
+ `cutadapt_error_rate`: maximum allowed adapter-matching error rate
+ `cutadapt_minimum_overlap`: minimum adapter overlap required for trimming
+ `bowtie2_genome`: path to the Bowtie2 index for the reference genome, such as mm10 or hg38
+ `effective_genome_size`: effective genome size used by DeepTools for coverage normalization
+ `genome_size`: genome size string used by MACS2, such as `mm` or `hs`
+ `binSize`: bin size for BigWig coverage generation
+ `macs2_qvalue`: q-value threshold for MACS2 peak calling
+ `use_fastq_qc`: whether to run FastQC, MultiQC, and FastQ Screen
+ `use_duplicate_cap`: whether to cap duplicate fragments
+ `duplicate_cap_max`: maximum number of identical fragments retained at each genomic position
+ `use_downsampling`: whether to randomly downsample the duplicate-capped or aligned BAM when creating final analysis BAMs
+ `downsample_target_mode`: method used to select the downsampling target
+ `downsample_target_fragments`: manual target used when `downsample_target_mode: "manual"`
+ `downsample_minimum_acceptable_fragments`: minimum acceptable depth used when `downsample_target_mode: "lowest_with_floor"`
+ `downsample_seed`: random seed used for reproducible downsampling
+ `fastqc`, `fastq_screen`, `multiqc`, `bowtie2`, `samtools`, `deeptools`, `bedtools`, `macs2`, `R`, `bioconductor`: module names and versions for use on an HPC

---

### **Adapter Trimming Settings**

Adapter trimming can be toggled on or off:

```yaml
use_trimming: true
```

When `use_trimming: true`, paired-end reads are processed with Cutadapt using the configured R1 and R2 adapter sequences. The trimmed FASTQ files are used for Bowtie2 alignment, and Bowtie2 runs in `--local` alignment mode.

```yaml
use_trimming: false
```

When `use_trimming: false`, the original FASTQ files are aligned directly, and Bowtie2 runs in `--end-to-end` alignment mode.

The trimming parameters are controlled by:

```yaml
cutadapt_adapter_r1: "CTGTCTCTTATACACATCT"
cutadapt_adapter_r2: "CTGTCTCTTATACACATCT"
cutadapt_minimum_length: 1
cutadapt_error_rate: 0.10
cutadapt_minimum_overlap: 5
```

Cutadapt reports are written to `results/qc/cutadapt/`, and trimmed FASTQ files are written to `results/trimming/`.

---

### **Raw FASTQ QC Toggle**

Raw FASTQ QC can be toggled on or off:

```yaml
use_fastq_qc: true
```

When `use_fastq_qc: true`, the workflow runs:
+ FastQC
+ MultiQC
+ FastQ Screen

When `use_fastq_qc: false`, these rules are skipped unless their output files are requested directly.

```yaml
use_fastq_qc: false
```

This is useful when rerunning the workflow multiple times using the same FASTQs, because raw FASTQ QC usually does not need to be repeated after it has already been initially checked.

---

### **Duplicate Capping Settings**

Duplicate capping can be toggled on or off:

```yaml
use_duplicate_cap: true
duplicate_cap_max: 5
```

When enabled, the workflow groups paired end fragments by exact genomic fragment coordinates and retains up to `duplicate_cap_max` identical fragments per coordinate.

For example, with:

```yaml
duplicate_cap_max: 5
```

The workflow allows up to 5 identical fragments at the same position. Additional identical fragments are removed before final analysis BAM generation.

To disable duplicate capping:

```yaml
use_duplicate_cap: false
duplicate_cap_max: 5
```

*Note: When duplicate capping is disabled, `duplicate_cap_max` is ignored.*

**Important:** This step is best described as **duplicate capping**, not complete duplicate removal. A duplicate cap of `5` does not remove all duplicate fragments. It keeps up to 5 identical fragments per fragment position and removes only duplicate fragments beyond that cap.

---

### **Random Downsampling Settings**

Downsampling can be toggled on or off:

```yaml
use_downsampling: true
```

When enabled, the workflow resolves a target fragment depth and randomly subsamples samples above that target toward the resolved depth. Samples below the target are not upsampled and are left unchanged.

When disabled:

```yaml
use_downsampling: false
```

the workflow does not downsample any BAM files. However, the workflow still records the empirical lowest final analysis fragment count as a BigWig scaling reference for `analysisScaled.bw`.

The downsampling target is selected using:

```yaml
downsample_target_mode: "manual"
```

Valid options are:

+ `manual`
+ `lowest`
+ `lowest_with_floor`

Downsampling uses `samtools view --subsample`. This method performs probabilistic subsampling so the resulting fragment count should be close to, but may not exactly equal, the resolved target. The actual final fragment count is recorded in the per-sample downsampling metrics. For reproducible testing or direct run to run comparisons, set a fixed seed:

```yaml
downsample_seed: 12345
```

*Note: This makes random downsampling reproducible when the same input BAM, target, and seed are used. A fixed seed is recommended when comparing duplicate cap or downsampling settings across repeated test runs.*


For routine operation, the workflow can instead generate a new seed:

```yaml
downsample_seed: "random"
```

*Note: When `downsample_seed: "random"` is used, the workflow generates and records a new seed for each sample that is actually downsampled. Repeated runs may therefore select different fragments.*

---

### **Downsampling: Manual Target**

**manual** mode uses the target specified by `downsample_target_fragments`.

```yaml
use_downsampling: true
downsample_target_mode: "manual"
downsample_target_fragments: 15000000
downsample_seed: "random"
```

The resolved target is set by the config:
```text
15,000,000
```

Example post duplicate cap fragment depths:
```text
Sample 1 = 30,000,000
Sample 2 = 28,000,000
Sample 3 = 18,000,000
Sample 4 = 24,000,000
Sample 5 =  8,000,000
```

Example downsampled final analysis depths are:
```text
Sample 1: 30,000,000 -> 15,000,000
Sample 2: 28,000,000 -> 15,000,000
Sample 3: 18,000,000 -> 15,000,000
Sample 4: 24,000,000 -> 15,000,000
Sample 5:  8,000,000 ->  8,000,000
```

Sample 5 is below the target, so it is not changed.

*Note: This mode is useful when a specific target depth is desired across experiments or when comparing several runs using a consistent manually selected threshold.*

---

### **Downsampling: Lowest Sample With Floor**

**lowest_with_floor** mode ignores samples below a minimum acceptable depth when choosing the downsampling target. It uses the lowest depth sample that is at or above the specified floor in `downsample_minimum_acceptable_fragments`.

```yaml
use_downsampling: true
downsample_target_mode: "lowest_with_floor"
downsample_minimum_acceptable_fragments: 15000000
downsample_seed: "random"
```

Example post duplicate cap fragment depths:
```text
Sample 1 = 30,000,000
Sample 2 = 28,000,000
Sample 3 = 18,000,000
Sample 4 = 24,000,000
Sample 5 =  8,000,000
```

Sample 5 is below the 15,000,000 floor and is not used to choose the target. The lowest sample at or above the floor is Sample 3:
```text
Resolved target = 18,000,000
```

Example downsampled final analysis depths are:
```text
Sample 1: 30,000,000 -> 18,000,000
Sample 2: 28,000,000 -> 18,000,000
Sample 3: 18,000,000 -> 18,000,000
Sample 4: 24,000,000 -> 18,000,000
Sample 5:  8,000,000 ->  8,000,000
```

Sample 5 is below the target, so it is not changed. Sample 3 set the target, so it is also unchanged. 

*Note: This mode prevents one very low depth sample from forcing all other samples down to an unnecessarily low depth when at least one sample meets the floor. If no sample is at or above the floor, the workflow falls back to the overall lowest available fragment depth.*

---

### **Downsampling: Lowest Sample**

**lowest** mode uses the lowest available post cap or aligned fragment depth as the target, as `downsample_target_fragments` and `downsample_minimum_acceptable_fragments` are not used to define the target.

```yaml
use_downsampling: true
downsample_target_mode: "lowest"
downsample_seed: "random"
```

Example post duplicate cap fragment depths:
```text
Sample 1 = 30,000,000
Sample 2 = 28,000,000
Sample 3 = 18,000,000
Sample 4 = 24,000,000
Sample 5 =  8,000,000
```

The resolved target is set by the low sample:
```text
8,000,000
```

Example downsampled final analysis depths are:
```text
Sample 1: 30,000,000 -> 8,000,000
Sample 2: 28,000,000 -> 8,000,000
Sample 3: 18,000,000 -> 8,000,000
Sample 4: 24,000,000 -> 8,000,000
Sample 5:  8,000,000 -> 8,000,000
```

Sample 5 is the lowest sample so it set the target for all other samples. 

*Note: This mode is strict and can force high depth samples down to the depth of the lowest depth sample. It can be useful when the lowest sample is still high quality, but it may be too aggressive if one sample has unusually low depth.*

---

### **How Final Analysis BAMs Are Defined**

The final analysis BAM depends on the duplicate capping and downsampling settings.

+ If duplicate capping is enabled and downsampling is enabled:
  + Final BAM = duplicate capped BAM after downsampling

+ If duplicate capping is enabled and downsampling is disabled:
  + Final BAM = duplicate capped BAM

+ If duplicate capping is disabled and downsampling is enabled:
  + Final BAM = aligned sorted BAM after downsampling

+ If duplicate capping is disabled and downsampling is disabled:
  + Final BAM = aligned sorted BAM

These final analysis BAMs are used for:

+ BigWig generation
+ MACS2 peak calling
+ Final fragment length analysis
+ FRiP score calculation
+ Fragment count correlation analysis

The workflow does **not** create artificial reads or upsample BAM files. Samples below a downsampling target are left unchanged.

---

### **Changing Genomes**

The workflow includes preconfigured settings for the mouse mm10 and mm39 genome builds, along with the human hg38 build in the `config/config.yml`.

+ To switch from mm10 to hg38:
  + Update `bowtie2_genome` path to the new hg38 Bowtie2 index
  + Change `effective_genome_size` to the appropriate value, such as `2913022398` for hg38
  + Set `genome_size` to `hs` for human in MACS2 peak calling
 
+ To switch from mm10 to mm39:
  + Update `bowtie2_genome` path to the new mm39 Bowtie2 index
  + Change `effective_genome_size` to the appropriate value, such as `2654621783` for mm39
  + Keep `genome_size` set to `mm` for mouse in MACS2 peak calling
 
+ Other organisms or genome assemblies can also be used:
  + Build or obtain a compatible Bowtie2 index
  + Update `bowtie2_genome` to the location of the new index
  + Update `effective_genome_size` for DeepTools normalization
  + Set the MACS2 `genome_size` option to one of the built-in codes:
    + `hs` for human
    + `mm` for mouse
    + `ce` for *C. elegans*
    + `dm` for *D. melanogaster*
  + For other organisms, provide an appropriate numeric effective genome size
  + Ensure that all downstream reference files, chromosome names, and BED files use the same genome assembly

### **Tool Versions and Resources**

+ The `config/config.yml` file defines module versions to load with `--use-envmodules`
+ The `config/cluster_config.yml` file defines resource usage per rule, including memory, time, and cores
+ The workflow is designed for Slurm managed HPC execution using Snakemake cluster submission

---

## 4) Tools & Modules

This workflow uses the following tools through environment modules on an HPC system:

+ **FastQC**: for assessing raw FASTQ file quality
+ **FastQ Screen**: for detecting sample contamination via database alignment
+ **MultiQC**: for aggregating QC outputs into a single HTML report
+ **Cutadapt**: for adapter trimming
+ **Bowtie2**: for aligning paired end reads to a reference genome
+ **Samtools**: for converting, sorting, indexing, filtering, downsampling, and counting SAM/BAM files
+ **DeepTools**: `bamCoverage` generates BigWig tracks from final analysis BAMs
+ **MACS2**: for narrow peak calling using paired end BAM mode
+ **Bedtools**: for fragment extraction, cleanup, and binning
+ **R**: for generating summary plots, including alignment stats, fragment length, duplicate/downsampling summaries, correlation, and peak summaries
+ **Bioconductor**: for R based QC visualizations
+ **Python**: for duplicate capping logic, metric generation, and workflow support scripts

**Note**: All tool paths, versions, and parameters are centrally managed in `config.yml`, making the pipeline reproducible, portable, and easy to maintain.

+ The workflow was developed and tested using **Snakemake 7.32.4**
+ The included `run_workflow.sbatch` script uses Snakemake's generic `--cluster` and `--cluster-config` submission interface
+ Other Snakemake versions, particularly version 8 or newer, may require migration to a Slurm executor plugin

---

## 5) Example Data

A compact dataset is included within the repository for testing purposes, along with example scripts for analyzing publicly available CUT&Tag datasets. This pipeline extends the original protocol and provides a reproducible framework for routine alignment, QC, duplicate capping, random downsampling, signal generation, and sample level assessment.

---

## 6) Explanation of `samples.csv`

**IMPORTANT:** Always update `config/samples.csv` with your sample names, FASTQ file paths, histone or target labels, and replicate information before running.

The `samples.csv` must contain the following columns:

| sample               | fastq1                                   | fastq2                                   | histone            | replicate    |
|----------------------|------------------------------------------|------------------------------------------|--------------------|--------------|
| H3K27ac_1_control    | path/to/H3K27ac_1_control_R1.fastq.gz    | path/to/H3K27ac_1_control_R2.fastq.gz    | H3K27ac_Control    | 1            |
| H3K27ac_2_control    | path/to/H3K27ac_2_control_R1.fastq.gz    | path/to/H3K27ac_2_control_R2.fastq.gz    | H3K27ac_Control    | 2            |
| H3K27me3_1_treatment | path/to/H3K27me3_1_treatment_R1.fastq.gz | path/to/H3K27me3_1_treatment_R2.fastq.gz | H3K27me3_Treatment | 1            |
| H3K27me3_2_treatment | path/to/H3K27me3_2_treatment_R1.fastq.gz | path/to/H3K27me3_2_treatment_R2.fastq.gz | H3K27me3_Treatment | 2            |

+ **sample**: Unique name used for all output files
+ **fastq1 / fastq2**: Paths to R1 and R2 FASTQ files
+ **histone**: Mark, target, or sample group name used in plots
+ **replicate**: Replicate number or label, such as `1`, `2`, `A`, or `B`

Include both the histone or target and the sample type in the `histone` column when possible. This determines how samples are grouped in plots and can improve interpretation of replicate level QC.

Use informative `sample` names that match your design and avoid periods in FASTQ/sample names before the file extension.

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

| Category                           | Output Folder/Format                                                                  |
|------------------------------------|---------------------------------------------------------------------------------------|
| **Optional Raw FASTQ QC**          | `results/qc/fastqc/`, `results/qc/multiqc/`, `results/qc/fastq_screen/`               |
| **Optional Trimmed FASTQs**        | `results/trimming/`, with reports in `results/qc/cutadapt/`                           |
| **Alignments**                     | `results/alignment/sam/`, `results/alignment/bam/`                                    |
| **Duplicate Capping Metrics**      | `results/qc/duplicates/`                                                              |
| **Downsampling Metrics**           | `results/qc/downsampling/`                                                            |
| **Duplicate/Downsampling Summary** | `results/qc/duplicates_downsampling/`, `results/plots/duplicate_downsampling_summary.png` |
| **Final Analysis BAMs**            | `results/alignment/bam/analysis/` when downsampling is enabled; otherwise selected final BAM path depends on config |
| **BigWig Tracks**                  | `results/alignment/bigwig/`                                                           |
| **Peak Calls**                     | `results/peakCalling/` `.narrowPeak` files                                            |
| **Fragment Length Files**          | Raw and final fragment length text files                                              |
| **FRiP Score & Coverage**          | FRiP score tables and binned fragment count files                                     |
| **Plots**                          | `results/plots/` for alignment, duplicate/downsampling, fragment length, correlation, peak count, and FRiP summaries |

### **BigWig Track Types**

This pipeline generates three types of BigWig tracks for genome browser visualization:

+ **Analysis BigWig**
  + Direct coverage from the final analysis BAM which could be
    + aligned sorted raw BAM
    + aligned sorted BAM after downsampling
    + duplicate capped BAM
    + duplicate capped BAM after downsampling
  + Output pattern: `results/alignment/bigwig/{sample}_analysis.bw`

+ **CPM normalized BigWig**
  + Generated using `bamCoverage --normalizeUsing CPM`
  + Output pattern: `results/alignment/bigwig/{sample}_analysisCPM.bw`
  + Useful for comparing signal after normalization to counts per million mapped reads

+ **Target scaled BigWig**
  + Generated using `bamCoverage --scaleFactor`
  + Output pattern: `results/alignment/bigwig/{sample}_analysisScaled.bw`
  + Uses a scale factor calculated as:

```text
scale_factor = target_fragments / final_analysis_fragments
```

The meaning of `target_fragments` depends on the workflow settings:

+ If downsampling is enabled:
  + `target_fragments` is the resolved downsampling target

+ If downsampling is disabled:
  + `target_fragments` is the empirical lowest final analysis fragment count and is used only as a BigWig scaling reference

**Important:** Target scaled BigWigs are for visualization and exploratory comparison. The workflow does not upsample BAM files or duplicate reads to artificially increase sequencing depth.

### **Example BigWig Scaling Behavior**

Example final analysis fragment depths:

```text
Sample 1 = 30,000,000
Sample 2 = 28,000,000
Sample 3 = 18,000,000
Sample 4 = 24,000,000
Sample 5 =  8,000,000
```

If downsampling is disabled, no BAM downsampling target is used. Instead, the BigWig scaling reference is resolved to the empirical lowest final analysis fragment count:

```text
Scaling reference = 8,000,000
```

The target scaled BigWig scale factors are:
```text
Sample 1: 8,000,000 / 30,000,000 = 0.267
Sample 2: 8,000,000 / 28,000,000 = 0.286
Sample 3: 8,000,000 / 18,000,000 = 0.444
Sample 4: 8,000,000 / 24,000,000 = 0.333
Sample 5: 8,000,000 /  8,000,000 = 1.000
```

Sample 5 is unchanged, and the higher depth samples are scaled down in the `analysisScaled.bw` files.

If using `lowest_with_floor` with the same example depths:

```text
Sample 1 = 30,000,000
Sample 2 = 28,000,000
Sample 3 = 18,000,000
Sample 4 = 24,000,000
Sample 5 =  8,000,000
```

and the resolved target is:
```text
18,000,000
```

The target-scaled BigWig scale factors are:
```text
Sample 1: 18,000,000 / 30,000,000 = 0.600
Sample 2: 18,000,000 / 28,000,000 = 0.643
Sample 3: 18,000,000 / 18,000,000 = 1.000
Sample 4: 18,000,000 / 24,000,000 = 0.750
Sample 5: 18,000,000 /  8,000,000 = 2.250
```

This scaling affects only the BigWig visualization tracks, not the BAM files.

---

## 8) Example Output Plots

Below are example plots generated by this pipeline.

| 1. **Alignment Summary Plot**                                                         |
| :-----------------------------------------------------------------------------------: |
| <img src="/images/alignment_summary_plot.png" width="800">                            |
| *Summary of alignment stats*                                                          |

| 2. **Fragment Length Plot**                                                           | 3. **Fragment Count Correlation Plot**                                                |
| :-----------------------------------------------------------------------------------: | :-----------------------------------------------------------------------------------: |
| <img src="/images/fragment_length_plot.png" width="600">                              | <img src="/images/fragCount_correlation_plot.png" width="600">                        |
| *Final analysis fragment length summary plots*                                        | *Correlation of fragment counts across shared occupied 500-bp bins*                   |

| 4. **Duplicate and Downsampling Summary Plot**                                        | 5. **FRiP and Peak Count Summary Plot**                                               |
| :-----------------------------------------------------------------------------------: | :-----------------------------------------------------------------------------------: |
| <img src="/images/duplicate_downsampling_summary.png" width="600">                    | <img src="/images/peak_summary_plot.png" width="600">                                 |
| *Duplicate burden and retained fragment summary*                                      | *FRiP and called peak summary*                                                        |

---

## 9) Testing Several Duplicate Caps

To compare the effect of duplicate capping, run the workflow in separate directories using different values of `duplicate_cap_max`.

Example runs:

```text
capOff_run
cap1_run
cap3_run
cap5_run
cap7_run
```

Example configs:

```yaml
use_duplicate_cap: false
duplicate_cap_max: 5       # <- This value doesn't matter when set to false
```

```yaml
use_duplicate_cap: true
duplicate_cap_max: 1
```

```yaml
use_duplicate_cap: true
duplicate_cap_max: 3
```

```yaml
use_duplicate_cap: true
duplicate_cap_max: 5
```

```yaml
use_duplicate_cap: true
duplicate_cap_max: 7
```

Duplicate capping is useful as a sensitivity analysis because duplicate fragments can reflect both technical amplification and true biological enrichment. Very strict caps, such as `duplicate_cap_max: 1`, reduce the influence of PCR duplicates but may remove real signal from highly enriched CUT&Tag sites. More permissive caps, such as `duplicate_cap_max: 5` or `duplicate_cap_max: 7`, retain more fragment complexity and signal intensity but may allow technical duplication to have a stronger effect.

Comparing several duplicate caps can help evaluate how duplicate burden affects retained fragments, signal tracks, FRiP scores, peak counts, and replicate correlations. A reasonable cap should reduce excessive duplicate-driven signal while preserving expected enrichment patterns, replicate consistency, and interpretable peak calls.

*Note: The optimal duplicate cap may depend on library complexity, sequencing depth, antibody/target behavior, and the expected signal distribution for the mark being profiled. This workflow does not assume a single universal duplicate cap for all datasets.*

---

## 10) Important Interpretation Notes

Throughout duplicate metrics, downsampling calculations, scale-factor calculations, and fragment-length summaries, one mapped fragment is represented by the primary read-1 alignment from a properly paired template. Unmapped, secondary, and supplementary alignments are excluded.

### **Duplicate Capping**

Duplicate capping limits the number of identical fragments retained at each genomic fragment position. It is intended to reduce the influence of PCR over amplification while avoiding complete duplicate removal.

+ `duplicate_cap_max: 1` is strict and behaves similarly to keeping only one fragment per exact position
+ `duplicate_cap_max: 3` is intermediate
+ `duplicate_cap_max: 5` is more permissive
+ `use_duplicate_cap: false` leaves duplicate fragments unchanged

The best cap may depend on sample quality, sequencing depth, mark type, and the amount of PCR duplication.

---

### **Downsampling**

Downsampling removes observed fragments from samples above the resolved target. It does not create artificial reads.

+ Samples above the target are downsampled
+ Samples at or below the target are left unchanged
+ Downsampling uses an integer as the configured seed to produce reproducible selection.
+ Another option in this workflow is using `"random"` which generates a new seed each time the workflow is run and downsampling is performed
+ Downsampling affects the final analysis BAM when enabled

Downsampling can help make samples more comparable when sequencing depth differs substantially, but aggressive downsampling can also remove useful signal.

---

### **BigWig Scaling**

This workflow produces three BigWig types:

+ Raw (possibly capped and/or downsampled) coverage
+ CPM normalized coverage
+ Target scaled coverage

The target scaled BigWig is for browser visualization and exploratory comparisons. It does not change the BAM file and does not rescue low depth samples for downstream statistical analysis.

When downsampling is disabled, the target scaled BigWig uses the empirical lowest final analysis fragment count as the scaling reference.

When downsampling is enabled, the target scaled BigWig uses the resolved downsampling target as the scaling reference.

---

### **Per Sample Peak Calling**

This workflow performs per sample peak calling with MACS2. These peaks are useful for sample level QC and FRiP score calculation.

This workflow does not define final consensus peak sets or merge biological replicates. Those steps are handled by the downstream companion workflow.

FRiP is calculated using the cleaned fragment BED generated from the final analysis BAM. This includes same-chromosome fragment pairs shorter than 1,000 bp that pass the workflow's BEDPE fragment-processing filters.

---

### **Fragment Count Correlation**

The correlation plot is calculated from fragment midpoint counts in 500-bp bins generated from cleaned final analysis fragments. The binned fragment count files contain only bins with at least one fragment, and missing bins are not explicitly represented as zero-count bins.

When the sample files are combined, correlation is calculated using complete observations. Therefore, only 500-bp bins containing fragments in every sample included in the correlation matrix contribute to the reported correlation values. Bins occupied in only a subset of samples and genome-wide bins with no fragments are excluded.

Restricting the calculation to shared occupied bins conditions the analysis on regions containing signal across all samples and may inflate correlation values relative to methods that evaluate a broader set of genomic bins. The correlation plot should therefore be interpreted as an exploratory sample-level QC measure of similarity among shared signal-containing regions, rather than as an unbiased genome-wide estimate of concordance or a formal measure of replicate reproducibility.

---

## 11) Instructions to run on Slurm managed HPC

### 11A. Download version controlled repository

```bash
wget https://github.com/KirklandLab/CutandTag_Alignment_QC/releases/download/v2.0.0/CutandTag_Alignment_QC-2.0.0.tar.gz
tar -xzf CutandTag_Alignment_QC-2.0.0.tar.gz
rm CutandTag_Alignment_QC-2.0.0.tar.gz
cd CutandTag_Alignment_QC-2.0.0
```

### 11B. Load modules

```bash
module purge
module load slurm python/3.10 pandas/2.2.3 numpy/1.22.3 matplotlib/3.7.1
```

### 11C. Modify samples and config files

```bash
vim config/samples.csv
vim config/config.yml
```

### 11D. Dry run

```bash
snakemake -npr
```

### 11E. Run on HPC with config options

Submit the full workflow using the included Slurm submission script:
```bash
sbatch run_workflow.sbatch
```

### 11F. Unlock a workflow directory if needed

If a previous Snakemake run was interrupted, unlock the directory before rerunning:
```bash
snakemake --unlock
```

### 11G. Run only a specific target and its upstream rules

To run the workflow to a specific file and its required upstream rules, provide the path and filename after the script:
```bash
sbatch run_workflow.sbatch results/plots/alignment_summary_plot.png
```
Also the rule itself can be targeted by name:
```bash
sbatch run_workflow.sbatch generate_alignment_plots
```

*Note: When no target is provided, Snakemake runs the default `rule all` and completes the full workflow.*

---

## 12) Citation

If you use this workflow in your research, please cite both the original CUT&Tag protocol and this pipeline:

+ Zheng, Y., Ahmad, K., & Henikoff, S. (2019). CUT&Tag for efficient epigenomic profiling of small samples and single cells. *Protocols.io*. https://dx.doi.org/10.17504/protocols.io.bjk2kkye

+ **Boyd, K.A.** (2025). *CutandTag_Alignment_QC: A reproducible Snakemake workflow for alignment, quality control, and signal generation from Cut-and-Tag sequencing data*. *Zenodo*. https://doi.org/10.5281/zenodo.15232228

[![DOI](https://zenodo.org/badge/873121124.svg)](https://doi.org/10.5281/zenodo.15232228)

---

## 13) Authorship & Contributions

**Kevin A. Boyd** – Designed and implemented the Snakemake workflow for a Slurm managed HPC environment, modularized the pipeline structure, implemented the processing steps, integrated duplicate capping, downsampling, signal normalization, target scaled BigWig support, quality control plots, and documentation.

**Jacob Kirkland** – Principal Investigator; provided project direction, conceptual guidance, and experimental data for workflow development and validation.

This work was developed by Kevin A. Boyd for the Kirkland Lab, with scientific direction and resources provided by Jacob Kirkland, as part of a COBRE funded collaborative effort at OMRF (Oklahoma Medical Research Foundation). While the pipeline was initially developed for use within the Kirkland Lab, it is broadly applicable to CUT&Tag data analysis in other research settings.

---

## 14) License

This project is licensed under the **Apache 2.0**. See the [LICENSE](LICENSE) file for details.

[![License: Apache 2.0](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

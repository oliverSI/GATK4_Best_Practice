# GATK4 Best Practice Nextflow Pipeline
**GATK4 Best Practice Nextflow Pipeline** is a variant calling pipeline for human whole genome sequencing based on the GATK best practice.

## Requirements
* <a href="https://www.nextflow.io/">NextFlow</a>
* <a href="https://www.docker.com/">Docker</a>


## Running the pipeline
There is no need to download the code explicitly. Nextflow seamlessly intergrates with GitHub with the following command.
  ```
  nextflow run oliverSI/GATK4_Best_Practice --fastq1 read_R1.fastq.gz --fastq2 read_R2.fastq.gz
  ```
  ```
  N E X T F L O W  ~  version 0.25.6
Launching `main.nf` [cheesy_roentgen] - revision: 2ba478ba4b
====================================================================
GATK4 Best Practice Nextflow Pipeline (v0.1)                        
====================================================================
[warm up] executor > local
[89/1d86be] Cached process > get_phase1_SNPs
[74/369556] Cached process > get_reference
[d8/8ace49] Cached process > get_omni
[6b/a2f67e] Cached process > get_dbSNP
[14/831c24] Cached process > get_golden_indel
[ea/1836b3] Cached process > get_hapmap
[6a/b862c1] Cached process > get_BWA_index
[ad/adb689] Cached process > BWA
[15/7d3a72] Cached process > BWA_sort
[9a/56cc98] Cached process > MarkDuplicates
[11/7b55c6] Cached process > BaseRecalibrator
[e0/984db2] Cached process > ApplyBQSR
[43/c70b8e] Cached process > HaplotypeCaller
[c2/30fac6] Cached process > GenotypeGVCFs
[84/74b170] Cached process > VariantRecalibrator_SNPs
[31/37e241] Cached process > ApplyVQSR_SNPs
[9e/a29311] Cached process > VariantRecalibrator_INDELs
[29/98dfcd] Cached process > ApplyVQSR_INDELs
[3b/7de781] Submitted process > copy
```
## Pipeline parameters
```
../nextflow-0.25.5/nextflow run main.nf --help
```
```
N E X T F L O W  ~  version 0.25.5
Launching `main.nf` [silly_baekeland] - revision: 82d1c9f7ca
====================================================================
GATK4 Best Practice Nextflow Pipeline (v0.1)                        
====================================================================
 
USAGE: 
 
nextflow run oliverSI/GATK4_Best_Practice --fastq1 read_R1.fastq.gz --fastq2 read_R2.fastq.gz
 
Mandatory arguments:
    --fastq1        FILE               Fastq(.gz) file for read1
    --fastq2        FILE               Fastq(.gz) file for read2
 
Optional arguments:
    --outdir        DIR                Output directory(default: ./Results)
    --samplename    STRING             Sample name(dafault: fastq1 basename)
    --rg            STRING             Read group tag(dafault: fastq1 basename)
 
====================================================================
```

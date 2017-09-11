# GATK4 Best Practice Nextflow Pipeline
**GATK4 Best Practice Nextflow Pipeline** is a bioinformatics pipeline for variant calling based on the GATK best practice.

## Requirements
* <a href="https://www.nextflow.io/">NextFlow</a>
* <a href="https://www.docker.com/">Docker</a>


## Running the pipeline
There is no need to download the code explicitly. Nextflow seamlessly intergrates with GitHub with the following command.
  ```
  $ nextflow run oliverSI/GATK4_Best_Practice --fastq1 read_R1.fastq.gz --fastq2 read_R2.fastq.gz
  ```

## Pipeline parameters
```
$ ../nextflow-0.25.5/nextflow run main.nf --help
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

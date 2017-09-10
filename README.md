# GATK4 Best Practice Nextflow Pipeline
**GATK4 Best Practice Nextflow Pipeline** is a bioinformatics pipeline for variant calling based on the GATK best practice.

## Requirements
* <a href="https://www.nextflow.io/">NextFlow</a>
* <a href="https://www.docker.com/">Docker</a>


## Running the pipeline
There is no need to download the code explicitly. Nextflow seamlessly intergrates with GitHub with the following command.
  ```
  nextflow run oliverSI/GATK4_Best_Practice --fastq1 read_R1.fastq.gz --fastq2 read_R2.fastq.gz
  ```

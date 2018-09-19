# NGS_SVs_Pipelines

Snakemake pipelines for NGS Structural Variants (SVs) and Copy Numbers Variation (CNVs) detections, annotation and visualization. The software package includes two pipelines, one pipeline is for handling illumina short-read data and the other one is 
for analyzing long-read sequencign data from PacBio and 10x Genomics technologies.

The pipeline features include NGS data preprocessing and QC, SVs detections, CNV calling, consensus SVs/CNVs calling from different software tools. The pipeline also includes a set of tools for structural variants annotation and visualization in order to help the result interpration.


## Prerequisites

 * Software packages
 
   Install all the programs listed in program.py 
   
 * Reference genomes
   Install all the reference files listed in the reference.py
   

## Software Installation

   download and install all Prerequisites and Referencee genome files according to the specific software tool's installation instruction.
   download and install the snakemake file from this software package.

## Running the NGS SVs Pipelines


    > run_submit.sh config.py

## Contact

  CCRSF_IFX@nih.gov



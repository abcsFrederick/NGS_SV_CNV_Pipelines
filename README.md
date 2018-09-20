# NGS_SV_CNV_Pipelines

Snakemake pipelines for NGS Structural Variants (SVs) and Copy Numbers Variation (CNVs) detections, annotation and visualization. The software package includes two pipelines, one pipeline is for handling Illumina short-read data and the other one is 
for analyzing long-read sequencing data from PacBio and 10x Genomics platforms.

The pipeline includes features such as NGS data preprocessing and QC, SVs detections, CNV calling, consensus SVs calling from different software tools. The pipeline also includes a set of tools for structural variants annotation and visualization in order to help the result interpration.


## Prerequisites

 * Software packages
 
   Install all the programs listed in program.py 
   
 * Reference genomes
   Install all the reference files listed in the reference.py
   

## Software Installation

   Please download and install all prerequisites and referencee genome files according to the specific software tool's installation instruction. Please download and install the snakemake file from this software package.

## Running the NGS SV Pipelines


    > submit.sh SNAKEMAKE

## Contact

  CCRSF_IFX@nih.gov



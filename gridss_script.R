library(stringr)
library(VariantAnnotation)
library(StructuralVariantAnnotation)
library(tools)

args <- commandArgs(TRUE)
show(args)
show(str(args))
inputfile <- args[1]

extension <- file_ext(inputfile)
newfile <- paste(substring(inputfile, 1, nchar(inputfile)-nchar(extension)-1), 'tumor.specific', extension, sep=".")

vcf <- readVcf(inputfile, "hg38")
vcf <- vcf[rowRanges(vcf)$FILTER %in% c(".", "PASS"),]
somatic_vcf <- vcf[geno(vcf)$QUAL[,2] == 0,]

writeVcf(somatic_vcf,file=newfile)

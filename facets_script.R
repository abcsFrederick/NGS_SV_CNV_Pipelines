library(facets)

args <- commandArgs(TRUE)

#show(args)

inputFile = args[1]
#prefix = args[2]
#outputPath = args[3]

datafile = (inputFile)
rcmat = readSnpMatrix(datafile)
xx = preProcSample(rcmat,ndepth=15, het.thresh=0.25, snp.nbhd=1000, cval=75,deltaCN=0, gbuild=c("hg38"),hetscale=TRUE, unmatched=FALSE, ndepthmax=1000)
oo = procSample(xx, cval=300, min.nhet=15)
fit=emcncf(oo)

pdf("facets.pdf")
plotSample(x=oo,emfit=fit)
dev.off()

write.table(fit$cncf,"facets.cnv.tsv",sep="\t",row.names=FALSE)
write.table(fit$ploidy,"facets.poidy.tsv",sep="\t",row.names=FALSE)
write.table(fit$purity,"facets.purity.tsv",sep="\t",row.names=FALSE)

pdf("facets.model.pdf")
logRlogORspider(oo$out, oo$dipLogR)
dev.off()
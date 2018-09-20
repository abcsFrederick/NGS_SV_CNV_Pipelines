library(sequenza)

args <- commandArgs(TRUE)

show(args)

inputFile = args[1]
#prefix = args[2]
#outputPath = args[3]

data.file = (inputFile)
seqzdata = sequenza.extract(data.file, min.reads = 15, min.reads.normal= 10)
CP.example = sequenza.fit(seqzdata)
sequenza.results(sequenza.extract = seqzdata,cp.table = CP.example,sample.id = "seqz")
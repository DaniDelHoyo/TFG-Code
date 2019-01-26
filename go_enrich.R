#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#Go enrichment

hitInSample=as.integer(args[[1]])
hitInPop = as.integer(args[[2]])
totalSize=as.integer(args[[3]])
sampleSize = as.integer(args[[4]])
failInPop = totalSize - hitInPop 


fisher.test(matrix(c(hitInSample, hitInPop-hitInSample,sampleSize-hitInSample,
                     failInPop-sampleSize +hitInSample), 2, 2), 
            alternative='greater')$p.value
